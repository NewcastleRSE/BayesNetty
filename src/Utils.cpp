/************************************************************************
 * BayesNetty, version 1.1
 * Copyright 2015-present,
 * Richard Howey
 * Institute of Genetic Medicine, Newcastle University
 *
 * richard.howey@ncl.ac.uk
 * http://www.staff.ncl.ac.uk/richard.howey/
 *
 * This file is part of BayesNetty, the SNP interaction analysis program.
 *
 * BayesNetty is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BayesNetty is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with BayesNetty.  If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/


/*! \file Utils.cpp
    \brief This file contains the source some useful functions for calculating statistics. .
    
*/

#include <map>
#include <list>
#include <math.h>
#include <string>
#include <sstream>

using namespace std; // initiates the "std" or "standard" namespace

#include "Data.h"
#include "Utils.h"
#include "cdflib.h"
#include "main.h"

//! Converts an integer to a string.
string toString(const unsigned int & i)
{
	ostringstream aStringStream;
	aStringStream << i;

	return aStringStream.str();
};

//! Converts an integer to a string.
string toString2DP(const double & d, unsigned int minLength)
{
	ostringstream aStringStream;
	aStringStream << std::setprecision(4);
	aStringStream << d;

	string ans = aStringStream.str();
	while(ans.length() < minLength) ans += " ";


	return ans;
};

//! Inverts a matrix. the inverse should be set to the identity when given to this function.
void getInverseMatrix(list< list<double> > & matrix, list< list<double> > & inverse)
{
	if(matrix.size() == 0) return;

	//loop thro' rows of matrix, m, and the inverse, i
	list< list<double> >::iterator mrow = matrix.begin();
	list< list<double> >::iterator irow = inverse.begin();
	list< list<double> >::iterator mrow2;
	list< list<double> >::iterator irow2;
	list<double>::iterator mcol;
	list<double>::iterator icol;
	list<double>::iterator mcol2;
	list<double>::iterator icol2;
	double factor;
	unsigned int rowNo = 1;
	unsigned int colNo = 1;

	for( ; mrow != matrix.end(); ++mrow, ++irow)
	{
		//set column to the first column in the row
		mcol = mrow->begin();
		icol = irow->begin();
		colNo = 1;

		//advance the column until the the row no. is equal to the column no.
		while(colNo != rowNo)
		{
			mcol++;			
			++colNo;
		};

		//divide the row in (m and i) by the value in the matrix, m, at (rowNo, colNo)
		factor = 1.0/(*mcol); //divide all elements instead?
		*mcol = 1;		
		mcol++;		

		//scale the remaining elements in the row - if there are any
		while(mcol != mrow->end())
		{
			*mcol *= factor;			
			mcol++;			
		};

		//scale all of the elements in the inverse for this row
		while(icol != irow->end())
		{			
			*icol *= factor;			
			icol++;
		};

		//subtract the row in question away from the remaining rows scaled  s.t. column = mrow will be zero below this row in matrix m
		mrow2 = mrow;
		irow2 = irow;
		mrow2++;
		irow2++;

		//loop thro' remaining rows
		while(mrow2 != matrix.end())
		{
			//set column iterators to begining of the rows
			mcol2 = mrow2->begin();
			icol2 = irow2->begin();
			mcol = mrow->begin();
			icol = irow->begin();

			//advance column of matrix, m, to the same as the main row being considered, with value rowNo
			colNo = 1;
			while(colNo != rowNo)
			{
				mcol++;
				mcol2++;				
				++colNo;
			};

			factor = *mcol2; //factor to multiple row, rowNo by to take away from subseq. rows
			*mcol2 -= (*mcol)*factor;//0;
			mcol++;
			mcol2++;

			//subtract scaled row for the rest of the matrix, m
			while(mcol2 != mrow2->end())
			{
				*mcol2 -= (*mcol)*factor;				
				mcol++;
				mcol2++;				
			};

			//now perform the same row operation on the inverse matrix, i
			while(icol2 != irow2->end())
			{
				*icol2 -= (*icol)*factor;
				icol++;
				icol2++;				
			};

			mrow2++;
			irow2++;
		};//end of performing row operations to set column (=rowNo) to zeroes below row = rowNo

		++rowNo;
	};//end of performing row operations to set lower left of matrix, m, to zeroes

	//Now reduce the upper right of matrix, m, to zero
	list< list<double> >::reverse_iterator mrowre = matrix.rbegin();
	list< list<double> >::reverse_iterator irowre = inverse.rbegin();
	list< list<double> >::reverse_iterator mrowre2 = matrix.rbegin();
	list< list<double> >::reverse_iterator irowre2 = inverse.rbegin();
	list<double>::reverse_iterator mcolre2;
	list<double>::reverse_iterator mcolre;

	rowNo = matrix.size();

	for( ; mrowre != matrix.rend(); ++mrowre, ++irowre)
	{

		mrowre2 = mrowre;
		irowre2 = irowre;
		mrowre2++;
		irowre2++;

		//loop tho' the remaining rows backwards - if there are any
		while(mrowre2 != matrix.rend())
		{			
			//set column iterators to begining of the rows
			mcolre2 = mrowre2->rbegin();
			icol2 = irowre2->begin();
			mcolre = mrowre->rbegin();
			icol = irowre->begin();

			//advance column of matrix, m, to the same as the main row being considered, with value rowNo
			colNo = mrowre2->size();//size will be 4
			while(colNo != rowNo)
			{
				mcolre++;
				mcolre2++;				
				--colNo;
			};

			factor = *mcolre2; //factor to multiple row, rowNo by to take away from subseq. rows
			*mcolre2 -= (*mcolre)*factor;
			mcolre++;
			mcolre2++;

			//subtract scaled row from the rest of the matrix, m
			while(mcolre2 != mrowre2->rend())//could stop at when col < row
			{
				*mcolre2 -= (*mcolre)*factor;				
				mcolre++;
				mcolre2++;				
			};

			//now perform the same row operation on the inverse matrix, i
			while(icol2 != irowre2->end())
			{
				*icol2 -= (*icol)*factor;
				icol++;
				icol2++;				
			};			

			mrowre2++;
			irowre2++;
		};

		--rowNo;
	};

};

//! Multiplies two matrices together, matrixProd should be empty when given to this function.
void getMatrixMulti(const list< list<double> > & matrix1, const list< list<double> > & matrix2, list< list<double> > & matrixProd)
{
	if(freeClearMemory) list< list<double> >().swap(matrixProd); else matrixProd.clear();

	if(matrix1.size() == 0) return;

	list< list<double>::const_iterator > mat2ColIters; //iterators for the matrix 2, one for each row

	for(list< list<double> >::const_iterator m2row = matrix2.begin(); m2row != matrix2.end(); ++m2row)
	{
		mat2ColIters.push_back((*m2row).begin()); 
	};

	list< list<double>::const_iterator >::iterator m2rowIt;

	double aValue;
	list<double> aProdRow;
	unsigned int prodMatNoCols = matrix2.begin()->size();

	//loop thro' rows
	for(list< list<double> >::const_iterator m1row = matrix1.begin(); m1row != matrix1.end(); ++m1row)
	{
		aProdRow.clear();
		for(unsigned int prodCol = 1; prodCol <= prodMatNoCols; ++prodCol)
		{
			aValue = 0;
		
			//loop thro coloumns of matrix1 and rows of matrix 2
			m2rowIt = mat2ColIters.begin();
			for(list<double>::const_iterator m1col = m1row->begin(); m1col != m1row->end(); ++m1col, ++m2rowIt)
			{
				//loop thro rows of matrix2			
				aValue += (*m1col)*(*(*m2rowIt));

				++(*m2rowIt); //move on for next time
			};

			aProdRow.push_back(aValue);
		};
		
		matrixProd.push_back(aProdRow);

		//reset columns of matrix2 to beginning
		m2rowIt = mat2ColIters.begin(); 
		for(list< list<double> >::const_iterator m2row = matrix2.begin(); m2row != matrix2.end(); ++m2row, ++m2rowIt)
		{
			*m2rowIt = (*m2row).begin(); //move on to next column for next time
		};

	}; //end of matrix 1 row

};

//! Multiplies a matrix with a vector.
void getMatrixVecMulti(const list< list<double> > & matrix, const list<double> & vec, list<double> & vecAns)
{
	if(freeClearMemory) list<double>().swap(vecAns); else vecAns.clear();
	if(vec.size() == 0) return;

	list<double>::const_iterator vecIt;
	double aValue;

	for(list< list<double> >::const_iterator m1row = matrix.begin(); m1row != matrix.end(); ++m1row)
	{
		aValue = 0;
		vecIt = vec.begin();
		for(list<double>::const_iterator m1col = m1row->begin(); m1col != m1row->end(); ++m1col, ++vecIt)
		{
			aValue += (*m1col)*(*vecIt);
		};
		vecAns.push_back(aValue);
	};
};

//! Multiplies a vector with a matrix.
void getVecMatrixMulti(const list<double> & vec, const list< list<double> > & matrix, list<double> & vecAns)
{
	if(freeClearMemory) list<double>().swap(vecAns); else vecAns.clear();

	if(vec.size() == 0) return;
	

	list< list<double>::const_iterator > mat2ColIters; //iterators for the matrix 2, one for each row

	for(list< list<double> >::const_iterator m2row = matrix.begin(); m2row != matrix.end(); ++m2row)
	{
		mat2ColIters.push_back((*m2row).begin()); 
	};

	list< list<double>::const_iterator >::iterator m2rowIt;

	double aValue;
	unsigned int prodMatNoCols = matrix.begin()->size();
	
	if(freeClearMemory) list<double>().swap(vecAns); else vecAns.clear();
	for(unsigned int prodCol = 1; prodCol <= prodMatNoCols; ++prodCol)
	{
		aValue = 0;
		
		//loop thro coloumns of vec and rows of matrix
		m2rowIt = mat2ColIters.begin();
		for(list<double>::const_iterator m1col = vec.begin(); m1col != vec.end(); ++m1col, ++m2rowIt)
		{
			//loop thro rows of matrix2			
			aValue += (*m1col)*(*(*m2rowIt));

			++(*m2rowIt); //move on for next time
		};

		vecAns.push_back(aValue);
	};
		
};

//! Multiplies a matrix transpose with itself M^t M, give M^T to the function.
void getMatrixTransMatrixMulti(const list< list<double> > & matrixT, list< list<double> > & matrixAns)
{
	
	if(freeClearMemory) list< list<double> >().swap(matrixAns); else matrixAns.clear();
	if(matrixT.size() == 0) return;

	list<double>::const_iterator col1, row2;
	double aValue;

	list<double> aRowAns;

	for(list< list<double> >::const_iterator row1 = matrixT.begin(); row1 != matrixT.end(); ++row1)
	{
		
		aRowAns.clear();

		//loop thro columns of matrix
		for(list< list<double> >::const_iterator col2 = matrixT.begin(); col2 != matrixT.end(); ++col2)
		{
			aValue = 0;
			row2 = col2->begin();
			for(col1 = row1->begin(); col1 != row1->end(); ++col1, ++row2)
			{
				aValue += (*col1)*(*row2);
			};
			aRowAns.push_back(aValue);
		};
		matrixAns.push_back(aRowAns);
	};

};

//! Multiplies a matrix transpose with another matrix.
void getMatrixTransMatrixMulti(const list< list<double> > & matrix1, const list< list<double> > & matrix2, list< list<double> > & matrixProd)
{	
	if(freeClearMemory) list< list<double> >().swap(matrixProd); else matrixProd.clear();
	if(matrix1.size() == 0) return;

	list< list<double>::const_iterator > mat1ColIters; //iterators for the matrix 1, one for each row
	for(list< list<double> >::const_iterator m1row = matrix1.begin(); m1row != matrix1.end(); ++m1row)
	{
		mat1ColIters.push_back((*m1row).begin()); 
	};
	list< list<double>::const_iterator >::const_iterator m1rowIt;


	list< list<double>::const_iterator > mat2ColIters; //iterators for the matrix 2, one for each row
	for(list< list<double> >::const_iterator m2row = matrix2.begin(); m2row != matrix2.end(); ++m2row)
	{
		mat2ColIters.push_back((*m2row).begin()); 
	};
	list< list<double>::const_iterator >::iterator m2rowIt;

	double aValue;
	list<double> aProdRow;
	unsigned int prodMatNoRows = matrix1.begin()->size(); //no rows in t(mat1) is no of columns in mat1
	unsigned int prodMatNoCols = matrix1.size(); //no of columns is no of rows in matrix1 or no of rows of matrix 2 

	//for(list< list<double> >::const_iterator m1row = matrix1.begin(); m1row != matrix1.end(); ++m1row)
	for(unsigned int prodRow = 1; prodRow <= prodMatNoRows; ++prodRow) //loop thro rows
	{
		aProdRow.clear();

		//calculate one value for the matrix
		for(unsigned int prodCol = 1; prodCol <= prodMatNoCols; ++prodCol) //loop thro cols
		{
			aValue = 0;
		
			//loop thro rows of matrix1 and rows of matrix 2
			
			m1rowIt = mat1ColIters.begin();
			m2rowIt = mat2ColIters.begin();
			
			for(unsigned int prodCol2 = 1; prodCol2 <= prodMatNoCols; ++prodCol2, ++m1rowIt, ++m2rowIt)
			{
				//loop thro rows of matrix2			
				aValue += (*(*m1rowIt))*(*(*m2rowIt));

				++(*m2rowIt); //move on to next col for next time
			
			};

			aProdRow.push_back(aValue);
		};
		
		matrixProd.push_back(aProdRow);

		//advance to next row
		for(list< list<double>::const_iterator >::iterator m1r = mat1ColIters.begin(); m1r != mat1ColIters.end(); ++m1r)
		{
			++(*m1r); //move on to next column for next time
		};

		//reset columns of matrix2 to beginning
		m2rowIt = mat2ColIters.begin(); 
		for(list< list<double> >::const_iterator m2row = matrix2.begin(); m2row != matrix2.end(); ++m2row, ++m2rowIt)
		{
			*m2rowIt = (*m2row).begin(); //move on to next column for next time
		};

	}; //end of matrix 1 row

};

//! Multiplies a matrix transposed with a vector.
void getMatrixTransVecMulti(const list< list<double> > & matrix, const list<double> & vec, list<double> & vecAns)
{
	if(freeClearMemory) list<double>().swap(vecAns); else vecAns.clear();
	if(vec.empty()) return;

	list< list<double>::const_iterator > matColIters; //iterators for the matrix, one for each row

	for(list< list<double> >::const_iterator mrow = matrix.begin(); mrow != matrix.end(); ++mrow)
	{
		matColIters.push_back((*mrow).begin()); 
	};

	//loop thro' rows of transposed matrix
	list< list<double>::const_iterator >::iterator mrowIt;
	unsigned int prodMatNoCols = matrix.begin()->size();
	double aValue;
	
	
	for(unsigned int prodCol = 1; prodCol <= prodMatNoCols; ++prodCol)
	{
		aValue = 0;
		
		//loop thro columns of vec and rows of matrix
		mrowIt = matColIters.begin();
		for(list<double>::const_iterator rowVec = vec.begin(); rowVec != vec.end(); ++rowVec, ++mrowIt)
		{
			//loop thro rows of matrix2			
			aValue += (*(*mrowIt))*(*rowVec);

			++(*mrowIt); //move on for next time
		};

		vecAns.push_back(aValue);
	};

};

//! Adds two matrix together.
void getMatrixAddMatrix(const list< list<double> > & matrix1, const list< list<double> > & matrix2, list< list<double> > & matrixAns)
{
	
	if(freeClearMemory) list< list<double> >().swap(matrixAns); else matrixAns.clear();
	if(matrix1.empty()) return;

	list<double> aRow;
	list<double>::const_iterator col2;

	list< list<double> >::const_iterator row2 = matrix2.begin();
	for(list< list<double> >::const_iterator row1 = matrix1.begin(); row1 != matrix1.end(); ++row1, ++row2)
	{
		col2 = row2->begin();
		aRow.clear();
		for(list<double>::const_iterator col1 = row1->begin(); col1 != row1->end(); ++col1, ++col2)
		{
			aRow.push_back(*col1 + *col2);
		};

		matrixAns.push_back(aRow);
	};

};

//! Adds two vectors.
void getVecVecAdd(const list<double> & vec1, const list<double> & vec2, list<double> & vecAns)
{
	if(freeClearMemory) list<double>().swap(vecAns); else vecAns.clear();
	if(vec1.empty()) return;

	list<double>::const_iterator v2 = vec2.begin();
	for(list<double>::const_iterator v1 = vec1.begin(); v1 != vec1.end(); ++v1, ++v2)
	{
		vecAns.push_back((*v1+*v2));
	};
};

//! Subtracts two vectors.
void getVecVecSub(const list<double> & vec1, const list<double> & vec2, list<double> & vecAns)
{
	if(freeClearMemory) list<double>().swap(vecAns); else vecAns.clear();
	if(vec1.empty()) return;

	list<double>::const_iterator v2 = vec2.begin();
	for(list<double>::const_iterator v1 = vec1.begin(); v1 != vec1.end(); ++v1, ++v2)
	{
		vecAns.push_back((*v1-*v2));
	};

};

//! Multiplies a two vectors.
void getVecVecMulti(const list<double> & vec1, const list<double> & vec2, double & ans)
{
	ans = 0;
	if(vec1.empty()) return;

	list<double>::const_iterator v2 = vec2.begin();
	for(list<double>::const_iterator v1 = vec1.begin(); v1 != vec1.end(); ++v1, ++v2)
	{
		ans += (*v1)*(*v2);
	};
};

//! Computes determinant of square (perhaps symetric) matrix.
double getDetSquareMatrix(const list< list<double> > & sqMatrix)
{
	unsigned int matSize = sqMatrix.size();
	if(matSize == 0) return 0;
	else if(matSize == 1)
	{
		double num = *sqMatrix.begin()->begin();
		if(num > 0) return num; else return -num;
	};

	if(matSize == 2)
	{
		double a11, a12, a21, a22;
		list< list<double> >::const_iterator r = sqMatrix.begin();
		list<double>::const_iterator c = (*r).begin();
		a11 = (*c);	++c; a12 = (*c);
		++r; c = (*r).begin();
		a21 = (*c);	++c; a22 = (*c);

		return a11*a22 - a12*a21; 
	};


	double ans = 1;
	list< list<double> > matrix = sqMatrix;

	//loop thro' rows of matrix, m, and the inverse, i
	list< list<double> >::iterator mrow = matrix.begin();
	list< list<double> >::iterator mrow2;
	list<double>::iterator mcol;
	list<double>::iterator mcol2;
	double factor;
	unsigned int rowNo = 1;
	unsigned int colNo = 1;

	for( ; mrow != matrix.end(); ++mrow)
	{
		//set column to the first column in the row
		mcol = mrow->begin();		
		colNo = 1;

		//advance the column until the the row no. is equal to the column no.
		while(colNo != rowNo)
		{
			mcol++;			
			++colNo;
		};

		//divide the row in (m and i) by the value in the matrix, m, at (rowNo, colNo)
		factor = 1.0/(*mcol); //divide all elements instead?
		ans *= (*mcol); //calculate determinant
		*mcol = 1;		
		mcol++;		

		//scale the remaining elements in the row - if there are any
		while(mcol != mrow->end())
		{
			*mcol *= factor;			
			mcol++;			
		};

		//subtract the row in question away from the remaining rows scaled  s.t. column = mrow will be zero below this row in matrix m
		mrow2 = mrow;
		mrow2++;

		//loop thro' remaining rows
		while(mrow2 != matrix.end())
		{
			//set column iterators to begining of the rows
			mcol2 = mrow2->begin();		
			mcol = mrow->begin();		

			//advance column of matrix, m, to the same as the main row being considered, with value rowNo
			colNo = 1;
			while(colNo != rowNo)
			{
				mcol++;
				mcol2++;				
				++colNo;
			};

			factor = *mcol2; //factor to multiple row, rowNo by to take away from subseq. rows
			*mcol2 -= (*mcol)*factor;//0;
			mcol++;
			mcol2++;

			//subtract scaled row for the rest of the matrix, m
			while(mcol2 != mrow2->end())
			{
				*mcol2 -= (*mcol)*factor;				
				mcol++;
				mcol2++;				
			};

			mrow2++;
		};//end of performing row operations to set column (=rowNo) to zeroes below row = rowNo

		++rowNo;
	};//end of performing row operations to set lower left of matrix, m, to zeroes

	if(ans*0 != 0) return 0;

	return ans;
};

//! Solves matrix equation.
map<unsigned int, double> getSolnMatrixEqun(map<unsigned int, map<unsigned int, double> > & matrix, const map<unsigned int, double> & vect)
{
	if(vect.empty()) return vect;
	map<unsigned int, double> ans = vect;
	
	//loop thro' rows of matrix, m, and the inverse, i
	map<unsigned int, map<unsigned int, double> >::iterator mrow = matrix.begin();	
	map<unsigned int, map<unsigned int, double> >::iterator mrow2;	
	map<unsigned int, double>::iterator mcol;	
	map<unsigned int, double>::iterator mcol2;	
	map<unsigned int, double>::iterator ansIt = ans.begin();
	map<unsigned int, double>::iterator ansIt2;
	
	double factor;
	unsigned int rowNo = 1;
	unsigned int colNo = 1;

	for( ; mrow != matrix.end(); ++mrow, ++ansIt)
	{
		//set column to the first column in the row
		mcol = mrow->second.begin();		
		colNo = 1;

		//advance the column until the the row no. is equal to the column no.
		while(colNo != rowNo)
		{
			mcol++;			
			++colNo;
		};

		//divide the row in (m and i) by the value in the matrix, m, at (rowNo, colNo)
		factor = 1.0/(mcol->second); 
		mcol->second = 1;		
		mcol++;		

		//scale the remaining elements in the row - if there are any
		while(mcol != mrow->second.end())
		{
			mcol->second *= factor;			
			mcol++;			
		};

		//scale answer vector by the same factor
		ansIt->second *= factor;

		//subtract the row in question away from the remaining rows scaled  s.t. column = mrow will be zero below this row in matrix m
		mrow2 = mrow;
		ansIt2 = ansIt;		
		mrow2++;
		ansIt2++;
		
		//loop thro' remaining rows
		while(mrow2 != matrix.end())
		{
			//set column iterators to begining of the rows
			mcol2 = mrow2->second.begin();			
			mcol = mrow->second.begin();			

			//advance column of matrix, m, to the same as the main row being considered, with value rowNo
			colNo = 1;
			while(colNo != rowNo)
			{
				mcol++;
				mcol2++;				
				++colNo;
			};

			factor = mcol2->second; //factor to multiple row, rowNo by to take away from subseq. rows
			mcol2->second -= (mcol->second)*factor;
			mcol++;
			mcol2++;

			//subtract scaled row for the rest of the matrix, m
			while(mcol2 != mrow2->second.end())
			{
				mcol2->second -= (mcol->second)*factor;				
				mcol++;
				mcol2++;				
			};

			//subtract scaled value from ans vector
			ansIt2->second -= (ansIt->second)*factor;			

			mrow2++;
			ansIt2++;			
		};//end of performing row operations to set column (=rowNo) to zeroes below row = rowNo

		++rowNo;	
	};//end of performing row operations to set lower left of matrix, m, to zeroes

	//Now reduce the upper right of matrix, m, to zero
	map<unsigned int, map<unsigned int, double> >::reverse_iterator mrowre = matrix.rbegin();	
	map<unsigned int, map<unsigned int, double> >::reverse_iterator mrowre2 = matrix.rbegin();	
	map<unsigned int, double>::reverse_iterator ansItre = ans.rbegin();
	map<unsigned int, double>::reverse_iterator ansItre2;
	map<unsigned int, double>::reverse_iterator mcolre2;
	map<unsigned int, double>::reverse_iterator mcolre;

	rowNo = matrix.size();

	for( ; mrowre != matrix.rend(); ++mrowre, ++ansItre)
	{

		mrowre2 = mrowre;		
		ansItre2 = ansItre;
		mrowre2++;
		ansItre2++;		

		//loop tho' the remaining rows backwards - if there are any
		while(mrowre2 != matrix.rend())
		{			
			//set column iterators to begining of the rows
			mcolre2 = mrowre2->second.rbegin();			
			mcolre = mrowre->second.rbegin();		

			//advance column of matrix, m, to the same as the main row being considered, with value rowNo
			colNo = mrowre2->second.size();//size will be 4
			while(colNo != rowNo)
			{
				mcolre++;
				mcolre2++;				
				--colNo;
			};

			factor = mcolre2->second; //factor to multiple row, rowNo by to take away from subseq. rows
			mcolre2->second -= (mcolre->second)*factor;//0;
			mcolre++;
			mcolre2++;

			//subtract scaled row from the rest of the matrix, m
			while(mcolre2 != mrowre2->second.rend())//could stop at when col < row
			{
				mcolre2->second -= (mcolre->second)*factor;				
				mcolre++;
				mcolre2++;				
			};

			//subtract scaled value from ans vector
			ansItre2->second -= (ansItre->second)*factor;			

			mrowre2++;
			ansItre2++;			
		};

		--rowNo;
	};

	return ans;
};

//! Copies matrix from one to another.
void copyMatrix(const list< list<double> > & matrix, list< list<double> > & matrixCopy)
{
	
	if(freeClearMemory) list< list<double> >().swap(matrixCopy); else matrixCopy.clear();
	if(matrix.empty()) return;

	list<double> aRow;
	for(list< list<double> >::const_iterator m1 = matrix.begin(); m1 != matrix.end(); ++m1)
	{
		aRow.clear();
		for(list<double>::const_iterator m2 = m1->begin(); m2 != m1->end(); ++m2)
		{
			aRow.push_back(*m2);
		};
		matrixCopy.push_back(aRow);
	};

};

//! Calculates the p-value from a Chi square value with 1 df.
double getPvalueChiSq1DF(const double & chisq)
{	
	if(chisq*0 != 0) return 1.0;
	double a = sqrt(chisq)*oneOverSqRoot2;
	int ind = 0;
	return erfc1(&ind, &a);
};

//! Calculates the p-value for a given chi square value with df=df.
double getPvalueChiSq(const double & chisq, const double & df)
{
	double x = chisq;
	double df0 = df;
	double cum, ccum;

	cumchi(&x, &df0, &cum, &ccum);

	return ccum;
};

//! Calculates the q-value (1-p) for a given chi square value with df=df.
double getQvalueChiSq(const double & chisq, const double & df)
{
	double x = chisq;
	double df0 = df;
	double cum, ccum;

	cumchi(&x, &df0, &cum, &ccum);

	return cum;
};

//! Calculates the p-value from a f-statistic with d1 and d2 dfs.
double getPvalueFStat(double & fstat, const unsigned int & d1, const unsigned int & d2)
{	
	if(fstat <= 0) return 1.0;
	double x = (double)(d1*fstat)/(double)(d1*fstat + d2);
	double y = 1 - x;
	double a = (double)(d1)/2.0;
	double b = (double)(d2)/2.0;
	double w, w1;
	int ierr;

	bratio(&a, &b, &x, &y, &w, &w1, &ierr);

	return w1;
};

//! Calculates the p-value for a given standard normal z.
double getPvalueZSqd(double & zsqd)
{
	double integral;
	double compIntegral; // = 1 - integral	
	double z = sqrt(zsqd);

	cumnor(&z, &integral, &compIntegral); //computes integral of -inf to x of standard normal distribution

	return 2*compIntegral; //2 sided z-score p-value
};

//! Calculate the Z value corresponding to a p-value.
double calculateZSqdFromPvalue(const double & pval)
{
	double p = pval*0.5;
	double q = 1 - p;
	double z = dinvnr(&p, &q);
	return z*z;
};

// Calculate the Chi sq value (with 1 df) corresponding to a p-value using normal distribution.
double calculateChiSqFromPvalue(const double & pval)
{
	return calculateZSqdFromPvalue(pval);
};
