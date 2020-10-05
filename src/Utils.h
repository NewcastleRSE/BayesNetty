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


/*! \file utils.h
    \brief This file defines some useful functions for calculating statistics. 
    
*/

#ifndef __UTILS
#define __UTILS

//! Converts an integer to a string
string toString(const unsigned int & i);

//! Converts a double to a string with 2 d.p. and pads out to 6 chars.
string toString2DP(const double & d, unsigned int minLength = 6);

//! Inverts a matrix. The inverse should be set to the identity when given to this function.
void getInverseMatrix(list< list<double> > & matrix, list< list<double> > & inverse);

//! Multiplies two matrices together, matrixProd should be empty when given to this function.
void getMatrixMulti(const list< list<double> > & matrix1, const list< list<double> > & matrix2, list< list<double> > & matrixProd);

//! Multiplies a matrix with a vector.
void getMatrixVecMulti(const list< list<double> > & matrix, const list<double> & vec, list<double> & vecAns);

//! Multiplies a matrix transpose with itself.
void getMatrixTransMatrixMulti(const list< list<double> > & matrix, list< list<double> > & matrixAns);

//! Multiplies a matrix transpose with another matrix.
void getMatrixTransMatrixMulti(const list< list<double> > & matrix1, const list< list<double> > & matrix2, list< list<double> > & matrixAns);

//! Multiplies a matrix transposed with a vector.
void getMatrixTransVecMulti(const list< list<double> > & matrix, const list<double> & vec, list<double> & vecAns);

//! Adds two matrix together.
void getMatrixAddMatrix(const list< list<double> > & matrix1, const list< list<double> > & matrix2, list< list<double> > & matrixAns);

//! Adds two vectors.
void getVecVecAdd(const list<double> & vec1, const list<double> & vec2, list<double> & ans);

//! Subtracts two vectors.
void getVecVecSub(const list<double> & vec1, const list<double> & vec2, list<double> & ans);

//! Multiplies a vector with a matrix.
void getVecMatrixMulti(const list<double> & vec, const list< list<double> > & matrix, list<double> & vecAns);

//! Multiplies a two vectors.
void getVecVecMulti(const list<double> & vec1, const list<double> & vec2, double & ans);

//! Computes determinant of square (perhaps symetric) matrix.
double getDetSquareMatrix(const list< list<double> > & sqMatrix);

//! Solves matrix equn Xz = y, returns z = X^{-1} y.
map<unsigned int, double> getSolnMatrixEqun(map<unsigned int, map<unsigned int, double> > & matrix, const map<unsigned int, double> & vect);

//! Copies matrix from one to another.
void copyMatrix(const list< list<double> > & matrix, list< list<double> > & matrixCopy);

//! Constant used in chi square calculations.
const double oneOverSqRoot2 = 0.7071067811865475244008443621048490392848359376884740365883;

//! Pi
const double pi = 3.14159265359;

//! Calculates the p-value for a given chi square value with df=1.
double getPvalueChiSq1DF(const double & chisq);

//! Calculates the p-value for a given chi square value with df=df.
double getPvalueChiSq(const double & chisq, const double & df);

//! Calculates the q-value (1-p) for a given chi square value with df=df.
double getQvalueChiSq(const double & chisq, const double & df);

//! Calculates the p-value from a f-statistic with d1 and d2 dfs.
double getPvalueFStat(double & fstat, const unsigned int & d1, const unsigned int & d2);

//! Calculates the p-value for a given standard normal z.
double getPvalueZSqd(double & zsqd);

//! Calculate the Z value corresponding to a p-value.
double calculateZSqdFromPvalue(const double & pval);

//! Calculate the Chi sq value (with 1 df) corresponding to a p-value.
double calculateChiSqFromPvalue(const double & pval);

#endif
