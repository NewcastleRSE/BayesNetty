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

/*! \file Model.cpp
    \brief This file contains the source of models used for logistic regression.
    
*/

#include <map>
#include <iostream>
#include <math.h>

using namespace std; // initiates the "std" or "standard" namespace

#include "Model.h"
#include "Data.h"
#include "Utils.h"

//! Sets the values of all of the parameters.
void Model::setNewParameters(map<unsigned int, double> & paras)
{
	for(map<unsigned int, double>::const_iterator p = paras.begin(); p != paras.end(); ++p)
	{
		setParameter(p->first, p->second);
	};
};

//! Fits 2 SNP model to quantitive traits.
bool LinearRegModel::fitModelCovar(double & rss)
{
	//Solve equn X^T y = (X^T X)betaHat, to find betaHat, where X is the design matrix 
	
	//construct vector X^T y and matrix (X^T X)
	map<unsigned int, double> vectorXTy;
	map<unsigned int, map<unsigned int, double> > matrixXTX;
	list<bool>::const_iterator cc = covariateData->notInGroupOrMissingData.begin();
	
	//list of iterators, one for each cts parent
	list<list<double>::const_iterator> ctsPas;

	list<double>::const_iterator aPa;

	//setup iterators for cts parents
	for(list<CtsData *>::const_iterator ctsPa = covariateData->covariateDataAllSubjects.begin(); ctsPa != covariateData->covariateDataAllSubjects.end(); ++ctsPa)
	{
		aPa = (*ctsPa)->values.begin();
		ctsPas.push_back(aPa);
	};

	double totalQT = 0;
	list<double> covarQT, sumCV;
	list<double>::iterator cvQT, sCV;
	for(list<CtsData *>::const_iterator ctsPa = covariateData->covariateDataAllSubjects.begin(); ctsPa != covariateData->covariateDataAllSubjects.end(); ++ctsPa)
	{
		covarQT.push_back(0); //set initial totals to 0
		sumCV.push_back(0);	
	};

	list<list<double> > covarCovarSums;
	list<double> aCovarCovarSum;
	
	for(list<CtsData *>::const_iterator ctsPa = covariateData->covariateDataAllSubjects.begin(); ctsPa != covariateData->covariateDataAllSubjects.end(); ++ctsPa)
	{
		aCovarCovarSum.push_back(0); //set initial totals to 0			
	};

	for(list<CtsData *>::const_iterator ctsPa = covariateData->covariateDataAllSubjects.begin(); ctsPa != covariateData->covariateDataAllSubjects.end(); ++ctsPa)
	{
		covarCovarSums.push_back(aCovarCovarSum);
	};

	list<list<double> >::iterator cvcvSums;
	list<double>::iterator cvcvSums2;
	
	unsigned int covariateParameterNo;

	unsigned int totalNotMissing = 0;

	for(list<double>::const_iterator i = covariateData->nodeValues->values.begin(); i != covariateData->nodeValues->values.end(); ++i, ++cc)
	{

		if(!*cc) // miss out indivs where there is missing SNP data or missing QT data
		{
			totalQT += *i;		
			totalNotMissing++;			
	
			//do covariate sums		
			cvQT = covarQT.begin();
			sCV = sumCV.begin();		
			cvcvSums = covarCovarSums.begin();	

			//need to loop thro' data for each parent
			for(list<list<double>::const_iterator>::const_iterator aCtsPa = ctsPas.begin(); aCtsPa != ctsPas.end(); ++aCtsPa)
			{				
				*cvQT += (*i)*(**aCtsPa);
				*sCV += (**aCtsPa);

				//loop thro' covar*covar sums
				cvcvSums2 = cvcvSums->begin();				
				for(list<list<double>::const_iterator>::const_iterator aCtsPa2 = ctsPas.begin(); aCtsPa2 != ctsPas.end(); ++aCtsPa2, ++cvcvSums2)
				{					
					*cvcvSums2 += (**aCtsPa)*(**aCtsPa2);
				};

				++cvQT; ++sCV; 
				++cvcvSums;
			};

		};

		//advance cts parent iterators to next subject data
		for(list<list<double>::const_iterator>::iterator aCtsPa = ctsPas.begin(); aCtsPa != ctsPas.end(); ++aCtsPa)
		{				
			(*aCtsPa)++;
		};
	};

	//construct vector XTy
	vectorXTy[1] = totalQT;

	//do covariate elements, row 2 onwards
	for(covariateParameterNo = 2, cvQT = covarQT.begin(); cvQT != covarQT.end(); ++cvQT, ++covariateParameterNo)
	{
		vectorXTy[covariateParameterNo] = *cvQT;
	};

	//construct matrix XTX
	//row 1
	map<unsigned int, double> rowOne;
	rowOne[1] = totalNotMissing;

	//do covariate elements, row 2 onwards
	for(covariateParameterNo = 2, sCV = sumCV.begin(); sCV != sumCV.end(); ++sCV, ++covariateParameterNo)
	{
		rowOne[covariateParameterNo] = *sCV;	
	};

	matrixXTX[1] = rowOne;

	//do covariate rows, row 5 onwards
	sCV = sumCV.begin();
	map<unsigned int, double> rowCovariate;
	unsigned int covariateParameterNo2;
	cvcvSums = covarCovarSums.begin();

	for(covariateParameterNo = 2; sCV != sumCV.end(); ++covariateParameterNo, ++cvcvSums)
	{
		rowCovariate[1] = *sCV; 
	
		for(covariateParameterNo2 = 2, cvcvSums2 = cvcvSums->begin(); cvcvSums2 != cvcvSums->end(); ++cvcvSums2, ++covariateParameterNo2)
		{
			rowCovariate[covariateParameterNo2] = *cvcvSums2;
		};
		matrixXTX[covariateParameterNo] = rowCovariate;
		++sCV;
	};

	//solve matrix equation to get betaHat
	map<unsigned int, double> ans = getSolnMatrixEqun(matrixXTX, vectorXTy);

	//intercept = 0;

	//set parameters
	for(map<unsigned int, double>::const_iterator par = ans.begin(); par != ans.end(); ++par)
	{
		setParameter(par->first, par->second);	
	};

	//set residual sum of squares	
	cc = covariateData->notInGroupOrMissingData.begin();

	double diff;
	double xiTbetaHat; //model fit

	rss = 0;

	list<CtsData *>::const_iterator ctsPa = covariateData->covariateDataAllSubjects.begin();
	//reset iterators for cts parents
	for(list<list<double>::const_iterator>::iterator aCtsPa = ctsPas.begin(); aCtsPa != ctsPas.end(); ++aCtsPa, ++ctsPa)
	{				
		*aCtsPa = (*ctsPa)->values.begin();
	};

	//loop thro' subjects and calc difference of model fit with the observed data.
	for(list<double>::const_iterator i = covariateData->nodeValues->values.begin(); i != covariateData->nodeValues->values.end(); ++i, ++cc)
	{	
		// miss out indivs where there is missing SNP data or missing QT data
		if(!*cc)
		{
			xiTbetaHat = getParameter(1); //intercept;			
			covariateParameterNo = 2;
			for(list<list<double>::const_iterator>::const_iterator aCtsPa = ctsPas.begin(); aCtsPa != ctsPas.end(); ++aCtsPa)
			{
				xiTbetaHat += (**aCtsPa)*(getParameter(covariateParameterNo));
				covariateParameterNo++;
			};

			diff = *i - xiTbetaHat;
			rss += diff*diff;
		};
		
		//advance cts parent iterators to next subject data
		for(list<list<double>::const_iterator>::iterator aCtsPa = ctsPas.begin(); aCtsPa != ctsPas.end(); ++aCtsPa)
		{				
			(*aCtsPa)++;
		};
	};

	double denom = ((double)(totalNotMissing - covariateData->covariateDataAllSubjects.size() - 1));
	if(denom==0) variance = 0;
	else variance = rss/denom; // = rss/(n-p)


	//use max like as in deal, biased estimator of variance
	if(variance < 0) variance = rss/((double)(totalNotMissing)); // = rss/(n)

	mean = totalQT/((double)(totalNotMissing));

	//uncomment to use all data size (incorrectly) as in deal
	//variance = rss/((double)(covariateData->nodeValues->values.size()));

	return rss*0 == 0 && variance != 0;
};
