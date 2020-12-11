/************************************************************************
 * BayesNetty, version 1.1.1
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


/*! \file Model.h
    \brief This file contains the models used for linear regression.
    
*/

#ifndef __MODEL
#define __MODEL


#include <map>

using namespace std;

#include "Data.h"

//! General Class for Models.
class Model
{
protected:

	map<unsigned int, double> parameters;
	
	//used for fitting covariates
	CovariateData * covariateData;	

public:

	Model() : parameters(), covariateData(new CovariateData())  {};
	//setup parameter values and initValues of variables, to be over written in subclass of each model
	//where the variables created and added to the vector for the variables
	Model(map<unsigned int, double> paras) : parameters(paras), covariateData(new CovariateData()) {};
	
	virtual ~Model() {delete covariateData;};

	double getParameter(const unsigned int & no) const {return parameters.find(no)->second;};
	map<unsigned int, double> getParameters() const {return parameters;};
	void setNewParameters(map<unsigned int, double> & paras);
	CovariateData * getCovariateData() {return covariateData;};
	unsigned int getNoCovariates() {return covariateData->covariateDataAllSubjects.size();};

	virtual void setParameter(const unsigned int & no, const double & value) {parameters[no] = value;};	
	virtual double negLogLikelihood() {return 0;};
	virtual bool fitModel(double & rss, const bool & fitAnc, const bool & fitPart, const bool & fitInter) {return 0;};
};

//! Linear regression model.
class LinearRegModel : public Model
{
private:

	double variance;  //fitted variance from linear regression
	double mean;

public:

	LinearRegModel() : variance(0), mean(0) {};
	
	~LinearRegModel()
	{
		
	};

	bool fitModelCovar(double & rss);
	double negLogLikelihood() {return 0;};
	double getVariance() {return variance;};
	double getMean() {return mean;};
};

#endif
