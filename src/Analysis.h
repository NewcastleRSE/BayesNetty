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


/*! \file Analysis.h
    \brief This file organises the Tasks and runs them.
    
*/

#ifndef __ANALYSIS
#define __ANALYSIS

#include <string>
#include <map>
#include <list>

using namespace std; // initiates the "std" or "standard" namespace

#include "main.h"
#include "Task.h"


//! Returns a string of the run time.
string getTime(const double & t);

//! Organises the correect analysis to perform.
class Analysis
{
private:

	list<CalculateMarkovBlanketTask *> calculateMarkovBlanketTasks;
	list<InputDataTask *> inputDataTasks;
	list<InputNetworkTask *> inputNetworkTasks;
	list<CalculatePosteriorTask *> calculatePosteriorTasks;
	list<CalculateNetworkScoreTask *> calculateNetworkScoreTasks;
	list<SearchNetworkModelsTask *> searchNetworkModelsTasks;
	list<AverageNetworksTask *> averageNetworksTasks;
	list<CompareNetworksTask *> compareNetworksTasks;
	list<SimulateNetworkDataTask *> simulateNetworkDataTasks;
	list<OutputNetworkTask *> outputNetworkTasks;
	list<OutputPriorsTask *> outputPriorsTasks;
	list<OutputPosteriorsTask *> outputPosteriorsTasks;
	list<ImputeNetworkDataTask *> imputeNetworkDataTasks;
	list<ImputeEstimateRecallPrecisionTask *> imputeEstimateRecallPrecisionTasks;
	list<CalculateRecallPrecisionTask *> calculateRecallPrecisionTasks;

	list<Task *> orderedTasks; //ordered task so that joint priors are done before posteriors etc
	set<string> taskNames;

	TaskCentral * taskCentral;

	//for setting up options
	ifstream parameterFile;
	string parameterFilename;
	bool usingParameterFile;
	int argcount;
	int argc;
	list<string> paraOptions;
	list<string>::const_iterator po;

	unsigned int randomSeed;

public:

	Analysis() : usingParameterFile(false), randomSeed(0), argcount(0), argc(0)
	{ 
		taskCentral = new TaskCentral();
	};

	//! Delete analysis things
	~Analysis()
	{
		for(list<Task *>::iterator ot = orderedTasks.begin(); ot != orderedTasks.end(); ++ot) delete *ot;
		orderedTasks.clear();		
		paraOptions.clear();
		delete taskCentral;
	};

	//set up Tasks from options methods
	void processParameterFile(string & paraFilename);
	void setCommandArgs(int & argcount0, int & argc0, char * argv0[]);	
	void updateArgcount(int & argcount0) {++argcount; ++po; argcount0 = argcount;};
	bool processOption(string & option);
	string getStringOptionValue();
	unsigned int getUIntOptionValue();
	double getDoubleOptionValue();
	void addTaskName(const string & taskName);
	void setRandomSeed(const unsigned int & seed) {randomSeed = seed; taskCentral->setRandomSeed(seed);};
	unsigned int getRandomSeed() {return randomSeed;};
	
	void doTasks();
	CalculateMarkovBlanketTask * getLatestCalculateMarkovBlanketTask();
	InputDataTask * getLatestInputDataTask();
	InputNetworkTask * getLatestInputNetworkTask();
	CalculatePosteriorTask * getLatestCalculatePosteriorTask();
	CalculateNetworkScoreTask * getLatestCalculateNetworkScoreTask();
	SearchNetworkModelsTask * getLatestSearchNetworkModelsTask();
	AverageNetworksTask * getLatestAverageNetworksTask();
	CompareNetworksTask * getLatestCompareNetworksTask();
	SimulateNetworkDataTask * getLatestSimulateNetworkDataTask();
	OutputNetworkTask * getLatestOutputNetworkTask();
	OutputPriorsTask * getLatestOutputPriorsTask();
	OutputPosteriorsTask * getLatestOutputPosteriorsTask();
	ImputeNetworkDataTask * getLatestImputeNetworkDataTask();
	ImputeEstimateRecallPrecisionTask * getLatestImputeEstimateRecallPrecisionTask();
	CalculateRecallPrecisionTask * getLatestCalculateRecallPrecisionTask();
};


#endif
