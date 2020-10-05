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


/*! \file Analysis.cpp
    \brief This file contains the methods the various analyse.
    
*/
#include <iostream>
#include <sstream>
#include <fstream>
#include <list>

using namespace std; // initiates the "std" or "standard" namespace

#include "main.h"
#include "Analysis.h"


//! Converts an integer to a string
string toString(int & i)
{
	ostringstream aStringStream;
	aStringStream << i;

	return aStringStream.str();
};

//! Returns a string of the run time
string getTime(const double & t)
{
	double time = t;
	int days = 0;
	int hours = 0;
	int minutes = 0;
	int seconds = 0;

	string ans = "";
	days = (int) (time / 86400); time -= (double)days*86400;
	hours = (int) (time / 3600); time -= (double)hours*3600;
	minutes = (int) (time / 60); time -= (double)minutes*60;
	seconds = (int) time;

	if(days == 1) ans += "1 day";
	else if(days > 0) ans += toString(days) + " days";

	if(hours > 0)
	{
		if(days != 0)
		{
			if(minutes == 0 && seconds == 0) ans += " and ";
			else ans += ", ";
		};

		if(hours == 1) ans += "1 hour";
		else ans += toString(hours) + " hours";
	};

	if(minutes > 0)
	{
		if(ans != "")
		{
			if(seconds == 0) ans += " and ";
			else ans += ", ";
		};

		if(minutes == 1) ans += "1 minute";
		else ans += toString(minutes) + " minutes";
	};

	if(seconds > 0)
	{
		if(ans != "")
		{
			ans += " and ";			
		};

		if(seconds == 1) ans += "1 second";
		else ans += toString(seconds) + " seconds";
	};

	if(ans == "") ans = "less than one second";

	return ans;
};

//! Processes a parameter file and creates tasks to do.
void Analysis::processParameterFile(string & paraFilename)
{
	parameterFilename = paraFilename; //set for use in error reporting
	
	usingParameterFile = true;
	parameterFile.open(parameterFilename.c_str());
	if(!parameterFile.is_open())
	{
		header();
		string message = "Cannot read parameter file: " + parameterFilename + "!";
		exitErr(message);
	};

	string word;

	parameterFile >> word;

	while(!parameterFile.eof()){		
		
		if(word.length() >= 2 && word.substr(0,1) == "-" && !processOption(word))
		{
			header();
    		string message = "Unrecognised option, " + word + ", in parameter file: " + parameterFilename + "!";
			exitErr(message);
		};

		parameterFile >> word;
	};

	parameterFile.close();
	usingParameterFile = false;
};

//! Get the latest added task for this type of task.
CalculateMarkovBlanketTask * Analysis::getLatestCalculateMarkovBlanketTask()
{
	list<CalculateMarkovBlanketTask *>::const_reverse_iterator t = calculateMarkovBlanketTasks.rbegin();

	//create new task if none exists
	if(t == calculateMarkovBlanketTasks.rend())
	{
		CalculateMarkovBlanketTask * aCalculateMarkovBlanketTask = new CalculateMarkovBlanketTask();
		calculateMarkovBlanketTasks.push_back(aCalculateMarkovBlanketTask);
		orderedTasks.push_back(aCalculateMarkovBlanketTask);
		return aCalculateMarkovBlanketTask;
	};

	return *t;
};

//! Get the latest added task for this type of task.
InputDataTask * Analysis::getLatestInputDataTask()
{
	list<InputDataTask *>::const_reverse_iterator t = inputDataTasks.rbegin();

	//create new task if none exists
	if(t == inputDataTasks.rend())
	{
		InputDataTask * aInputDataTask = new InputDataTask();
		inputDataTasks.push_back(aInputDataTask);
		orderedTasks.push_back(aInputDataTask);
		return aInputDataTask;
	};

	return *t;
};

//! Get the latest added task for this type of task.
InputNetworkTask * Analysis::getLatestInputNetworkTask()
{
	list<InputNetworkTask *>::const_reverse_iterator t = inputNetworkTasks.rbegin();

	//create new task if none exists
	if(t == inputNetworkTasks.rend())
	{
		InputNetworkTask * aInputNetworkTask = new InputNetworkTask();
		inputNetworkTasks.push_back(aInputNetworkTask);
		orderedTasks.push_back(aInputNetworkTask);
		return aInputNetworkTask;
	};

	return *t;
};

//! Get the latest added task for this type of task.
CalculatePosteriorTask * Analysis::getLatestCalculatePosteriorTask()
{
	list<CalculatePosteriorTask *>::const_reverse_iterator t = calculatePosteriorTasks.rbegin();

	//create new task if none exists
	if(t == calculatePosteriorTasks.rend())
	{
		CalculatePosteriorTask * aCalculatePosteriorTask = new CalculatePosteriorTask();
		calculatePosteriorTasks.push_back(aCalculatePosteriorTask);
		orderedTasks.push_back(aCalculatePosteriorTask);
		return aCalculatePosteriorTask;
	};

	return *t;
};

//! Get the latest added task for this type of task.
CalculateNetworkScoreTask * Analysis::getLatestCalculateNetworkScoreTask()
{
	list<CalculateNetworkScoreTask *>::const_reverse_iterator t = calculateNetworkScoreTasks.rbegin();

	//create new task if none exists
	if(t == calculateNetworkScoreTasks.rend())
	{
		CalculateNetworkScoreTask * aCalculateNetworkScoreTask = new CalculateNetworkScoreTask();
		calculateNetworkScoreTasks.push_back(aCalculateNetworkScoreTask);
		orderedTasks.push_back(aCalculateNetworkScoreTask);
		return aCalculateNetworkScoreTask;
	};

	return *t;
};

//! Gets the latest task of this type.
SearchNetworkModelsTask * Analysis::getLatestSearchNetworkModelsTask()
{
	list<SearchNetworkModelsTask *>::const_reverse_iterator t = searchNetworkModelsTasks.rbegin();

	//create new task if none exists
	if(t == searchNetworkModelsTasks.rend())
	{
		SearchNetworkModelsTask * aSearchNetworkModelsTask = new SearchNetworkModelsTask();
		searchNetworkModelsTasks.push_back(aSearchNetworkModelsTask);
		orderedTasks.push_back(aSearchNetworkModelsTask);
		return aSearchNetworkModelsTask;
	};

	return *t;
};

//! Gets the latest task of this type.
AverageNetworksTask * Analysis::getLatestAverageNetworksTask()
{
	list<AverageNetworksTask *>::const_reverse_iterator t = averageNetworksTasks.rbegin();

	//create new task if none exists
	if(t == averageNetworksTasks.rend())
	{
		AverageNetworksTask * aAverageNetworksTask = new AverageNetworksTask();
		averageNetworksTasks.push_back(aAverageNetworksTask);
		orderedTasks.push_back(aAverageNetworksTask);
		return aAverageNetworksTask;
	};

	return *t;
};

//! Gets the latest task of this type.
CompareNetworksTask * Analysis::getLatestCompareNetworksTask()
{
	list<CompareNetworksTask *>::const_reverse_iterator t = compareNetworksTasks.rbegin();

	//create new task if none exists
	if(t == compareNetworksTasks.rend())
	{
		CompareNetworksTask * aCompareNetworksTask = new CompareNetworksTask();
		compareNetworksTasks.push_back(aCompareNetworksTask);
		orderedTasks.push_back(aCompareNetworksTask);
		return aCompareNetworksTask;
	};

	return *t;
};

//! Gets the latest task of this type.
SimulateNetworkDataTask * Analysis::getLatestSimulateNetworkDataTask()
{
	list<SimulateNetworkDataTask *>::const_reverse_iterator t = simulateNetworkDataTasks.rbegin();

	//create new task if none exists
	if(t == simulateNetworkDataTasks.rend())
	{
		SimulateNetworkDataTask * aSimulateNetworkDataTask = new SimulateNetworkDataTask();
		simulateNetworkDataTasks.push_back(aSimulateNetworkDataTask);
		orderedTasks.push_back(aSimulateNetworkDataTask);
		return aSimulateNetworkDataTask;
	};

	return *t;
};

//! Gets the latest task of this type.
OutputNetworkTask * Analysis::getLatestOutputNetworkTask()
{
	list<OutputNetworkTask *>::const_reverse_iterator t = outputNetworkTasks.rbegin();

	//create new task if none exists
	if(t == outputNetworkTasks.rend())
	{
		OutputNetworkTask * aOutputNetworkTask = new OutputNetworkTask();
		outputNetworkTasks.push_back(aOutputNetworkTask);
		orderedTasks.push_back(aOutputNetworkTask);
		return aOutputNetworkTask;
	};

	return *t;
};

//! Gets the latest task of this type.
OutputPriorsTask * Analysis::getLatestOutputPriorsTask()
{
	list<OutputPriorsTask *>::const_reverse_iterator t = outputPriorsTasks.rbegin();

	//create new task if none exists
	if(t == outputPriorsTasks.rend())
	{
		OutputPriorsTask * aOutputPriorsTask = new OutputPriorsTask();
		outputPriorsTasks.push_back(aOutputPriorsTask);
		orderedTasks.push_back(aOutputPriorsTask);
		return aOutputPriorsTask;
	};

	return *t;
};

//! Gets the latest task of this type.
OutputPosteriorsTask * Analysis::getLatestOutputPosteriorsTask()
{
	list<OutputPosteriorsTask *>::const_reverse_iterator t = outputPosteriorsTasks.rbegin();

	//create new task if none exists
	if(t == outputPosteriorsTasks.rend())
	{
		OutputPosteriorsTask * aOutputPosteriorsTask = new OutputPosteriorsTask();
		outputPosteriorsTasks.push_back(aOutputPosteriorsTask);
		orderedTasks.push_back(aOutputPosteriorsTask);
		return aOutputPosteriorsTask;
	};

	return *t;
};

//! Gets the latest task of this type.
ImputeNetworkDataTask * Analysis::getLatestImputeNetworkDataTask()
{
	list<ImputeNetworkDataTask *>::const_reverse_iterator t = imputeNetworkDataTasks.rbegin();

	//create new task if none exists
	if(t == imputeNetworkDataTasks.rend())
	{
		ImputeNetworkDataTask * aImputeNetworkDataTask = new ImputeNetworkDataTask();
		imputeNetworkDataTasks.push_back(aImputeNetworkDataTask);
		orderedTasks.push_back(aImputeNetworkDataTask);
		return aImputeNetworkDataTask;
	};

	return *t;
};

//! Gets the latest task of this type.
ImputeEstimateRecallPrecisionTask * Analysis::getLatestImputeEstimateRecallPrecisionTask()
{
	list<ImputeEstimateRecallPrecisionTask *>::const_reverse_iterator t = imputeEstimateRecallPrecisionTasks.rbegin();

	//create new task if none exists
	if(t == imputeEstimateRecallPrecisionTasks.rend())
	{
		ImputeEstimateRecallPrecisionTask * aImputeEstimateRecallPrecisionTask = new ImputeEstimateRecallPrecisionTask();
		imputeEstimateRecallPrecisionTasks.push_back(aImputeEstimateRecallPrecisionTask);
		orderedTasks.push_back(aImputeEstimateRecallPrecisionTask);
		return aImputeEstimateRecallPrecisionTask;
	};

	return *t;
};

//! Gets the latest task of this type.
CalculateRecallPrecisionTask * Analysis::getLatestCalculateRecallPrecisionTask()
{
	list<CalculateRecallPrecisionTask *>::const_reverse_iterator t = calculateRecallPrecisionTasks.rbegin();

	//create new task if none exists
	if(t == calculateRecallPrecisionTasks.rend())
	{
		CalculateRecallPrecisionTask * aCalculateRecallPrecisionTask = new CalculateRecallPrecisionTask();
		calculateRecallPrecisionTasks.push_back(aCalculateRecallPrecisionTask);
		orderedTasks.push_back(aCalculateRecallPrecisionTask);
		return aCalculateRecallPrecisionTask;
	};

	return *t;
};

//! Get the next unsigned int in the parameter file or command line.
unsigned int Analysis::getUIntOptionValue()
{
	unsigned int anUnInt;

	if(usingParameterFile)
	{
		if(parameterFile.eof())
		{
			header();
			string mess = "Missing option(s)! Came to the end of the parameter, " + parameterFilename + ", too soon!";
			exitErr(mess);
		};
		
		parameterFile >> anUnInt;
	}
	else
	{
		argcount++; po++;
		if(argcount >= argc || po == paraOptions.end()) {header(); exitErr("Missing option(s)! Came to the end of the command line too soon!");};
		anUnInt = atoi((*po).c_str());
	};

	return anUnInt;
};

//! Get the next double in the parameter file or command line.
double Analysis::getDoubleOptionValue()
{
	double aDouble;

	if(usingParameterFile)
	{
		if(parameterFile.eof())
		{
			header();
			string mess = "Missing option(s)! Came to the end of the parameter, " + parameterFilename + ", too soon!";
			exitErr(mess);
		};
		
		parameterFile >> aDouble;
	}
	else
	{
		argcount++; po++;
		if(argcount >= argc || po == paraOptions.end()) {header(); exitErr("Missing option(s)! Came to the end of the command line too soon!");};		
		aDouble = atof((*po).c_str());	
	};

	return aDouble;
};

//! Get the next string in the parameter file or command line.
string Analysis::getStringOptionValue()
{
	string aString;

	if(usingParameterFile)
	{
		parameterFile >> aString;

		if(parameterFile.eof())
		{
			header();
			string mess = "Missing option(s)! Came to the end of the parameter file, " + parameterFilename + ", too soon!";
			exitErr(mess);
		};	
	}
	else
	{
		argcount++; po++;
		if(argcount > argc || po == paraOptions.end()) {header(); exitErr("Missing option(s)! Came to the end of the command line too soon!");};		
		aString = *po;
	};

	return aString;
};

//! Copy command line arguments.
void Analysis::setCommandArgs(int & argcount0, int & argc0, char * argv0[])
{
	argcount = argcount0;
	argc = argc0; 	
	string aOp;

	for (int i = 1; i <= argc; ++i)
	{	
		paraOptions.push_back(argv0[i - 1]); //copy contents
	};

	po = paraOptions.begin(); ++po;
};


//! Process option in parameter file.
bool Analysis::processOption(string & option)
{
	bool validOption1 = true;
	bool validOption2 = true;
	bool validOption3 = true;

	if(option == "-log") logFilename = getStringOptionValue();			
	else if(option == "-so") outputToScreen = false;
	else if(option == "-free-memory") freeClearMemory = true;
	else if(option == "-seed") {setRandomSeed(getUIntOptionValue());}
	
	else if(option == "-input-data") {InputDataTask * idt = new InputDataTask(); inputDataTasks.push_back(idt); orderedTasks.push_back(idt);}
	else if(option == "-input-data-name") getLatestInputDataTask()->setTaskName(getStringOptionValue());
	else if(option == "-input-data-file") getLatestInputDataTask()->addDataFile(getStringOptionValue());
	else if(option == "-input-data-include-file") getLatestInputDataTask()->setIncludeVariables(getStringOptionValue());
	else if(option == "-input-data-exclude-file") getLatestInputDataTask()->setExcludeVariables(getStringOptionValue());
	else if(option == "-input-data-cts") getLatestInputDataTask()->setDataType(1);
	else if(option == "-input-data-discrete") getLatestInputDataTask()->setDataType(2);
	else if(option == "-input-data-discrete-snp") getLatestInputDataTask()->setDataType(3);
	else if(option == "-input-data-cts-snp") getLatestInputDataTask()->setDataType(4);
	else if(option == "-input-data-factor") getLatestInputDataTask()->setDataType(5);
	else if(option == "-input-data-factor-snp") getLatestInputDataTask()->setDataType(6);
	else if(option == "-input-data-cts-snp2") getLatestInputDataTask()->setDataType(7);
	else if(option == "-input-data-discrete-snp2") getLatestInputDataTask()->setDataType(8);
	else if(option == "-input-data-factor-snp2") getLatestInputDataTask()->setDataType(9);
	else if(option == "-input-data-cts-missing-value") getLatestInputDataTask()->setCtsMissingValue(getDoubleOptionValue());
	else if(option == "-input-data-discrete-missing-value" || option == "-input-data-factor-missing-value" ) getLatestInputDataTask()->setDiscreteMissingValue(getStringOptionValue());
	else if(option == "-input-data-ids") getLatestInputDataTask()->setNumberOfIDs(getUIntOptionValue());
	else if(option == "-input-data-csv") getLatestInputDataTask()->setIsCSV(true);
	else if(option == "-input-data-prefix-name") getLatestInputDataTask()->setPrefixName(getStringOptionValue());

	else if(option == "-input-network") {InputNetworkTask * intk = new InputNetworkTask(); inputNetworkTasks.push_back(intk); orderedTasks.push_back(intk);}
	else if(option == "-input-network-name") getLatestInputNetworkTask()->setTaskName(getStringOptionValue());
	else if(option == "-input-network-type") getLatestInputNetworkTask()->setNetworkType(getStringOptionValue());
	else if(option == "-input-network-file") getLatestInputNetworkTask()->addNetworkFile(getStringOptionValue());
	else if(option == "-input-network-file1") getLatestInputNetworkTask()->addNetworkFile(getStringOptionValue());
	else if(option == "-input-network-file2") getLatestInputNetworkTask()->addNetworkFile2(getStringOptionValue());
	else if(option == "-input-network-igraph-file-prefix") getLatestInputNetworkTask()->addNetworkFileIGraph(getStringOptionValue());
	else if(option == "-input-network-empty") getLatestInputNetworkTask()->setEmptyNetwork(true);
	else if(option == "-input-network-using-imputed-data") getLatestInputNetworkTask()->setUsingImputedData();
	else if(option == "-input-network-whitelist-file") getLatestInputNetworkTask()->addWhitelistDataFile(getStringOptionValue());
	else if(option == "-input-network-blacklist-file") getLatestInputNetworkTask()->addBlacklistDataFile(getStringOptionValue());
	else if(option == "-input-network-imaginary-sample-size") getLatestInputNetworkTask()->setImaginarySampleSize(getDoubleOptionValue());
	else if(option == "-input-network-score") getLatestInputNetworkTask()->setScore(getStringOptionValue());
	else if(option == "-input-network-score-fix") getLatestInputNetworkTask()->setScoreFix(getStringOptionValue());
	else if(option == "-input-network-no-parents-node") getLatestInputNetworkTask()->addNoParentsNode(getStringOptionValue());
	else if(option == "-input-network-no-children-node") getLatestInputNetworkTask()->addNoChildrenNode(getStringOptionValue());
	else if(option == "-input-network-blacklist-edge-type")
	{
		string nodeType1 = getStringOptionValue();
		string nodeType2 = getStringOptionValue();
		getLatestInputNetworkTask()->setNotAllowedEdgeType(nodeType1, nodeType2);
	}
	else if(/*option == "-input-network-cost-edge" ||*/ option == "-input-network-prob-edge")
	{
		string node1 = getStringOptionValue();
		string node2 = getStringOptionValue();
		double cost = getDoubleOptionValue();
		if(cost < 0 || cost > 1)
		{
			exitErr("The edge probability must be set to between 0 and 1!");
		};
		getLatestInputNetworkTask()->setCostEdge(node1, node2, cost);
		getLatestInputNetworkTask()->setCostEdge(node2, node1, 1.0 - cost);
		getLatestInputNetworkTask()->setScore("BICprob");
	}
	else if(/*option == "-input-network-cost-edge-type" ||*/ option == "-input-network-prob-edge-type")
	{
		string nodeType1 = getStringOptionValue();
		string nodeType2 = getStringOptionValue();
		double cost = getDoubleOptionValue();
		if(cost < 0 || cost > 1)
		{
			exitErr("The edge probability must be set to between 0 and 1!");
		};
		getLatestInputNetworkTask()->setCostEdgeType(nodeType1, nodeType2, cost);
		getLatestInputNetworkTask()->setCostEdgeType(nodeType2, nodeType1, 1.0 - cost);
		getLatestInputNetworkTask()->setScore("BICprob");
	}
	else
	{
		validOption1 = false;
	};

	if(option == "-calc-posterior") {CalculatePosteriorTask * cpt = new CalculatePosteriorTask(); calculatePosteriorTasks.push_back(cpt); orderedTasks.push_back(cpt);}
	else if(option == "-calc-posterior-name") getLatestCalculatePosteriorTask()->setTaskName(getStringOptionValue());
	else if(option == "-calc-posterior-network-name") getLatestCalculatePosteriorTask()->setNetworkName(getStringOptionValue());

	else if(option == "-calc-network-score") {CalculateNetworkScoreTask * cnst = new CalculateNetworkScoreTask(); calculateNetworkScoreTasks.push_back(cnst); orderedTasks.push_back(cnst);}
	else if(option == "-calc-network-score-name") getLatestCalculateNetworkScoreTask()->setTaskName(getStringOptionValue());
	else if(option == "-calc-network-score-network-name") getLatestCalculateNetworkScoreTask()->setNetworkName(getStringOptionValue());
	else if(option == "-calc-network-score-file") getLatestCalculateNetworkScoreTask()->setFilename(getStringOptionValue());
	else if(option == "-calc-network-score-all-scores") getLatestCalculateNetworkScoreTask()->setAllScoresFilename(getStringOptionValue());

	else if(option == "-markov-blanket") {CalculateMarkovBlanketTask * cnst = new CalculateMarkovBlanketTask(); calculateMarkovBlanketTasks.push_back(cnst); orderedTasks.push_back(cnst);}
	else if(option == "-markov-blanket-name") getLatestCalculateMarkovBlanketTask()->setTaskName(getStringOptionValue());
	else if(option == "-markov-blanket-network-name") getLatestCalculateMarkovBlanketTask()->setNetworkName(getStringOptionValue());
	else if(option == "-markov-blanket-node-name") getLatestCalculateMarkovBlanketTask()->setNodeName(getStringOptionValue());

	else if(option == "-search-models") {SearchNetworkModelsTask * snmt = new SearchNetworkModelsTask(); searchNetworkModelsTasks.push_back(snmt); orderedTasks.push_back(snmt);}
	else if(option == "-search-models-name") getLatestSearchNetworkModelsTask()->setTaskName(getStringOptionValue());
	else if(option == "-search-models-network-name") getLatestSearchNetworkModelsTask()->setNetworkName(getStringOptionValue());
	else if(option == "-search-models-file") getLatestSearchNetworkModelsTask()->setFilename(getStringOptionValue());
	else if(option == "-search-models-random-restarts") getLatestSearchNetworkModelsTask()->setRandomRestarts(getUIntOptionValue());
	else if(option == "-search-models-jitter-restarts") getLatestSearchNetworkModelsTask()->setJitterRestarts(getUIntOptionValue());
	
	else if(option == "-average-networks") {AverageNetworksTask * snmt = new AverageNetworksTask(); averageNetworksTasks.push_back(snmt); orderedTasks.push_back(snmt);}
	else if(option == "-average-networks-name") getLatestAverageNetworksTask()->setTaskName(getStringOptionValue());
	else if(option == "-average-networks-network-name") getLatestAverageNetworksTask()->setNetworkName(getStringOptionValue());
	else if(option == "-average-networks-file") getLatestAverageNetworksTask()->setFilename(getStringOptionValue());
	else if(option == "-average-networks-fam-file") getLatestAverageNetworksTask()->setFamilyFile(getStringOptionValue());
	else if(option == "-average-networks-igraph-file-prefix") getLatestAverageNetworksTask()->setRFilenamePrefix(getStringOptionValue());
	else if(option == "-average-networks-threshold") getLatestAverageNetworksTask()->setArcThreshold(getDoubleOptionValue());
	else if(option == "-average-networks-bootstraps") getLatestAverageNetworksTask()->setNoBootstraps(getUIntOptionValue());
	else if(option == "-average-networks-random-restarts") getLatestAverageNetworksTask()->setRandomRestarts(getUIntOptionValue());
	else if(option == "-average-networks-jitter-restarts") getLatestAverageNetworksTask()->setJitterRestarts(getUIntOptionValue());
	else if(option == "-average-networks-likelihood-file") getLatestAverageNetworksTask()->setLikelihoodFilename(getStringOptionValue());
	else if(option == "-average-networks-use-score-method") getLatestAverageNetworksTask()->setUseNetworkScoreMethod(true);
	else if(option == "-average-networks-use-weight-method") getLatestAverageNetworksTask()->setUseNetworkWeightMethod(true);
	else
	{
		validOption2 = false;
	};

	//else if(option == "-compare-networks") {CompareNetworksTask * cnt = new CompareNetworksTask(); compareNetworksTasks.push_back(cnt); orderedTasks.push_back(cnt);}
	//else if(option == "-compare-networks-name") getLatestCompareNetworksTask()->setTaskName(getStringOptionValue());

	if(option == "-simulate-network-data") {SimulateNetworkDataTask * pnt = new SimulateNetworkDataTask(); simulateNetworkDataTasks.push_back(pnt); orderedTasks.push_back(pnt);}
	else if(option == "-simulate-network-data-name") getLatestSimulateNetworkDataTask()->setTaskName(getStringOptionValue());
	else if(option == "-simulate-network-data-network-name") getLatestSimulateNetworkDataTask()->setNetworkName(getStringOptionValue());
	else if(option == "-simulate-network-data-no-sims") getLatestSimulateNetworkDataTask()->setNoSims(getUIntOptionValue());
	else if(option == "-simulate-network-data-parameter-file") getLatestSimulateNetworkDataTask()->setParameterFile(getStringOptionValue());
	else if(option == "-simulate-network-data-snp-integers") getLatestSimulateNetworkDataTask()->setSNPIntegers(true);
	else if(option == "-simulate-network-data-whitelist-file") getLatestSimulateNetworkDataTask()->addWhitelistDataFile(getStringOptionValue());
	else if(option == "-simulate-network-data-blacklist-file") getLatestSimulateNetworkDataTask()->addBlacklistDataFile(getStringOptionValue());
	else if(option == "-simulate-network-data-imaginary-sample-size") getLatestSimulateNetworkDataTask()->setImaginarySampleSize(getDoubleOptionValue());
	else if(option == "-simulate-network-data-missing-first") getLatestSimulateNetworkDataTask()->setMissingFirst(getStringOptionValue(), getUIntOptionValue());
	else if(option == "-simulate-network-data-missing-last") getLatestSimulateNetworkDataTask()->setMissingLast(getStringOptionValue(), getUIntOptionValue());
	else if(option == "-simulate-network-data-missing-prob") getLatestSimulateNetworkDataTask()->setMissingProb(getStringOptionValue(), getDoubleOptionValue());
	else if(option == "-simulate-network-data-score") getLatestSimulateNetworkDataTask()->setScore(getStringOptionValue());
	else if(option == "-simulate-network-data-score-fix") getLatestSimulateNetworkDataTask()->setScoreFix(getStringOptionValue());

	else if(option == "-output-network") {OutputNetworkTask * ont = new OutputNetworkTask(); outputNetworkTasks.push_back(ont); orderedTasks.push_back(ont);}
	else if(option == "-output-network-name") getLatestOutputNetworkTask()->setTaskName(getStringOptionValue());
	else if(option == "-output-network-network-name") getLatestOutputNetworkTask()->setNetworkName(getStringOptionValue());
	else if(option == "-output-network-file") getLatestOutputNetworkTask()->setFilename(getStringOptionValue());
	else if(option == "-output-network-file1") getLatestOutputNetworkTask()->setFilename(getStringOptionValue());
	else if(option == "-output-network-file2") getLatestOutputNetworkTask()->setFilename2(getStringOptionValue());
	else if(option == "-output-network-equivalent-networks-file") getLatestOutputNetworkTask()->setEquivNetsFilename(getStringOptionValue());
	else if(option == "-output-network-igraph-file-prefix") getLatestOutputNetworkTask()->setFilenamePrefix(getStringOptionValue());
	else if(option == "-output-network-node-data-file-prefix") getLatestOutputNetworkTask()->setNodeDataFilePrefix(getStringOptionValue());
	else if(option == "-output-network-node-data-bed-file") getLatestOutputNetworkTask()->setOutputBedFile(true);
	else if(option == "-output-network-node-data-start-indiv") getLatestOutputNetworkTask()->setStartIndivNo(getUIntOptionValue());
	else if(option == "-output-network-node-data-end-indiv") getLatestOutputNetworkTask()->setEndIndivNo(getUIntOptionValue());
	else if(option == "-output-network-node-data-job") {
		getLatestOutputNetworkTask()->setJobNo(getUIntOptionValue());
		getLatestOutputNetworkTask()->setJobTotal(getUIntOptionValue());
	}

	else if(option == "-output-priors") {OutputPriorsTask * oprt = new OutputPriorsTask(); outputPriorsTasks.push_back(oprt); orderedTasks.push_back(oprt);}
	else if(option == "-output-priors-name") getLatestOutputPriorsTask()->setTaskName(getStringOptionValue());
	else if(option == "-output-priors-network-name") getLatestOutputPriorsTask()->setNetworkName(getStringOptionValue());
	else if(option == "-output-priors-file") getLatestOutputPriorsTask()->setFilename(getStringOptionValue());
	
	else if(option == "-output-posteriors") {OutputPosteriorsTask * opot = new OutputPosteriorsTask(); outputPosteriorsTasks.push_back(opot); orderedTasks.push_back(opot);}
	else if(option == "-output-posteriors-name") getLatestOutputPosteriorsTask()->setTaskName(getStringOptionValue());
	else if(option == "-output-posteriors-network-name") getLatestOutputPosteriorsTask()->setNetworkName(getStringOptionValue());
	else if(option == "-output-posteriors-file") getLatestOutputPosteriorsTask()->setFilename(getStringOptionValue());
	
	else if(option == "-impute-network-data") {ImputeNetworkDataTask * indt = new ImputeNetworkDataTask(); imputeNetworkDataTasks.push_back(indt); orderedTasks.push_back(indt);}
	else if(option == "-impute-network-data-name") getLatestImputeNetworkDataTask()->setTaskName(getStringOptionValue());
	else if(option == "-impute-network-data-network-name") getLatestImputeNetworkDataTask()->setNetworkName(getStringOptionValue());
	else if(option == "-impute-network-data-max-missing") getLatestImputeNetworkDataTask()->setMaxMissing(getUIntOptionValue());
	else if(option == "-impute-network-data-min-non-missing-edges") getLatestImputeNetworkDataTask()->setMinNonMissingEdges(getDoubleOptionValue());
	else if(option == "-impute-network-data-complete-training") getLatestImputeNetworkDataTask()->setUseCompleteData();
	else if(option == "-impute-network-data-random-training") getLatestImputeNetworkDataTask()->setUseRandomData();
	else if(option == "-impute-network-data-use-prev-network") getLatestImputeNetworkDataTask()->setUsePrevNetwork();
	else if(option == "-impute-network-data-random-restarts") getLatestImputeNetworkDataTask()->setRandomRestarts(getUIntOptionValue());
	else if(option == "-impute-network-data-jitter-restarts") getLatestImputeNetworkDataTask()->setJitterRestarts(getUIntOptionValue());
	else if(option == "-impute-network-data-start-indiv") getLatestImputeNetworkDataTask()->setStartIndivNo(getUIntOptionValue());
	else if(option == "-impute-network-data-end-indiv") getLatestImputeNetworkDataTask()->setEndIndivNo(getUIntOptionValue());
	else if(option == "-impute-network-data-job") {
		getLatestImputeNetworkDataTask()->setJobNo(getUIntOptionValue());
		getLatestImputeNetworkDataTask()->setJobTotal(getUIntOptionValue());
	}
	else if(option == "-impute-estimate-recall-precision") {ImputeEstimateRecallPrecisionTask * iert = new ImputeEstimateRecallPrecisionTask(); imputeEstimateRecallPrecisionTasks.push_back(iert); orderedTasks.push_back(iert);}
	else if(option == "-impute-estimate-recall-precision-name") getLatestImputeEstimateRecallPrecisionTask()->setTaskName(getStringOptionValue());
	else if(option == "-impute-estimate-recall-precision-network-name") getLatestImputeEstimateRecallPrecisionTask()->setNetworkName(getStringOptionValue());
	else if(option == "-impute-estimate-recall-precision-random-restarts") getLatestImputeEstimateRecallPrecisionTask()->setRandomRestarts(getUIntOptionValue());
	else if(option == "-impute-estimate-recall-precision-jitter-restarts") getLatestImputeEstimateRecallPrecisionTask()->setJitterRestarts(getUIntOptionValue());
	else if(option == "-impute-estimate-recall-precision-skip-imputation") getLatestImputeEstimateRecallPrecisionTask()->setDoImps(false);
	else if(option == "-impute-estimate-recall-precision-iterations") getLatestImputeEstimateRecallPrecisionTask()->setIterations(getUIntOptionValue());
	else if(option == "-impute-estimate-recall-precision-output-sim-data-file") getLatestImputeEstimateRecallPrecisionTask()->setSimDataFilename(getStringOptionValue());
	else if(option == "-impute-estimate-recall-precision-output-sim-network-file") getLatestImputeEstimateRecallPrecisionTask()->setSimNetworkFilename(getStringOptionValue());
	else if(option == "-calculate-recall-precision") {CalculateRecallPrecisionTask * crpt = new CalculateRecallPrecisionTask(); calculateRecallPrecisionTasks.push_back(crpt); orderedTasks.push_back(crpt);}
	else if(option == "-calculate-recall-precision-network-name") getLatestCalculateRecallPrecisionTask()->setNetworkName(getStringOptionValue());
	else if(option == "-calculate-recall-precision-true-network-name") getLatestCalculateRecallPrecisionTask()->setTrueNetworkName(getStringOptionValue());
	else if(option == "-calculate-recall-precision-file") getLatestCalculateRecallPrecisionTask()->setFilename(getStringOptionValue());
	else
	{
		validOption3 = false;
	};

	return validOption1 || validOption2 || validOption3;
};

//! Adds task name to list.
void Analysis::addTaskName(const string & taskName)
{
	if(taskNames.find(taskName) != taskNames.end())
	{
		outErr("Task name \""); outErr(taskName); outErr("\" is repeated, task names must be unique!");		
		exitErr("");
	};

	taskNames.insert(taskName);
};

//! Do all of the tasks.
void Analysis::doTasks()
{
	int taskNo = 1;
	string taskName;

	for(list<Task *>::const_iterator t = orderedTasks.begin(); t != orderedTasks.end(); ++t, ++taskNo)
	{
		out("--------------------------------------------------\n");
		
		if((*t)->getTaskName() == "") (*t)->setTaskName("Task-"+toString(taskNo)); //set default task name		
		addTaskName((*t)->getTaskName()); //check task name is unique

		(*t)->outputTaskHeader();
		(*t)->setTaskCentral(taskCentral); //allow Tasks to link together to access stuff
		(*t)->initialiseTask();
		(*t)->doTask();

		(*t)->outputTaskDetails();
		out("--------------------------------------------------\n");
	};

};
