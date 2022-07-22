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


/*! \file Task.cpp
    \brief This file contains the methods for all the tasks.
    
*/
#include <iostream>
#include <sstream>

#include <errno.h>

using namespace std; // initiates the "std" or "standard" namespace

#include "Task.h"
#include "main.h"
#include "Network.h"
#include "Utils.h"
#include "cdflib.h"
#include <math.h>
#include <stack>
#include <vector>

//! Adds a network.
void TaskCentral::addNetwork(const string & name, Network *network)
{
	map<string, Network *>::const_iterator nt = allNetworks.find(name);

	if(nt != allNetworks.end())
	{
		string mess = "Attempt to add network \"" + name + "\" twice!\n";
		exitErr(mess);
	};
	
	allNetworks[name] = network;
	orderedNetworks.push_back(name);
};

//! Removes a network.
void TaskCentral::removeNetwork(const string & name)
{
	map<string, Network *>::const_iterator nt = allNetworks.find(name);

	if(nt == allNetworks.end())
	{
		string mess = "Attempt to remove network \"" + name + "\" that does not exist!\n";
		exitErr(mess);
	};
	
	allNetworks.erase(name);
	orderedNetworks.remove(name);
};

//! Gets network by name.
Network * TaskCentral::getNetwork(const string & name)
{
	map<string, Network *>::const_iterator nt = allNetworks.find(name);

	if(nt == allNetworks.end())
	{
		string mess = "Attempt to access network \"" + name + "\", which does not exist!\n";
		exitErr(mess);
	};

	return nt->second;
};

//! Gets the latest network.
Network * TaskCentral::getLatestNetwork(string & networkName)
{
	list<string>::const_reverse_iterator nt = orderedNetworks.rbegin();

	if(nt == orderedNetworks.rend())
	{
		Network * noEdgeNet = getNoEdgeNetwork();
		networkName = "defaultNetwork";
		
		addNetwork(networkName, noEdgeNet);

		unsigned int totalDisNodes, totalFactorNodes, totalCtsNodes, totalNoDataNodes, totalEdges;
		noEdgeNet->getNumberNodesAndEdges(totalDisNodes, totalFactorNodes, totalCtsNodes, totalNoDataNodes, totalEdges);
		out("--------------------------------------------------\n");
		out("Loading "); out(networkName); out(" network\n");	
		out("Network type: "); out(noEdgeNet->getNetworkType()); out("\n");
		out("Network score type: "); out(noEdgeNet->getScoreTypeName()); out("\n");
		out("Total number of nodes: "); out(totalNoDataNodes + totalDisNodes + totalFactorNodes + totalCtsNodes); out(" (Discrete: ");  out(totalDisNodes); out(" | Factor: ");  out(totalFactorNodes); out(" | Continuous: "); out(totalCtsNodes);
			if(totalNoDataNodes != 0) {out(" | No data: "); out(totalNoDataNodes);} out(")"); out("\n");
		out("Total number of edges: "); out(totalEdges); out("\n");
		out("Network Structure: "); out(noEdgeNet->getNetworkString()); out("\n");
		if(noEdgeNet->hasLoop()) {out("The network contains a loop!\n");};
		if(noEdgeNet->getScoreTypeName() == "bayes" || noEdgeNet->getNetworkType() == "deal") {out("Imaginary sample size: "); out(noEdgeNet->getImaginarySampleSize()); out("\n");};
		addNetworkMissingData(noEdgeNet->updateNetworkMissingData());
		pair<unsigned int, unsigned int> missingNotMissing = noEdgeNet->getNoMissingNotMissing();
		out("Total data at each node: "); out(missingNotMissing.second); out("\n");
		out("Missing data at each node: "); out(missingNotMissing.first); out("\n");	
		out("--------------------------------------------------\n");

		return noEdgeNet;
	};

	networkName = *nt;

	return getNetwork(networkName);
};

//! Creates a network with all nodes and no edges.
Network * TaskCentral::getNoEdgeNetwork()
{
	unsigned int noObs = allNodeData->getAmountOfData();
	
	if(noObs==0) exitErr("No data loaded to create a network!");

	NetworkBNLearn * defaultNetwork = new NetworkBNLearn(allNodeData, defaultScoreType, 0);

	map<string, unsigned int> nameNumber = allNodeData->getNameNumber();

	for(map<string, unsigned int>::const_iterator nn = nameNumber.begin(); nn != nameNumber.end(); ++nn)
	{
		defaultNetwork->addNodeInit(nn->first);
	};

	double imaginarySampleSize = 12; 
	defaultNetwork->setImaginarySampleSize(imaginarySampleSize);

	set<pair<unsigned int, unsigned int> > allowedEdgeTypesFileNos; // fileNo1, fileNo2
	set<pair<unsigned int, unsigned int> > notAllowedEdgeTypesFileNos;
	
	defaultNetwork->setUpAllowedEdgeTypes(allowedEdgeTypesFileNos, notAllowedEdgeTypesFileNos);

	return defaultNetwork;
};

//! Creates a network with nodes given by a network with updated edges.
void TaskCentral::setupRestartNetwork(Network * restartNetwork, Network * network, const string & preNodeName)
{

	unsigned int noObs = allNodeData->getAmountOfData(); //number of individuals
	
	if(noObs==0) exitErr("No data loaded to create a network!");

	map<string, unsigned int> nameNumber = allNodeData->getNameNumber();
	
	//add nodes
	unsigned int nd = 1;
	unsigned int noCurNetworkNodes = 0;
	unsigned int noNodes = network->getNoNodes();
	unsigned int noAllDataNodes = allNodeData->getNoNodeData();
	string nodeName, displayNodeName;

	while(nd <= noAllDataNodes && noCurNetworkNodes < noNodes)
	{
			if(network->nodeExists(nd) && !network->isBlackNode(nd))
			{
				nodeName = network->getNetworkNode(nd)->getName();
				restartNetwork->addNodeInit(nodeName);

				displayNodeName = network->getNetworkNode(nd)->getDisplayName();
				restartNetwork->getNetworkNode(nd)->setDisplayName(displayNodeName);

				++noCurNetworkNodes;
			};
			++nd;
	};

	//copy missing data
	restartNetwork->setNetworkMissingData(network);

	double imaginarySampleSize = network->getImaginarySampleSize();
	restartNetwork->setImaginarySampleSize(imaginarySampleSize);

	restartNetwork->updateBlackWhite(network); //also updates allowed edges
	//if(preNodeName == "") 
	//else restartNetwork->updateBlackWhiteDifferentData(network, preNodeName); //restarat network always has the same name nodes

	//copy edge costs
	map<pair<unsigned int, unsigned int>, double> costEdges = network->getCostEdges();
	map<pair<unsigned int, unsigned int>, double> costEdgeTypes = network->getCostEdgeTypes();;
	restartNetwork->setCostEdges(costEdges, costEdgeTypes);	
};


//! Creates a network with nodes given by a network with updated edges.
void TaskCentral::updateRestartNetwork(Network * restartNetwork, Network * network, const unsigned int & minEdges, const unsigned int & maxEdges, const bool & doJitter)
{
	restartNetwork->removeAllEdges();

	//add nodes
	unsigned int nd = 1;
	unsigned int noCurNetworkNodes = 0;
	unsigned int noNodes = network->getNoNodes();
	unsigned int noAllDataNodes = allNodeData->getNoNodeData();

	//choose number of edges to add or change
	unsigned int noEdges = rand() % (maxEdges - minEdges + 1) + minEdges;

	Node * origNode;
	Node * copyingNode;

	if(doJitter)
	{
		//add the same edges as the prev network		
		nd = 1;
		noCurNetworkNodes = 0;

		while(nd <= noAllDataNodes && noCurNetworkNodes < noNodes)
		{
			if(network->nodeExists(nd))
			{
				origNode = network->getNetworkNode(nd);
				copyingNode = restartNetwork->getNetworkNode(nd);
				origNode->copyParents(copyingNode, restartNetwork);
				++noCurNetworkNodes;
			};
			++nd;
		};

		//update random edges
		for(unsigned int r = 1; r <= noEdges; ++r) restartNetwork->changeRandomEdge();

	}
	else
	{
		//add edges
		for(unsigned int r = 1; r <= noEdges; ++r) restartNetwork->addRandomEdge();
	};

};

//! Copies edges of network to another network.
void TaskCentral::copyNetworkEdges(Network * copyNetwork, Network * network)
{
	copyNetwork->removeAllEdges();

	//add nodes
	unsigned int nd = 1;
	unsigned int noCurNetworkNodes = 0;
	unsigned int noNodes = network->getNoNodes();
	unsigned int noAllDataNodes = allNodeData->getNoNodeData();

	Node * origNode;
	Node * copyingNode;

	//add the same edges as the prev network		
	nd = 1;
	noCurNetworkNodes = 0;

	while(nd <= noAllDataNodes && noCurNetworkNodes < noNodes)
	{
		if(network->nodeExists(nd))
		{
			origNode = network->getNetworkNode(nd);
			copyingNode = copyNetwork->getNetworkNode(nd);
			origNode->copyParents(copyingNode, copyNetwork);
			++noCurNetworkNodes;
		};
		++nd;
	};	

};

//! Returns the file no from the data name
unsigned int TaskCentral::getFileNoFromDataName(const string & fileDataName)
{
	map<string, unsigned int>::const_iterator fn = dataNameFileNo.find(fileDataName);

	if(fn != dataNameFileNo.end()) return fn->second;

	outErr("Cannot find data input with task name "); outErr(fileDataName); outErr("!");
	exitErr("");

	return 0; //never called
};

//! Sets up the start and end SNP nos based on job n of N type things.
void TaskCentral::setStartAndEndIndivsBasedOnJobNo(unsigned int & jobNo, unsigned int & jobTotal, unsigned int & startIndivNo, unsigned int & endIndivNo)
{
	double noIndivs = (double)(getAmountOfData());

	double stepSize = noIndivs/(double)(jobTotal);

	//check job options are acceptable
	if(jobNo > jobTotal)
	{
		outErr("The job number, "); outErr(jobNo); outErr(", may not be larger than the job total, "); outErr(jobTotal);
		
		exit(1);
	};

	//round up start Indiv to be multiple of stepSize, then plus 1
	startIndivNo = (unsigned int)(stepSize*(jobNo-1) + 1);

	endIndivNo = (unsigned int)(stepSize*jobNo);

	if(jobNo == jobTotal) endIndivNo = (unsigned int)noIndivs; //just to make sure last indiv is imputed
};

//! Initialises the calc. posterior task.
void CalculatePosteriorTask::initialiseTask()
{
	
};

//! Sets up SNP data from the .bim file
void DescriptionOfSNPs::setUpSNPDesciptionData(string & filename)
{
	//try and find the binary map file, .bim, and read in data
	unsigned int length = filename.length();
	string mapFilename = filename.substr(0, length-4) + ".bim";

	ifstream readMapFile;
	readMapFile.open(mapFilename.c_str());
	if(!readMapFile.is_open())
	{
		string mess = "Cannot read bim file: " + mapFilename + "!\n";
		exitErr(mess);
	};

	string snpIdentifier, geneticDistance;
	string prevSnpIdentifier = "";
	unsigned int chromosome, basePairPosition;	
	string alleleName1, alleleName2;
	unsigned int snpID = 1;

	//loop thro' subjects and store the cases
	do{
		
		readMapFile >> chromosome >> snpIdentifier >> geneticDistance >> basePairPosition >> alleleName1 >> alleleName2;
		
		if(snpIdentifier != prevSnpIdentifier)
		{				
				snpInfos[snpID] = new SNPInfo(chromosome, snpIdentifier, basePairPosition);
				
				snpID++;
		};
		
		prevSnpIdentifier = snpIdentifier;				
	}while(!readMapFile.eof());

	readMapFile.close();

	//set the iterator to the first SNP
	snpInfo = snpInfos.begin();
};

//! Outputs info for the SNP data.
void DescriptionOfSNPs::outputSNPInfo()
{
	out("Total number of SNPs: "); out(snpInfos.size()); out("\n");
};

//! Outputs the task header for the input data task.
void InputDataTask::outputTaskHeader()
{
	outputTaskName();
	out("Loading data\n");
};

//! Outputs task details for the input data task.
void InputDataTask::outputTaskDetails()
{
	
	if(taskDataType == 1)
	{
		out("Continuous data file: "); out(filename); out("\n");
	}
	else if(taskDataType == 2)
	{
		out("Discrete data file: "); out(filename); out("\n");
	}
	else if(taskDataType == 3)
	{
		out("SNP binary data file: "); out(filename); out("\n");
		out("SNP data treated as discrete data\n");
		descSNPs.outputSNPInfo();
		out("Total number of subjects: "); out(totalNoSubjects); out("\n");
	}
	else if(taskDataType == 4)
	{
		out("SNP binary data file: "); out(filename); out("\n");
		out("SNP data treated as continuous data (0, 1 or 2)\n");
		descSNPs.outputSNPInfo();
		out("Total number of subjects: "); out(totalNoSubjects); out("\n");
	}
	else if(taskDataType == 5)
	{
		out("Discrete data file: "); out(filename); out("\n");
		out("Data treated as factors\n");
	}
	else if(taskDataType == 6)
	{
		out("SNP binary data file: "); out(filename); out("\n");
		out("SNP data treated as factors\n");
		descSNPs.outputSNPInfo();
		out("Total number of subjects: "); out(totalNoSubjects); out("\n");
	}
	else if(taskDataType == 7)
	{
		out("Continuous SNP data file: "); out(filename); out("\n");
		out("SNP data treated as continuous data\n");
	}
	else if(taskDataType == 8)
	{
		out("Discrete SNP data file: "); out(filename); out("\n");
		out("SNP data treated as discrete data\n");
	}
	else if(taskDataType == 9)
	{
		out("Discrete SNP data file: "); out(filename); out("\n");
		out("SNP data treated as factors\n");
	}

	out("Number of ID columns: "); out(numberOfIDs); out("\n");

	if(includeVariables.empty() && excludeVariables.empty())
	{
		if(totalNodeDatas == 1)
		{
			out("Including the 1 and only variable in analysis\n");
		}
		else
		{
			out("Including (all) "); out(totalNodeDatas); out(" variables in analysis\n");
		};
	}
	else 
	{
		if(!includeVariables.empty()) out("Include variables file: ");
		else out("Exclude variables file: ");
		out(variableFilename); out("\n");
		if(totalNodeDatas == 1)
		{
			out("Including 1 variable in analysis\n");
		}
		else
		{
			out("Including "); out(totalNodeDatas); out(" variables in analysis\n");
		};
	};

	out("Each variable has "); out(taskCentral->getAmountOfData()); out(" data entries\n");

	if(taskDataType == 1 || taskDataType == 7)
	{
		if(!ctsMissingValueSet) out("Missing value: not set\n");
		else
		{
			if(ctsMissingValueStr != "") { out("Missing value: "); out(ctsMissingValueStr); out("\n"); }
			else { out("Missing value: "); out(ctsMissingValue); out("\n"); };
		};
	}
	else if(taskDataType == 2 || taskDataType == 5 || taskDataType == 8 || taskDataType == 9) {out("Missing value: "); out(discreteMissingValue); out("\n");};

};
	
//! Inputs data.
void InputDataTask::doTask()
{
	fileno = taskCentral->getFileNo();

	if(filename == "")
	{
		string mess = "Data file not set!\n Set with -input-data-file myfile.dat";		
		exitErr(mess);
	}
	else if(taskDataType == 0)
	{
		string mess = "Data file type not set!\n Choose: -input-data-cts, -input-data-discrete or -input-data-snp";		
		exitErr(mess);
	}
	else if(taskDataType == 1) loadCtsData();
	else if(taskDataType == 2) loadDiscreteData();
	else if(taskDataType == 5) loadDiscreteDataAsFactors();
	else if(taskDataType == 3 || taskDataType == 4) loadSNPData();
	else if(taskDataType == 6) loadSNPDataAsFactors();
	else if(taskDataType == 7) loadCtsData(true);
	else if(taskDataType == 8) loadDiscreteData(true);
	else if(taskDataType == 9) loadDiscreteDataAsFactors(true);

	taskCentral->updateDataNames(); //update node names to remove ":"s from names

	taskCentral->nextFileNo();

	taskCentral->setFileNoFromDataName(getTaskName(), fileno);
};
	
//! Initialises task for inputtting data.
void InputDataTask::initialiseTask()
{
				
};

//! Reads in variables to include.
void InputDataTask::setIncludeVariables(const string & filename)
{
	ifstream readData(filename.c_str());
	
	if(!readData.is_open())
	{
		string mess = "Cannot find include variables file: " + filename + "!";
		exitErr(mess);
	};

	string aVariable;

	readData >> aVariable;

	while(!readData.eof()){

		includeVariables.insert(aVariable);

		readData >> aVariable;
	};

	readData.close();

	checkIncludeExcludeSet();

	variableFilename = filename;
};

//! Checks if both include and exclude variables are set.
void InputDataTask::checkIncludeExcludeSet()
{
	if(!includeVariables.empty() && !excludeVariables.empty())
	{
		string mess = "You may not set include and exclude variables for one data file!";
		exitErr(mess);
	};
};

//! Reads in variables to exclude.
void InputDataTask::setExcludeVariables(const string & filename)
{
	ifstream readData(filename.c_str());
	
	if(!readData.is_open())
	{
		string mess = "Cannot find exclude variables file: " + filename + "!";
		exitErr(mess);
	};

	string aVariable;

	readData >> aVariable;

	while(!readData.eof()){

		excludeVariables.insert(aVariable);

		readData >> aVariable;
	};

	readData.close();

	checkIncludeExcludeSet();

	variableFilename = filename;
};

//! Checks if data is to be included as a (potential) node.
bool InputDataTask::isIncluded(const string & variable)
{
	if(includeVariables.empty() && excludeVariables.empty()) return true; //include everything by default
	else if(!excludeVariables.empty())
	{
		set<string>::const_iterator cv = excludeVariables.find(variable);

		return cv == excludeVariables.end();
	}
	else
	{
		set<string>::const_iterator cv = includeVariables.find(variable);

		return cv != includeVariables.end();
	};
};

//! Loads data from SNP data.
void InputDataTask::loadSNPData()
{
	//check no of IDs
	if(numberOfIDs != 2)
	{
		outErr("The number of column IDs must be set to 2 for (.bed) SNP data, but it is set to "); outErr(numberOfIDs); outErr("!\n");
		outErr("If you are using SNP data in a text file (0, 1, 2) then use option -input-data-cts-snp2\n");
		exit(1);
	};

	handleFamFile(); //sets the total number of subjects, and check IDs
	handleBimFile(); //gets SNP names and number of SNPs

	openBinaryFileFirst();

	//loop thro' SNPs and SNP info
	unsigned int snpNo = 1;
	unsigned int totalSNPs = descSNPs.getNoSNPs();
	string snpName; //rs number
	unsigned char snp;
	bool snpIsIncluded;
	DiscreteData * nothingDisData = new DiscreteData(snpName, fileno); //just so it compiles
	DiscreteData * aDiscreteData = nothingDisData;
	CtsData * nothingCtsData = new CtsData(snpName, fileno);
	CtsData * aCtsData = nothingCtsData;

	while(snpNo <= totalSNPs)
	{
		snpName = descSNPs.getSNPInfo()->name;
		snpIsIncluded = isIncluded(snpName);

		//if SNP is to be included then create SNP Node and add data to it
		if(snpIsIncluded)
		{
			if(taskDataType == 3) //discrete data for taskDataType == 3
			{
				aDiscreteData = new DiscreteData(snpName, fileno);
				aDiscreteData->isSNPData = true;
				taskCentral->addNodeData(aDiscreteData); //add to general list for access later
			}
			else if(taskDataType == 4) //discrete data for taskDataType == 4
			{
				aCtsData = new CtsData(snpName, fileno);	
				aCtsData->isSNPData = true;
				taskCentral->addNodeData(aCtsData);
			};
		
			++totalNodeDatas;
		};

		bitCount = 9; //ensure new byte is read in for the start of a new SNP

		string noMinorAllelesFactor; //factor for discrete data
		bool missing;
		double ctsData;

		//add SNP data, need to scroll thro' it even if not included
		for(unsigned int indiv = 1; indiv <= totalNoSubjects; ++indiv)
		{
			missing = false;
			snp = getNextNoOfMinorAlleles();

			//add data
			if(snpIsIncluded)
			{
				if(taskDataType == 3)
				{
					if(snp == 0) noMinorAllelesFactor = "0";
					else if(snp == 1) noMinorAllelesFactor = "1";
					else if(snp == 2) noMinorAllelesFactor = "2";
					else {noMinorAllelesFactor = discreteMissingValue; missing = true;};

					aDiscreteData->addData(noMinorAllelesFactor, missing);					
				}
				else if(taskDataType == 4)
				{
					if(snp == 3) {ctsData = ctsMissingValue; missing = true;}
					else ctsData = (double)snp;

					aCtsData->addData(ctsData, missing);
				};
			};
			
		};

		descSNPs.advanceIterator();
		++snpNo;
	};
	//end of loop

	delete nothingDisData;
	delete nothingCtsData;
	readSNPData.close();
};

//! Loads data from SNP data.
void InputDataTask::loadSNPDataAsFactors()
{
	//check no of IDs
	if(numberOfIDs != 2)
	{
		outErr("The number of column IDs must be set to 2 for (.bed) SNP data, but it is set to "); outErr(numberOfIDs); outErr("!\n");
		outErr("If you are using SNP data in a text file (0, 1, 2) then use option -input-data-factor-snp2\n");
		exit(1);
	};

	handleFamFile(); //sets the total number of subjects, and check IDs
	handleBimFile(); //gets SNP names and number of SNPs

	openBinaryFileFirst();

	//loop thro' SNPs and SNP info
	unsigned int snpNo = 1;
	unsigned int totalSNPs = descSNPs.getNoSNPs();
	string snpName; //rs number
	unsigned char snp;
	bool snpIsIncluded;
	CtsData * nothingCtsData = new CtsData(snpName, fileno); //just so it compiles
	CtsData * aCtsFactorData1 = nothingCtsData;
	CtsData * aCtsFactorData2 = nothingCtsData;
	string factorNodeName;
	list<string> factorNames;

	while(snpNo <= totalSNPs)
	{
		snpName = descSNPs.getSNPInfo()->name;
		snpIsIncluded = isIncluded(snpName);

		//if SNP is to be included then create SNP Node and add data to it
		if(snpIsIncluded)
		{
			//add head factor, for hetreozguous, 1 minor allele
			aCtsFactorData1 = new CtsData(snpName, fileno);	
			aCtsFactorData1->isSNPData = true;
			aCtsFactorData1->isFactorNode = true;
			aCtsFactorData1->nodeName = snpName;
			aCtsFactorData1->factorName = "1";			
			taskCentral->addNodeData(aCtsFactorData1);

			factorNodeName = snpName + " 2";
			aCtsFactorData2 = new CtsData(factorNodeName, fileno);
			aCtsFactorData2->isSNPData = true;
			aCtsFactorData2->isFactorChildNode = true;
			aCtsFactorData2->nodeName = snpName;
			aCtsFactorData2->factorName = "2";			
			taskCentral->addNodeData(aCtsFactorData2);			
			
			factorNames.clear();
			factorNames.push_back(factorNodeName);
			taskCentral->getAllNodeData()->addFactorDataGroup(snpName, factorNames); //add list of other factors other than head factor, i.e. "2" only

			totalNodeDatas += 2;			
		};

		bitCount = 9; //ensure new byte is read in for the start of a new SNP

		string noMinorAllelesFactor; //factor for discrete data

		//add SNP data, need to scroll thro' it even if not included
		for(unsigned int indiv = 1; indiv <= totalNoSubjects; ++indiv)
		{
		
			snp = getNextNoOfMinorAlleles();

			//add data
			if(snpIsIncluded)
			{
			
				if(snp == 0)
				{
					aCtsFactorData1->addData(0, false);
					aCtsFactorData2->addData(0, false);
				}
				else if(snp == 1)
				{
					aCtsFactorData1->addData(1, false);
					aCtsFactorData2->addData(0, false);
				}
				else if(snp == 2)
				{
					aCtsFactorData1->addData(0, false);
					aCtsFactorData2->addData(1, false);
				} 
				else
				{					
					aCtsFactorData1->addData(ctsMissingValue, true);
					aCtsFactorData2->addData(ctsMissingValue, true);
				};
				
			};
			
		};

		descSNPs.advanceIterator();
		++snpNo;
	};
	//end of loop

	delete nothingCtsData;
	readSNPData.close();
};

//! Open the binary file for the first time
void InputDataTask::openBinaryFileFirst()
{
	//try and find the binary pedigree file, .bed, and read in data for the first window
	readSNPData.open(filename.c_str(), ios::binary);
	
	if(!readSNPData.is_open())
	{
		string mess = "Cannot read binary pedigree file: " + filename + "!";
		exitErr(mess);
	};

	char buffer[3];
	readSNPData.read(buffer, 3);

	//check the plink magic numbers for the file type
	//3rd number indicates format of genotype data, 1 => subjects x SNPs, 0 => SNPs x subjects
	unsigned int magicNumber1 = buffer[0];
	unsigned int magicNumber2 = buffer[1];

	if(magicNumber1 != 108 || magicNumber2 != 27)
	{
		readSNPData.close();	

		outErr("Detected an old version .bed file!");
		exitErr("Please use PLINK to update the .bed file.");	
	};

	//determine binary file type
	unsigned int mode = buffer[2];
	if(mode == 0)
	{
		readSNPData.close();		
		
		outErr("The binary pedigree file must be in SNP-major mode!");
		exitErr("Please use PLINK to update the .bed file.");
	};
};

//! Gets no of binary alleles from .bed file.
unsigned char InputDataTask::getNextNoOfMinorAlleles()
{
	int allele1, allele2;
	unsigned char noMinorAlleles = 0;

	//read in the next piece of data
	if(bitCount == 9)
	{
		
		readSNPData.read(buffer, 1);
		if(readSNPData.eof())
		{			
			string mess = "Error: reached end of binary SNP file!";
			exitErr(mess);
		};
			
		aBit = buffer[0];

		bitCount = 1;
	};

	allele1 = aBit & one; //read the least significant bit				
	aBit = aBit >> 1; //shift bits to the right
	allele2 = aBit & one; //read the new least significant bit				
	aBit = aBit >> 1; //shift bits to the right for next time

	bitCount += 2;	

	//if genotype is encoded 1/0 then the genotype is missing so do not add it
	if(allele1 == 1 && allele2 == 1)
	{	
		noMinorAlleles = 0;
	}
	else if(allele1 == 0 && allele2 == 1)
	{	
		noMinorAlleles = 1;
	}
	else if(allele1 == 0 && allele2 == 0)
	{	
		noMinorAlleles = 2;
	}
	else
		noMinorAlleles = 3; //denotes missing genotype

	return noMinorAlleles;
};

//! Loads family data from .fam SNP data.
void InputDataTask::handleFamFile()
{
	//try and find the family file and read in data
	unsigned int length = filename.length();
	string famFilename = filename.substr(0, length-4) + ".fam";

	ifstream readFamilyFile;
	readFamilyFile.open(famFilename.c_str());
	if(!readFamilyFile.is_open())
	{
		string mess = "Cannot read family file: " + famFilename + "!";
		exitErr(mess);
	};

	string famID, indivID, FatherId, MotherID, sexID, famIndivID;
	string prevFamIndivID = "";
	double phenoType;
	unsigned int noCases = 0;
	unsigned int subjectNo = 0;
	unsigned int noMissingPhenos = 0;
	list<string> dataIDs;

	//loop thro' subjects 
	do{		
		readFamilyFile >> famID >> indivID >> FatherId >> MotherID >> sexID >> phenoType;
		famIndivID = famID + " " + indivID;
	
		//do not duplicate the last row
		if(famIndivID != prevFamIndivID) 
		{
			dataIDs.push_back(famID);
			dataIDs.push_back(indivID);		
			subjectNo++;			
		};

		prevFamIndivID = famIndivID;
	}while(!readFamilyFile.eof());

	readFamilyFile.close();

	taskCentral->checkIDs(famFilename, dataIDs);

	//set name of IDs for possible data output
	list<string> nameOfIDs;
	nameOfIDs.push_back("FAMID");
	nameOfIDs.push_back("ID");
	taskCentral->setNameOfIDs(nameOfIDs);

	totalNoSubjects = subjectNo;
};

//! Loads SNP data infomation from .bim SNP data.
void InputDataTask::handleBimFile()
{
	descSNPs.setUpSNPDesciptionData(filename);
};

//! Returns variable names from header of data file.
void InputDataTask::getVariableNames(list<string> & variableNames, const bool & secondLine)
{
	ifstream readData(filename.c_str());
	
	if(!readData.is_open())
	{
		string mess = "Cannot read data file: " + filename + "!";
		exitErr(mess);
	};

	string line;
	unsigned int pos = 0;
	string variableName;
	bool prevSpace = true;

	//get header
	getline(readData, line);
	if(secondLine) getline(readData, line);

	unsigned int length = line.length();

	if(length < 1)
	{
		string mess = "There is a problem with the header of the data file: " + filename + "!";
		exitErr(mess);
	};

	//count names then just use <<
	unsigned int noWords = 1;

	do{
		if((!isCSV && (line.substr(pos, 1) == " " || line.substr(pos, 1) == "\t")) ||
		   (isCSV && line.substr(pos, 1) == ","))
		{
			if(!prevSpace)
			{				
				noWords++;				
			};

			prevSpace = true;		
		}
		else
			prevSpace = false;
			
		++pos;		
	}while(pos < length);

	readData.close();

	if(numberOfIDs >= noWords)
	{
		if(numberOfIDs==1) {outErr("There is 1 column ID set, but only ");}
		else {outErr("There are "); outErr(numberOfIDs); outErr(" column IDs set, but only ");}

		if(noWords==1) {outErr("1 column in the data file "); outErr(filename); outErr("!\n\n");}
		else {outErr(noWords); outErr(" columns in the data file "); outErr(filename); outErr("!\n\n");}

		exit(1);
	};

	//add words
	readData.open(filename.c_str());

	for(unsigned int w = 1; w <= noWords; ++w) 
	{
		getStringInputData(readData, variableName);
		variableNames.push_back(prefixName+variableName);
	};

	//pop off the id names from the front
	for(unsigned int idNo = 1; idNo <= numberOfIDs; ++idNo)
	{
		variableNames.pop_front();
	};

	readData.close();
};

//! Checks input from file and deals with nans and NAs.
void checkInput(ifstream & readFile, bool & missing, string & badNumber)
{
	if(readFile.fail())
	{
		readFile.clear();
		readFile >> badNumber;
		missing = true;
	}
	else
	{
		missing = false;
	};
};

//! Loads string input data from file.
void InputDataTask::getStringInputData(ifstream & readData, string & aString)
{
	if(!isCSV)
	{
		readData >> aString;
	}
	else
	{		
		
		if(!getline(lineStream, aString, comma))
		{
			getline(readData, line);
			line = line.substr(0, line.length() - 1);// remove dodgy return 
			lineStream.str(line);
			lineStream.clear();
			
			getline(lineStream, aString, comma);
		};
	};
};

//! Loads number input data from file.
void InputDataTask::getNumberInputData(ifstream & readData, double & aNum, bool & missing, string & badNumber)
{
	if(!isCSV)
	{
		readData >> aNum;
		checkInput(readData, missing, badNumber);
	}
	else
	{
		if(!getline(lineStream, aNumber, comma))
		{
			getline(readData, line);
			line = line.substr(0, line.length() - 1); // remove dodgy return 
			lineStream.str(line);
			lineStream.clear();
			
			getline(lineStream, aNumber, comma);
		};

		char * e;
		errno = 0;
		aNum = strtod(aNumber.c_str(), &e);

		if (*e != '\0' ||  // error, we didn't consume the entire string
			errno != 0 )   // error, overflow or underflow
		{
			// fail
			missing = true;
			badNumber = aNumber;
		}
		else
		{
			missing = false;
		};

	};

};

//! Sets missing value for cts data.
void InputDataTask::setCtsMissingValue(const string & cmv)
{ 
	double dbl = 0.0;
	istringstream num(cmv);

	num >> dbl;

	if(!num.fail() && num.eof()) // This second test is important! This makes sure that the entire string was converted to a number
	{
		// success
		ctsMissingValue = dbl;
		ctsMissingValueStr = "";
	}
	else
	{
		// failure
		ctsMissingValueStr = cmv;
	};

	ctsMissingValueSet = true;
};

//! Loads data for cts data.
void InputDataTask::loadCtsData(const bool & isSNPData)
{
	//first read in header and create Nodes
	list<string> variableNames, secondLine;
	CtsData * aCtsData;
	list<CtsData *> ctsDatas;
	list<CtsData *>::iterator cd;
	list<string> dataIDs;

	double aNum;
	string notNumber;
	bool missingData = false;

	getVariableNames(variableNames);
	getVariableNames(secondLine, true);
	if(variableNames.size() != secondLine.size())
	{
		outErr("There are "); outErr(variableNames.size()); outErr(" variables in the header but only ");
		outErr(secondLine.size()); outErr(" variables in the next line for file "); outErr(filename); outErr("!");
		exitErr("");
	};

	secondLine.clear();

	ifstream readData(filename.c_str());

	if(!readData.is_open())
	{
		string mess = "Cannot read data file: " + filename + "!";
		exitErr(mess);
	};

	//set name of IDs for possible data output
	list<string> nameOfIDs;

	//read in IDs, to advance to variable names below
	for(unsigned int idNo = 1; idNo <= numberOfIDs; ++idNo)
	{
		getStringInputData(readData, notNumber);
		nameOfIDs.push_back(notNumber);
	};

	taskCentral->setNameOfIDs(nameOfIDs);

	//create CtsData objects
	list<string>::const_iterator vn = variableNames.begin();
	for( ; vn != variableNames.end(); ++vn)
	{
		//read in header, check if header exists
		getNumberInputData(readData, aNum, missingData, notNumber);

		if(!missingData)
		{
			string mess = "A header must be included in the data file!\n Found a number in the first line of " + filename + "!";
			exitErr(mess);
		};
	
		if(isIncluded(*vn))
		{
			aCtsData = new CtsData(*vn, fileno);
			if(isSNPData) aCtsData->isSNPData = true;

			taskCentral->addNodeData(aCtsData); //add to general list for access later
			ctsDatas.push_back(aCtsData);
			++totalNodeDatas;
		};
	};

	
	//loop thro subjects
	unsigned int subjectNo = 1; 
	cd = ctsDatas.begin();
	vn = variableNames.begin();

	//read in IDs, to advance to data to the right
	for(unsigned int idNo = 1; idNo <= numberOfIDs; ++idNo)
	{
		getStringInputData(readData, notNumber);
		dataIDs.push_back(notNumber);
	};
	
	getNumberInputData(readData, aNum, missingData, notNumber);

	while(!readData.eof())
	{
		missingData = missingData || (ctsMissingValueSet && aNum == ctsMissingValue);

		if(isIncluded(*vn))
		{
			//find corresponding cts Node
			while((*cd)->name != *vn)
			{
				++cd;
				if(cd == ctsDatas.end()) cd = ctsDatas.begin();
			};

			(*cd)->addData(aNum, missingData);
		};
	
		//update variable
		++vn;
		if(vn == variableNames.end()) //new line of data
		{
			vn = variableNames.begin();
			++subjectNo; //for the next subject

			//read in IDs, to advance to data 
			for(unsigned int idNo = 1; idNo <= numberOfIDs; ++idNo)
			{			
				getStringInputData(readData, notNumber);
				if(!readData.eof()) dataIDs.push_back(notNumber);
				else break;			
			};
		};
	
		//read in data for next time		
		if(!readData.eof())
		{
			getNumberInputData(readData, aNum, missingData, notNumber);
		};

	};

	readData.close();

	taskCentral->checkIDs(filename, dataIDs);
};

//! Loads data for discrete data.
void InputDataTask::loadDiscreteData(const bool & isSNPData)
{
	//first read in header and create Nodes
	list<string> variableNames, secondLine;
	DiscreteData * aDiscreteData;
	list<DiscreteData *> discreteDatas;
	list<DiscreteData *>::iterator dd;
	list<string> dataIDs;

	string aString;
	bool missingData = false;

	getVariableNames(variableNames); 
	getVariableNames(secondLine, true);
	if(variableNames.size() != secondLine.size())
	{
		outErr("There are "); outErr(variableNames.size()); outErr(" variables in the header but only ");
		outErr(secondLine.size()); outErr(" variables in the next line for file "); outErr(filename); outErr("!");
		exitErr("");
	};

	secondLine.clear();

	ifstream readData(filename.c_str());

	if(!readData.is_open())
	{
		string mess = "Cannot read data file: " + filename + "!";
		exitErr(mess);
	};

	//set name of IDs for possible data output
	list<string> nameOfIDs;

	//read in IDs, to advance to variable names below
	for(unsigned int idNo = 1; idNo <= numberOfIDs; ++idNo)
	{
		getStringInputData(readData, aString);
		nameOfIDs.push_back(aString);
	};

	taskCentral->setNameOfIDs(nameOfIDs);

	//create discrete Nodes
	list<string>::const_iterator vn = variableNames.begin();
	for( ; vn != variableNames.end(); ++vn)
	{
		//read in header, check if header exists
		getStringInputData(readData, aString);

		if(isIncluded(*vn))
		{		
			aDiscreteData = new DiscreteData(*vn, fileno);
			if(isSNPData) aDiscreteData->isSNPData = true;
			taskCentral->addNodeData(aDiscreteData); //add to general list for access later
			discreteDatas.push_back(aDiscreteData);
			++totalNodeDatas;
		};
	};

	
	//loop thro subjects
	unsigned int subjectNo = 1; 
	dd = discreteDatas.begin();
	vn = variableNames.begin();

	//read in IDs, to advance to data
	for(unsigned int idNo = 1; idNo <= numberOfIDs; ++idNo)
	{
		getStringInputData(readData, aString);
		dataIDs.push_back(aString);
	};

	getStringInputData(readData, aString);

	while(!readData.eof())
	{
		missingData = (aString == discreteMissingValue);

		if(isIncluded(*vn))
		{
			//find corresponding discrete Data
			while((*dd)->name != *vn)
			{
				++dd;
				if(dd == discreteDatas.end()) dd = discreteDatas.begin();
			};
			(*dd)->addData(aString, missingData);
		};

	
		//update variable
		++vn;
		if(vn == variableNames.end()) //new line of data
		{
			vn = variableNames.begin();
			++subjectNo; //for the next subject
			//read in IDs, to advance to data 
			for(unsigned int idNo = 1; idNo <= numberOfIDs; ++idNo)
			{
				getStringInputData(readData, aString);
				if(!readData.eof()) dataIDs.push_back(aString);
				else break;			
			};
		};
	
		//read in data for next time		
		if(!readData.eof())
		{
			getStringInputData(readData, aString);
		};

	};

	readData.close();

	taskCentral->checkIDs(filename, dataIDs);
};

//! Loads data for discrete data but records as cts factor data.
void InputDataTask::loadDiscreteDataAsFactors(const bool & isSNPData)
{
	//first read in header and create Nodes
	list<string> variableNames, secondLine;
	CtsData * aCtsFactorData;
	list<CtsData *> ctsFactorDataHeads;
	list<CtsData *>::iterator cd;
	list<string> dataIDs;
	map<string, map<string, unsigned int> > varStringToFactorNumber; //variable name, factor name, number
	map<string, map<string, unsigned int> >::iterator vstfn;
	map<string, unsigned int> stringToFactorNumber;
	map<string, unsigned int>::iterator stfn;
	map<string, unsigned int> stringToFactorNumberInit;
	unsigned int factorNumber;

	string aString;
	bool missingData = false;

	getVariableNames(variableNames); 
	getVariableNames(secondLine, true);
	if(variableNames.size() != secondLine.size())
	{
		outErr("There are "); outErr(variableNames.size()); outErr(" variables in the header but only ");
		outErr(secondLine.size()); outErr(" variables in the next line for file "); outErr(filename); outErr("!");
		exitErr("");
	};

	secondLine.clear();

	ifstream readData(filename.c_str());

	if(!readData.is_open())
	{
		string mess = "Cannot read data file: " + filename + "!";
		exitErr(mess);
	};

	//set name of IDs for possible data output
	list<string> nameOfIDs;

	//read in IDs, to advance to variable names below
	for(unsigned int idNo = 1; idNo <= numberOfIDs; ++idNo)
	{
		getStringInputData(readData, aString);
		nameOfIDs.push_back(aString);
	};

	taskCentral->setNameOfIDs(nameOfIDs);

	//create cts Nodes as factors for discrete data
	list<string>::const_iterator vn = variableNames.begin();
	for( ; vn != variableNames.end(); ++vn)
	{
		//read in header, check if header exists
		getStringInputData(readData, aString);

		if(isIncluded(*vn))
		{		
			aCtsFactorData = new CtsData(*vn, fileno);
			aCtsFactorData->isFactorNode = true;
			aCtsFactorData->nodeName = *vn;
			if(isSNPData) aCtsFactorData->isSNPData = true;

			taskCentral->addNodeData(aCtsFactorData); //add to general list for access later
			ctsFactorDataHeads.push_back(aCtsFactorData);
			++totalNodeDatas;
		};
	};

	//loop thro subjects
	unsigned int subjectNo = 1; 
	cd = ctsFactorDataHeads.begin();
	vn = variableNames.begin();

	//read in IDs, to advance to data
	for(unsigned int idNo = 1; idNo <= numberOfIDs; ++idNo)
	{
		getStringInputData(readData, aString);
		dataIDs.push_back(aString);
	};

	getStringInputData(readData, aString);

	while(!readData.eof())
	{
		missingData = (aString == discreteMissingValue);

		if(isIncluded(*vn))
		{
			//find corresponding cts Factor Data
			while((*cd)->name != *vn)
			{
				++cd;
				if(cd == ctsFactorDataHeads.end()) cd = ctsFactorDataHeads.begin();
			};
			
			if(!missingData)
			{
				//get factor number
				vstfn = varStringToFactorNumber.find(*vn);
				if(vstfn == varStringToFactorNumber.end())
				{
					//set up map for factors to number for this variable
					stringToFactorNumberInit.clear();
					factorNumber = 0;
					stringToFactorNumberInit[aString] = factorNumber;
					varStringToFactorNumber[*vn] = stringToFactorNumberInit;
				}
				else
				{
					//find number for this factor
					stfn = vstfn->second.find(aString);
					if(stfn == vstfn->second.end())
					{
						factorNumber = vstfn->second.size();
						vstfn->second[aString] = factorNumber;
					}
					else
					{
						factorNumber = stfn->second;
					};

				};

				(*cd)->setLevelName(factorNumber, aString);
			}
			else
			{
				factorNumber = 0;
			};

			(*cd)->addData(factorNumber, missingData);
		};

		//update variable
		++vn;
		if(vn == variableNames.end()) //new line of data
		{
			vn = variableNames.begin();
			++subjectNo; //for the next subject
			//read in IDs, to advance to data 
			for(unsigned int idNo = 1; idNo <= numberOfIDs; ++idNo)
			{
				getStringInputData(readData, aString);
				if(!readData.eof()) dataIDs.push_back(aString);
				else break;			
			};
		};
	
		//read in data for next time		
		if(!readData.eof())
		{
			getStringInputData(readData, aString);
		};

	}; //end of data loop

	readData.close();

	list<string> factorNames;
	map<string, map<unsigned int, CtsData *> > headFactorToOtherFactors; 
	map<string, map<unsigned int, CtsData *> >::iterator ihftof;
	map<unsigned int, CtsData *> aListOfFactors;
	map<string, list<string> > headFactorToFactorNames;
	list<string> aListOfFactorNames;

	string factorVariableName;
	cd = ctsFactorDataHeads.begin();
	unsigned int factorCount;

	//create new cts factor data objects
	for(vstfn = varStringToFactorNumber.begin(); vstfn != varStringToFactorNumber.end(); ++vstfn, ++cd)
	{
		
		stfn = vstfn->second.begin();
		stfn++;
		if(stfn != vstfn->second.end())
		{
			//add factor name for factor number 1, the head factor
			(*cd)->factorName = stfn->first;  
			stfn++;
		};
		
		aListOfFactors.clear();
		aListOfFactorNames.clear();

		for(;stfn != vstfn->second.end(); ++stfn)
		{
			factorVariableName = vstfn->first + " " + stfn->first;
			aCtsFactorData = new CtsData(factorVariableName, fileno);
			aCtsFactorData->isFactorChildNode = true;
			aCtsFactorData->nodeName = vstfn->first;
			aCtsFactorData->factorName = stfn->first;
			if(isSNPData) aCtsFactorData->isSNPData = true;

			taskCentral->addNodeData(aCtsFactorData); //add to general list for access later
			aListOfFactors[stfn->second] = aCtsFactorData;
			aListOfFactorNames.push_back(factorVariableName);
		};

		headFactorToOtherFactors[vstfn->first] = aListOfFactors;
		headFactorToFactorNames[vstfn->first] = aListOfFactorNames;
	};

	vstfn = varStringToFactorNumber.begin();
	list<bool>::const_iterator mis;

	for(list<CtsData *>::iterator cfdh = ctsFactorDataHeads.begin(); cfdh != ctsFactorDataHeads.end(); ++cfdh, ++vstfn)
	{
		
		ihftof = headFactorToOtherFactors.find((*cfdh)->name);
		if(ihftof != headFactorToOtherFactors.end()) aListOfFactors = ihftof->second;
		else
		{
			exitErr("Error: Problem processing factors");
		};

		//loop thro' data and add data to nodes
		mis = (*cfdh)->missingValues.begin();
		for(list<double>::iterator v = (*cfdh)->values.begin(); v != (*cfdh)->values.end(); ++v, ++mis)
		{
			factorCount = 1;
			//loop thro' factors and set data
			for(map<unsigned int, CtsData *>::iterator fac = aListOfFactors.begin(); fac != aListOfFactors.end(); ++fac)
			{
				++factorCount;
				if(*v == factorCount) fac->second->addData(1, *mis);
				else fac->second->addData(0, *mis);
			};

			if(*v > 1) *v = 0; //set head factor, so 1 remain 1s, 0 remain 0s else set to 0
		};
	};


	//keep record of which data sets belong together as a factor
	for(vn = variableNames.begin(); vn != variableNames.end(); ++vn)
	{
		//vstfn = varStringToFactorNumber.find(*vn);
		//vstfn->second.size();

		//set up CtsData objects for the other factors and change data to 0 and 1s for factor variable
		//loop thro' data
		factorNames = headFactorToFactorNames.find(*vn)->second;

		taskCentral->getAllNodeData()->addFactorDataGroup(*vn, factorNames);
	};

	taskCentral->checkIDs(filename, dataIDs);
};

//! Outputs task header for inputting a network.
void InputNetworkTask::outputTaskHeader()
{
	outputTaskName();
	out("Loading network\n");
};

//! Outputs task details for inputting a network.
void InputNetworkTask::outputTaskDetails()
{
		
	if(filename != "") {out("Network file: "); out(filename); out("\n");}
	else if(filename2 != "") {out("Network file in bnlearn format: "); out(filename2); out("\n");}
	else if(filenamePrefix != "") {out("Network files in igraph format: "); out(filenamePrefix+"-nodes.dat and "+filenamePrefix+"-edges.dat"); out("\n");}
	else out("Network set with no edges\n");

	if(!allowedEdgeTypes.empty())
	{
		out("Edge types allowed in the network:\n");
		for(set<pair<string, string> >::const_iterator ae = allowedEdgeTypes.begin(); ae != allowedEdgeTypes.end(); ++ae)
		{		
			 out("  "); out(ae->first); out(" --> "); out(ae->second); out("\n");
		};
	}
	else if(!notAllowedEdgeTypes.empty())
	{
		out("Blacklist edge types:\n");
		for(set<pair<string, string> >::const_iterator nae = notAllowedEdgeTypes.begin(); nae != notAllowedEdgeTypes.end(); ++nae)
		{
			out("  "); out(nae->first); out(" --> "); out(nae->second); out("\n");
		};
	};

	if(!noParentsNodes.empty())
	{
		out("Nodes with no parents: ");
		for(set<string>::const_iterator npn = noParentsNodes.begin(); npn != noParentsNodes.end(); ++npn)
		{
			out(*npn); out("\n");
		};
	};

	if(!noChildrenNodes.empty())
	{
		out("Nodes with no children: ");
		//add no children nodes
		for(set<string>::const_iterator ncn = noChildrenNodes.begin(); ncn != noChildrenNodes.end(); ++ncn)
		{
			out(*ncn); out("\n");
		};
	};

	out("Network type: "); out(network->getNetworkType()); out("\n");
	if(network->getNetworkType() != "deal") { out("Network score type: "); out(network->getScoreTypeName()); out("\n"); };
	if(network->getScoreFixName() != "none") { out("Network score fix: "); out(network->getScoreFixName()); out("\n"); };
	if((score == "cost" || score == "BICprob") && !costEdgeTypes.empty())
	{
		if(score == "cost") out("Cost of edge types:\n");
		else out("Probability of edge types:\n");

		for(map<pair<string, string>, double>::const_iterator cet = costEdgeTypes.begin(); cet != costEdgeTypes.end(); ++cet)
		{		
			 out("  "); out(cet->first.first); out(" --> "); out(cet->first.second); out(": "); out(cet->second); out("\n");
		};
	}
	if((score == "cost" || score == "BICprob") && (!costEdges.empty()))
	{
		if(score == "cost") out("Cost of edges:\n");
		else out("Probability of edges:\n");

		for(map<pair<string, string>, double>::const_iterator ce = costEdges.begin(); ce != costEdges.end(); ++ce)
		{		
			 out("  "); out(ce->first.first); out(" --> "); out(ce->first.second); out(": "); out(ce->second); out("\n");
		};
	}

	if(blacklistFilename != "") {out("Blacklist network file: "); out(blacklistFilename); out("\n");};
	if(whitelistFilename != "") {out("Whitelist network file: "); out(whitelistFilename); out("\n");};
	out("Total number of nodes: "); out(totalNoDataNodes + totalDisNodes + totalFactorNodes + totalCtsNodes); out(" (Discrete: ");  out(totalDisNodes); out(" | Factor: ");  out(totalFactorNodes); out(" | Continuous: "); out(totalCtsNodes);
			if(totalNoDataNodes != 0) {out(" | No data: "); out(totalNoDataNodes);} out(")"); out("\n");
	out("Total number of edges: "); out(totalEdges); out("\n");
	out("Network Structure: "); out(network->getNetworkString()); out("\n");
	if(network->hasLoop()) {out("The network contains a loop!\n");};
	if(network->getScoreTypeName() == "bayes" || network->getNetworkType() == "deal") { out("Imaginary sample size: "); out(imaginarySampleSize); out("\n"); };
	if(network->hasNodesWithNoData()) out("The network has nodes with no data\n");
	else
	{
		NetworkMissingData * networkMissingData = network->updateNetworkMissingData();
		pair<unsigned int, unsigned int> missingNotMissing = network->getNoMissingNotMissing();
		out("Total data at each node: "); out(missingNotMissing.second); out("\n");
		out("Missing data at each node: "); out(missingNotMissing.first); out("\n");
		delete networkMissingData;
	};
	
};

//! Gets the node names of an edge or one node name.
void getNodes(string & aLine, string & node1, string & node2)
{
	unsigned int length = aLine.length();
	
	//first trim spaces before and after nodes
	string delimiters = " \f\n\r\t\v";

	unsigned int last = aLine.find_last_not_of( delimiters );
	if(last > length)
	{
		node1 = "";
		node2 = "";
		return;
	};

    aLine = aLine.substr(0, last + 1 );
	aLine = aLine.substr(aLine.find_first_not_of( delimiters ));

	node1 = aLine.substr(0, aLine.find_first_of( delimiters ) );
	node2 = aLine.substr(aLine.find_first_of( delimiters ) + 1);
};

//! Sets network type.
void InputNetworkTask::setNetworkType(const string & nm)
{
	if (nm != "deal" && nm != "bnlearn")
	{
		exitErr("Network type must be set to \"deal\" or \"bnlearn\"!");
	};

	networkType = nm;
};

//! Returns score fix from string.
unsigned int getNetworkType(const string & score)
{
	unsigned int scoreType = 1;
	if(score == "loglike") scoreType = 0;
	else if(score == "BIC") scoreType = 1;
	else if(score == "AIC") scoreType = 2;
	else if(score == "bayes") {scoreType = 3; exitErr("Bayes not done!"); }
	else if(score == "cost") scoreType = 4;
	else if(score == "BICprob") scoreType = 5;
	else exitErr("Network score, " + score + ", not recognized! (Can be: BIC, AIC or loglike)");

	return scoreType;
};

//! Returns score type from string.
unsigned int getNetworkScoreFix(const string & scoreFix)
{
	unsigned int scoreFixNum = 1;
	if(scoreFix == "" || scoreFix == "none") scoreFixNum = 0;
	else if(scoreFix == "average") scoreFixNum = 1;
	else if(scoreFix == "skip") scoreFixNum = 2;	
	else exitErr("Network score fix, " + scoreFix + ", not recognized! (Can be: none, average or skip)");

	return scoreFixNum;
};

//! Loads data for discrete data.
void InputNetworkTask::loadNetworkData(const unsigned int & type)
{
	if(type == 0)
	{
		network = taskCentral->getNoEdgeNetwork();
		taskCentral->addNetwork(name, network);
		return;
	};

	bool networkCreated = false;

	//create the network
	if(network == 0)
	{
		unsigned int scoreType = getNetworkType(score);
		if(networkType == "deal") network = new NetworkDeal(taskCentral->getAllNodeData(), scoreType);
		else network = new NetworkBNLearn(taskCentral->getAllNodeData(), scoreType, getNetworkScoreFix(scoreFix));

		networkCreated = true;

	};

	if(filename != "" || whitelistFilename != "" || blacklistFilename != "") loadNetworkDataFormat1(type); 
	else if(filename2 != "") loadNetworkDataFormat2(); 
	else if(filenamePrefix != "") loadNetworkDataFormat3(); 

	if(networkCreated) taskCentral->addNetwork(name, network); //add to general list for access later
};

//Load network structure using the standard BayesNetty format.
void InputNetworkTask::loadNetworkDataFormat1(const unsigned int & type)
{

	string theFilename;
	if(type == 1) theFilename = filename;
	else if(type == 2) theFilename = whitelistFilename;
	else if(type == 3) theFilename = blacklistFilename;

	//read in network in format 1
	string aLine, node1, node2;

	//read in nodes and edges
	ifstream readNetwork(theFilename.c_str());

	if(!readNetwork.is_open())
	{
		string mess = "Cannot read network file: " + theFilename + "!";
		exitErr(mess);
	};

	//read in one line and add edge or node
	do
	{
		getline(readNetwork, aLine);
		getNodes(aLine, node1, node2);

		//add node or edge
		if(node1 == "") {}
		else if(node1 == node2)
		{
			if(type == 1)
			{				
				if(!network->nodeExistsInit(node1))	network->addNodeInit(node1);
				else
				{
					out("Warning: attempt to add node "); out(node1); out(" twice!\n");
				};
			}
			else if(type == 2) network->addWhiteNode(node1);
			else if(type == 3) network->addBlackNode(node1);			
		}
		else if(node2 != "")
		{
			if(type == 1)
			{
				if(network->nodeExistsInit(node1) && network->nodeExistsInit(node2))
				{
					if(!network->edgeExistsInit(node1, node2)) network->addEdgeInit(node1, node2);
					else
					{
						out("Warning: attempt to add edge "); out(node1); out("-->"); out(node2); out(" twice!\n");
					};
				}
				else
				{
					string mess = "Cannot add edge between nodes: " + node1 + " and " + node2 + "!";
					exitErr(mess);
				};
			}
			else if(type == 2) network->addWhiteEdge(node1, node2);
			else if(type == 3) network->addBlackEdge(node1, node2);
		};


	}while(!readNetwork.eof());

	readNetwork.close();
};

//Load network structure using the standard bnlearn format.
void InputNetworkTask::loadNetworkDataFormat2()
{
	//read in network firstly
	string networkString;
	ifstream readNetwork(filename2.c_str());

	if(!readNetwork.is_open())
	{
		string mess = "Cannot read network file: " + filename2 + "!";
		exitErr(mess);
	};

	readNetwork >> networkString;

	readNetwork.close();

	network->setNetwork(networkString);
};

//Load network structure using the igraph format.
void InputNetworkTask::loadNetworkDataFormat3()
{
	//read in nodes first
	string nodeFileName = filenamePrefix + "-nodes.dat";
	string id, name, type, fileno;
	map<string, string> idName;
	map<string, string>::const_iterator i, j;
	ifstream readNetworkNodes(nodeFileName.c_str());

	if(!readNetworkNodes.is_open())
	{
		string mess = "Cannot read network node file: " + nodeFileName + "!";
		exitErr(mess);
	};

	//read in header
	readNetworkNodes >> id >> id >> id >> id;
	
	//read in first line
	readNetworkNodes >> id >> name >> type >> fileno;

	while(!readNetworkNodes.eof())
	{
		
		i = idName.find(id);
		if(i != idName.end())
		{
			outErr("Node "); outErr(name); out(" and "); outErr(i->second); out(" have the same id "); outErr(id);  out("!\n");
			exit(1);
		} 
		else if(network->nodeExistsInit(name))
		{
			outErr("Node "); outErr(name); out(" is repeated in the network file!\n");
			exit(1);
		};
		
		idName[id] = name;
	
		network->addNodeInit(name);

		//read next node
		readNetworkNodes >> id >> name >> type >> fileno;
	};

	readNetworkNodes.close();

	//now read in edges 
	string edgeFileName = filenamePrefix + "-edges.dat";
	ifstream readNetworkEdges(edgeFileName.c_str());

	if(!readNetworkEdges.is_open())
	{
		string mess = "Cannot read network edge file: " + edgeFileName + "!";
		exitErr(mess);
	};

	string id1, id2, chisq, name1, name2;

	//read in header
	readNetworkEdges >> id >> id >> id;

	//read in first line
	readNetworkEdges >> id1 >> id2 >> chisq;

	while(!readNetworkEdges.eof())
	{
		i = idName.find(id1);
		j = idName.find(id2);
		
		if(i == idName.end())
		{
			outErr("Node with id "); outErr(id1); outErr(" not found in list of nodes!\n");
			exit(1);
		} 
		else if(j == idName.end())
		{
			outErr("Node with id "); outErr(id2); outErr(" not found in list of nodes!\n");
			exit(1);
		};

		//add the edge		
		network->addEdgeInit(i->second, j->second);

		//read in next line
		readNetworkEdges >> id1 >> id2 >> chisq;
	};

	readNetworkEdges.close();
};

//! Inputs network.
void InputNetworkTask::doTask()
{
	taskCentral->setDefaultScoreType(getNetworkType(score)); //ensure any default networks, empty etc, use netscore if set

	//setup nodes and edges
	if(emptyNetwork && filename != "") exitErr("You may not specify a network to have no edges and set a network file!");
	if((filename != "" && (filename2 != "" || filenamePrefix != "") || (filename2 != "" && filenamePrefix != ""))) exitErr("You may not specify two input network files for the same network!"); 
	if(!emptyNetwork && filename == "" && filename2 == "" && filenamePrefix == "") exitErr("You must specify an input network file or declare a no-edge network!"); 

	if(emptyNetwork) loadNetworkData(0);
	if(filename != "" || filename2 != "" || filenamePrefix != "") loadNetworkData(1);
	if(whitelistFilename != "") loadNetworkData(2);
	if(blacklistFilename != "") loadNetworkData(3);

	set<pair<unsigned int, unsigned int> > allowedEdgeTypesFileNos; // fileNo1, fileNo2
	set<pair<unsigned int, unsigned int> > notAllowedEdgeTypesFileNos;
		
	//set up corresponding file numbers
	for(set<pair<string, string> >::const_iterator ae = allowedEdgeTypes.begin(); ae != allowedEdgeTypes.end(); ++ae)
	{		
		allowedEdgeTypesFileNos.insert(make_pair(taskCentral->getFileNoFromDataName(ae->first), taskCentral->getFileNoFromDataName(ae->second)));
	};

	for(set<pair<string, string> >::const_iterator nae = notAllowedEdgeTypes.begin(); nae != notAllowedEdgeTypes.end(); ++nae)
	{
		notAllowedEdgeTypesFileNos.insert(make_pair(taskCentral->getFileNoFromDataName(nae->first), taskCentral->getFileNoFromDataName(nae->second)));
	};

	network->setUpAllowedEdgeTypes(allowedEdgeTypesFileNos, notAllowedEdgeTypesFileNos);

	//set up cost of edges in network
	for(map<pair<string, string>, double>::const_iterator cet = costEdgeTypes.begin(); cet != costEdgeTypes.end(); ++cet)
	{		
			network->setCostEdgeType(taskCentral->getFileNoFromDataName(cet->first.first), taskCentral->getFileNoFromDataName(cet->first.second), cet->second);
	};
	
	for(map<pair<string, string>, double>::const_iterator ce = costEdges.begin(); ce != costEdges.end(); ++ce)
	{		
			network->setCostEdge(ce->first.first, ce->first.second, ce->second);
	};

	//add no parents nodes
	for(set<string>::const_iterator npn = noParentsNodes.begin(); npn != noParentsNodes.end(); ++npn)
	{
		network->addNoParentsNode(*npn);
	};

	//add no children nodes
	for(set<string>::const_iterator ncn = noChildrenNodes.begin(); ncn != noChildrenNodes.end(); ++ncn)
	{
		network->addNoChildrenNode(*ncn);
	};

	//now setup network for use
	network->initialise();	
	network->getNumberNodesAndEdges(totalDisNodes, totalFactorNodes, totalCtsNodes, totalNoDataNodes, totalEdges);
	network->setImaginarySampleSize(imaginarySampleSize);
	
	//set if using imputed data for plotting the network later
	if(usingImputedData) network->setIsImputedData();

};
	
//! Initialises task for inputting a network.
void InputNetworkTask::initialiseTask()
{
	
	if((!costEdges.empty() || !costEdgeTypes.empty()) && score != "BICprob")
	{
		string msg = "If edge probabilities are used then the network score must be set to \"BICprob\", but you have set it to "+ score +"!\n"; 
		exitErr(msg);
	};

	if(score == "BICprob") //edge probabilities, BICprob
	{
		
		for(map<pair<string, string>, double>::const_iterator ce = costEdges.begin(); ce != costEdges.end(); ++ce)
		{

			if(ce->second > 1 || ce->second  < 0)
			{
				exitErr("Edge probabilites must be set between 0 and 1!\n");
			};

			//check if reverse exists, if not then add it
			map<pair<string, string>, double>::const_iterator ce1 = costEdges.find(make_pair(ce->first.second, ce->first.first));

			if(ce1 == costEdges.end())
			{
				costEdges[make_pair(ce->first.second, ce->first.first)] = 1 - ce->second;
			}
			else
			{
				if(ce->second + ce1->second > 1)
				{
					string msg = "Edge probabilities between " + ce->first.first + " and " + ce->first.second + " sum to over 1!\n"; 
					exitErr(msg);
				}
			}
		};

		for(map<pair<string, string>, double>::const_iterator cet = costEdgeTypes.begin(); cet != costEdgeTypes.end(); ++cet)
		{

			if(cet->second > 1 || cet->second  < 0)
			{
				exitErr("Edge probabilites must be set between 0 and 1!\n");
			};

			//check if reverse exists, if not then add it
			map<pair<string, string>, double>::const_iterator cet1 = costEdgeTypes.find(make_pair(cet->first.second, cet->first.first));

			if(cet1 == costEdgeTypes.end())
			{
				costEdgeTypes[make_pair(cet->first.second, cet->first.first)] = 1 - cet->second;
			}
			else
			{
				if(cet->second + cet1->second > 1)
				{
					string msg = "Edge type probabilities between " + cet->first.first + " and " + cet->first.second + " sum to over 1!\n"; 
					exitErr(msg);
				}
			}
		};


	};
};

//! Outputs task header for calculating posteriors.
void CalculatePosteriorTask::outputTaskHeader()
{
	outputTaskName();
	out("Calculating posterior\n");	
};

//! Outputs task details for calculating posteriors.
void CalculatePosteriorTask::outputTaskDetails()
{
	
	out("Network: "); out(networkName); out("\n");
	out("Network Structure: "); out(network->getNetworkString()); out("\n");
	if(hasLoop) {out("The network contains a loop!\n");};
};

//! Calculates posteriors.
void CalculatePosteriorTask::doTask()
{
	//choose network to perform analysis
	if(networkName != "") network =	taskCentral->getNetwork(networkName);
	else network = taskCentral->getLatestNetwork(networkName);

	//check all of the Nodes have the same amount of data
	taskCentral->checkData();
	taskCentral->addNetworkMissingData(network->updateNetworkMissingData());

	hasLoop = network->hasLoop();

	if(!hasLoop)
	{		
		network->calculatePriorsAndPosteriors();
	};
};
	
//! Initialises task for calculating a Markov blanket.
void CalculateMarkovBlanketTask::initialiseTask()
{
		
		
};

//! Outputs task header for calculating a Markov blanket.
void CalculateMarkovBlanketTask::outputTaskHeader()
{
	outputTaskName();
	out("Calculating Markov blanket\n");
};

//! Outputs task details for calculating a Markov blanket.
void CalculateMarkovBlanketTask::outputTaskDetails()
{
	out("Network: "); out(networkName); out("\n");
	out("Node: "); out(nodeName); out("\n");
	out("Network structure: "); out(network->getNetworkString()); out("\n");
	out("Markov blanket network structure: "); out(blanketNetwork->getNetworkString()); out("\n");	
};
	
//! Calculates a Markov blanket.
void CalculateMarkovBlanketTask::doTask()
{
	if(nodeName == "")
	{
		outErr("A node must be chosen with option -markov-blanket-node-name!\n");
		exit(1);
	};

	//choose network to perform analysis
	if(networkName != "") network =	taskCentral->getNetwork(networkName);
	else network = taskCentral->getLatestNetwork(networkName);

	blanketNetwork = network->getMarkovBlanketSubNetwork(nodeName); //includes setting up master prior and posterior

	taskCentral->addNetwork(name, blanketNetwork);
};

//! Initialises task for calculating a network score.
void CalculateNetworkScoreTask::initialiseTask()
{
		
		
};

//! Outputs task header for calculating a network score.
void CalculateNetworkScoreTask::outputTaskHeader()
{
	outputTaskName();
	out("Calculating network score\n");
};

//! Outputs task details for calculating a network score.
void CalculateNetworkScoreTask::outputTaskDetails()
{
	out("Network: "); out(networkName); out("\n");
	if(allScoresFilename == "")
	{
		out("Network structure: "); out(network->getNetworkString()); out("\n");
		out("Network score type: "); out(network->getScoreTypeName()); out("\n");
		if(filename != "") {out("Filename: "); out(filename); out("\n");};
		if(network->getScoreFixName() != "none") { out("Network score fix: "); out(network->getScoreFixName()); out("\n"); };
		if(hasLoop) {out("The network contains a loop!\n");}
		else {out("Network score = "); out(networkScore); out("\n"); };
	}
	else
	{
		out("All network scores filename: "); out(allScoresFilename); out("\n");
		out("Number of networks evaluated: "); out(noNetworksEval); out("\n");
		out("Network score type: "); out(network->getScoreTypeName()); out("\n");
		out("Best network structure: "); out(bestNetworkStructure); out("\n");
		out("Best network score = "); out(networkScore); out("\n");
	};
};
	
//! Calculates a network score.
void CalculateNetworkScoreTask::doTask()
{
	//choose network to perform analysis
	if(networkName != "") network =	taskCentral->getNetwork(networkName);
	else network = taskCentral->getLatestNetwork(networkName);

	//check all of the Nodes have the same amount of data
	taskCentral->checkData();
	taskCentral->addNetworkMissingData(network->updateNetworkMissingData());

	if(allScoresFilename != "")
	{
		noNetworksEval = calculateAllScores();
	}
	else
	{
		hasLoop = network->hasLoop();

		if(!hasLoop)
		{
			network->calculatePriorsAndPosteriorsBNL();

			networkScore = network->calcScore(); //includes setting up master prior and posterior

			if(filename != "")
			{
				ofstream scoresFile(filename.c_str());

				scoresFile << networkScore << "\n";	

				scoresFile.close();
			};
		};
	};


};

//! Calculates a network score for all possible networks.
unsigned int CalculateNetworkScoreTask::calculateAllScores()
{
	return network->calculateAllScores(allScoresFilename, networkScore, bestNetworkStructure);
};

//! Initialises task for searching for the best network.
void SearchNetworkModelsTask::initialiseTask()
{
		
};

//! Outputs task header for searching for the best network.
void SearchNetworkModelsTask::outputTaskHeader()
{
	outputTaskName();
	out("Searching network models\n");	
};

//! Outputs task details for searching for the best network.
void SearchNetworkModelsTask::outputTaskDetails()
{
	out("Network: "); out(networkName); out("\n");
	out("Search: "); out(searchType); out("\n");
	out("Random restarts: "); out(randomRestarts); out("\n");
	if(randomRestarts > 0)
	{
		out(" Random restart minimum edge/node ratio:"); out(minEdgeNodeRatio); out("\n");
		out(" Random restart maximum edge/node ratio:"); out(maxEdgeNodeRatio); out("\n");
		out(" Random restart minimum number of edges:"); out(minEdges); out("\n");
		out(" Random restart maximum number of edges:"); out(maxEdges); out("\n");
	};
	out("Random jitter restarts: "); out(jitterRestarts); out("\n");
	if(jitterRestarts > 0)
	{
		out(" Random jitter restart minimum edge/node ratio:"); out(minEdgeNodeRatioJitter); out("\n");
		out(" Random jitter restart maximum edge/node ratio:"); out(maxEdgeNodeRatioJitter); out("\n");
		out(" Random jitter restart minimum number of edges:"); out(minEdgesJitter); out("\n");
		out(" Random jitter restart maximum number of edges:"); out(maxEdgesJitter); out("\n");
	};
	out("Network Structure: "); out(bestNetwork->getNetworkString()); out("\n");
	out("Network score type: "); out(bestNetwork->getScoreTypeName()); out("\n");
	out("Network score = "); out(bestNetworkScore); out("\n");
	if(filename!="") {out("Network search output to file: "); out(filename); out("\n");};
};
	
//! Searches for the best network.
void SearchNetworkModelsTask::doTask()
{
	//choose network for search
	if(networkName != "") bestNetwork =	taskCentral->getNetwork(networkName);
	else bestNetwork = taskCentral->getLatestNetwork(networkName);

	searchNetworks = new GreedySearch();
	
	//check all of the Nodes have the same amount of data
	taskCentral->checkData(); //do not update missing data here
	taskCentral->addNetworkMissingData(bestNetwork->updateNetworkMissingData());

	double networkScore; //best = greatest loglike

	//output search if req
	bool outputSearch = (filename != "");
	if(outputSearch) searchNetworks->setSearchOutput(filename);

	if(bestNetwork->getNetworkType() == "deal")
	{
		restartNetwork = new NetworkDeal(taskCentral->getAllNodeData(), bestNetwork->getScoreType());
	}
	else {
		restartNetwork = new NetworkBNLearn(taskCentral->getAllNodeData(), bestNetwork->getScoreType(), bestNetwork->getScoreFix());
	};

	taskCentral->setupRestartNetwork(restartNetwork, bestNetwork, preNodeName);

	//keep record of the best network score and network
	bestNetworkScore = searchNetworks->doSearch(bestNetwork);

	unsigned int noNodes = bestNetwork->getNoNodes();
	minEdges = (unsigned int)(minEdgeNodeRatio*((double)noNodes) + 0.5);
	maxEdges = (unsigned int)(maxEdgeNodeRatio*((double)noNodes) + 0.5);
	minEdgesJitter = (unsigned int)(minEdgeNodeRatioJitter*((double)noNodes) + 0.5);
	maxEdgesJitter = (unsigned int)(maxEdgeNodeRatioJitter*((double)noNodes) + 0.5);


	//random restarts
	unsigned int randomRestartsToDo = randomRestarts;
	unsigned int jitterRestartsToDo = jitterRestarts;
	unsigned int randomRestartsToDoPrev = 0;
	unsigned int jitterRestartsToDoPrev = 0;

	bool foundValidNetwork = false;

	while(!(randomRestartsToDo == 0 && jitterRestartsToDo == 0))
	{

		for(unsigned int restartNo = 1; restartNo <= 20; ++restartNo)
		{

			//do jitters
			if(jitterRestartsToDo > 0)
			{
			
				taskCentral->updateRestartNetwork(restartNetwork, bestNetwork, minEdgesJitter, maxEdgesJitter, true);
		
				if(outputSearch && jitterRestartsToDoPrev != jitterRestartsToDo){
						searchNetworks->outputSearchText("0 jitter_restart_");
						searchNetworks->outputSearchNum(randomRestarts-randomRestartsToDo);
						searchNetworks->outputSearchText("_");
						searchNetworks->outputSearchNum(jitterRestarts-jitterRestartsToDo+1);
						searchNetworks->outputSearchText("\n");
				};

				networkScore = searchNetworks->doSearch(restartNetwork);

				if(networkScore*0 == 0)
				{
					jitterRestartsToDo--;
					break;
				};
				
				jitterRestartsToDoPrev = randomRestartsToDo;

			} //end of jitter
			else
			{
		
				taskCentral->updateRestartNetwork(restartNetwork, bestNetwork, minEdges, maxEdges, false);

				if(outputSearch && randomRestartsToDoPrev != randomRestartsToDo){searchNetworks->outputSearchText("0 random_restart_"); searchNetworks->outputSearchNum(randomRestarts-randomRestartsToDo+1); searchNetworks->outputSearchText("\n");};
			
				networkScore = searchNetworks->doSearch(restartNetwork);
		
				if(networkScore*0 == 0)
				{
					randomRestartsToDo--;
					
					if(jitterRestarts > 0)
					{
						jitterRestartsToDo = jitterRestarts;
						jitterRestartsToDoPrev = 0;
					};

					break; //if a valid search then stop, if a nan then try again up to 20 times
				};			

				randomRestartsToDoPrev = randomRestartsToDo;

			}; //end of random restart
	
		}; //end of trying to find restart network, jitters or random restarts

		//failed to find valid network after 20 tries, reduce counters anyway to avoid infinite loop
		if(networkScore*0 != 0)
		{
			if(jitterRestartsToDo > 0) jitterRestartsToDo--;
			else 
			{
				randomRestartsToDo--;
				if(jitterRestarts > 0) jitterRestartsToDo = jitterRestarts;
			};		
		}
		else if(networkScore > bestNetworkScore)
		{
			//set updates of best the same as the restart network
			taskCentral->copyNetworkEdges(bestNetwork, restartNetwork);
			bestNetworkScore = networkScore;			
		};

	};

	if(copyCoeffs)
	{
		bestNetwork->calcScore(); //reset coeffs to correct ones for use later			
	};

	if(outputSearch)
	{
		searchNetworks->outputSearchText("0 final_loglike\n");
		searchNetworks->outputSearchNum(bestNetworkScore);
		searchNetworks->outputSearchText(" ");
		searchNetworks->outputSearchText(bestNetwork->getNetworkString(0));
		searchNetworks->outputSearchText("\n");

		searchNetworks->closeFileSearch();
	};

	delete restartNetwork;
	delete searchNetworks; 
};

//! Initialises the task for averaging networks.
void AverageNetworksTask::initialiseTask()
{
		
		
};

//! Outputs task header for averaging networks.
void AverageNetworksTask::outputTaskHeader()
{
	outputTaskName();
	if(useNetworkScoreMethod) out("Calculating average network using network score method\n");
	else if(useNetworkWeightMethod)	out("Calculating average network using weighted method\n");
	else out("Calculating average network using bootstrapping\n");	

};

//!  Outputs task details for averaging networks.
void AverageNetworksTask::outputTaskDetails()
{
	out("Network: "); out(networkName); out("\n");
	if(!useNetworkScoreMethod)
	{
		out("Number of bootstrap iterations: "); out(noBootstraps); out("\n");
		out("Random restarts: "); out(randomRestarts); out("\n");
	
		out("Random jitter restarts: "); out(jitterRestarts); out("\n");
	};
	
	if(filename!="") {out("Average network output to file: "); out(filename); out("\n");};	
	if(likelihoodFilename!="") {out("Likelihoods of separate bootstrap networks output to file: "); out(likelihoodFilename); out("\n");};	
	if(igraphPrefix!="") {out("R code to plot average network: "); out(igraphPrefix+".R"); out("\n");};	

	if(arcThresholdSet) {out("Set edge threshold: "); out(arcThreshold); out("\n");}
	else {out("Estimated edge threshold: "); out(arcThreshold); out("\n");};
	if(thresholdFilename != "") { out("Estimated edge threshold output to file: "); out(thresholdFilename); out("\n"); };
	out("Network structure (after above threshold): "); out(network->getNetworkString()); out("\n");
	out("Network score type: "); out(network->getScoreTypeName()); out("\n");
	if(network->getScoreFixName() != "none") { out("Network score fix: "); out(network->getScoreFixName()); out("\n"); };
	if(network->hasLoop()) {out("The network contains a loop!\n");}
	if(!network->checkNetworkIsValid()) {out("The network is not valid!\n");}
	else {out("Network score = "); out(network->calcScore()); out("\n"); };
		
};

Network * BootstrapNetworkTaskHelp::createBootstrapNetwork(TaskCentral * taskCentral, Network * network, const string & name)
{
	return taskCentral->createSimNetworkFromNetwork(network, name, true);
};

//! Sets up bootstrap network.
void BootstrapNetworkTaskHelp::setupBootstrap(TaskCentral * taskCentral, Network * bootstrapNetwork,  Network * network, const bool & useRandomData, const unsigned int & percent)
{
	bootstrapNetwork->setIsBootstrapData();
	bootstrapNetwork->setImaginarySampleSize(network->getImaginarySampleSize());

	NetworkMissingData * networkMissing = network->getNetworkMissingData();

	unsigned int count = 1;
	unsigned int dataCount = 1;

	//do not want to choose missing data when bootstrapping
	for(list<bool>::const_iterator m = networkMissing->missing.begin(); m != networkMissing->missing.end(); ++m, ++dataCount)
	{
		if(!*m || useRandomData)
		{
			dataForBootstrap[count] = dataCount; //map count to non-missing data that may be chosen for bootstrap
			++count;
		};		

		noTimesDataForBootstrap[dataCount] = 0;
	};

	
	pair<unsigned int, unsigned int> missingAndNotMissing = network->getNoMissingNotMissing();
	
	totalMissingData = missingAndNotMissing.first;
	totalNonMissingData = missingAndNotMissing.second;

	unsigned int bootstrapTotalMissingData = totalMissingData;
	unsigned int bootstrapTotalNonMissingData = totalNonMissingData;

	if(percent != 0)
	{
		if(useRandomData)
		{ 
			bootstrapTotalNonMissingData = (unsigned int)((double)(totalNonMissingData + totalMissingData)*((double)(percent)/100.0) + 0.5); 
		}
		else
		{
			bootstrapTotalNonMissingData = (unsigned int)((double)(totalNonMissingData)*((double)(percent)/100.0) + 0.5);	
		};

		bootstrapTotalMissingData = (totalMissingData + totalNonMissingData) - bootstrapTotalNonMissingData;
	};

	if(famFilename != "")
	{
		bootstrapTotalNonMissingData = noFamiliesNotMiss;
		bootstrapTotalMissingData = (totalMissingData + totalNonMissingData) - noFamiliesNotMiss;
	};

	//get lists of nodes
	map<unsigned int, DiscreteNode *> origNetworkDisNodes = network->getDiscreteNodes();
	map<unsigned int, DiscreteNode *> bootNetworkDisNodes = bootstrapNetwork->getDiscreteNodes();
	map<unsigned int, CtsNode *> origNetworkCtsNodes = network->getCtsNodes();
	map<unsigned int, CtsNode *> bootNetworkCtsNodes = bootstrapNetwork->getCtsNodes();

	//get lists of data
	for(map<unsigned int, DiscreteNode *>::const_iterator ondn = origNetworkDisNodes.begin(); ondn != origNetworkDisNodes.end(); ++ondn)
	{
		origDiscreteData.push_back(ondn->second->getDiscreteData());
	};

	for(map<unsigned int, DiscreteNode *>::const_iterator bndn = bootNetworkDisNodes.begin(); bndn != bootNetworkDisNodes.end(); ++bndn)
	{
		bootDiscreteData.push_back(bndn->second->getDiscreteData());
	};

	for(map<unsigned int, CtsNode *>::const_iterator oncn = origNetworkCtsNodes.begin(); oncn != origNetworkCtsNodes.end(); ++oncn)
	{
		origCtsData.push_back(oncn->second->getCtsData());
	};
	
	for(map<unsigned int, CtsNode *>::const_iterator bncn = bootNetworkCtsNodes.begin(); bncn != bootNetworkCtsNodes.end(); ++bncn)
	{
		bootCtsData.push_back(bncn->second->getCtsData());
	};

	//set up initial missing data at end of data
	for(list<DiscreteData *>::iterator bd = bootDiscreteData.begin(); bd != bootDiscreteData.end(); ++bd)
	{	
		for(unsigned int i = 1; i <= bootstrapTotalNonMissingData; ++i)
		{
			(*bd)->values.push_back(0);
			(*bd)->missingValues.push_back(false);
		};	

		for(unsigned int j = 1; j <= bootstrapTotalMissingData; ++j)
		{
			(*bd)->values.push_back(0);
			(*bd)->missingValues.push_back(true);
		};
	};

	for(list<CtsData *>::iterator cd = bootCtsData.begin(); cd != bootCtsData.end(); ++cd)
	{	
		for(unsigned int i = 1; i <= bootstrapTotalNonMissingData; ++i)
		{
			(*cd)->values.push_back(0);
			(*cd)->missingValues.push_back(false);
		};	

		for(unsigned int j = 1; j <= bootstrapTotalMissingData; ++j)
		{
			(*cd)->values.push_back(0);
			(*cd)->missingValues.push_back(true);
		};
	};


	taskCentral->checkData();
	taskCentral->addNetworkMissingData(bootstrapNetwork->updateNetworkMissingData());
};

//! Set up st. dev.s of MEs.
void BootstrapNetworkTaskHelp::setupMEStDevs(const double & measurementError)
{
	list<CtsData *>::const_iterator ocd;
	map<string, double> variableStDevs;
	map<string, double>::const_iterator vsd;

	string variableName;
	string preVariableName;
	double stDev;

	if(measurementErrorFile != "")
	{
		ifstream readMEFile;
		readMEFile.open(measurementErrorFile.c_str());
		if(!readMEFile.is_open())
		{
			string mess = "Cannot read measurement error file: " + measurementErrorFile + "!";
			exitErr(mess);
		};

		//loop thro' st devs
		do{
			readMEFile >> variableName >> stDev;

			if(preVariableName != variableName)
			{
				variableStDevs[variableName] = stDev;
			};

			preVariableName = variableName;
		} while(!readMEFile.eof());

		readMEFile.close();

		for(ocd = origCtsData.begin(); ocd != origCtsData.end(); ++ocd)
		{
			vsd = variableStDevs.find((*ocd)->name);
			if(vsd != variableStDevs.end())
			{
				stDev = vsd->second;				
			}
			else
			{
				stDev = measurementError;
			};

			if(stDev != 0) useMeasurementError = true;
			measurementErrors.push_back(stDev);
		};
	}
	else
	{
		for(ocd = origCtsData.begin(); ocd != origCtsData.end(); ++ocd)
		{
			measurementErrors.push_back(measurementError);
		};
		if(measurementError != 0) useMeasurementError = true;
	};

};

//! Set up st. dev.s of MEs.
void BootstrapNetworkTaskHelp::setupMEStDevMultiples(const double & measurementErrorMultiple)
{
	list<CtsData *>::const_iterator ocd;
	map<string, double> variableMEStDevMultiples;
	map<string, double>::const_iterator vsd;

	string variableName;
	string preVariableName;
	double stDevMultiple;
	double stDev;

	if(measurementErrorFile != "")
	{
		ifstream readMEFile;
		readMEFile.open(measurementErrorFile.c_str());
		if(!readMEFile.is_open())
		{
			string mess = "Cannot read measurement error multiple file: " + measurementErrorFile + "!";
			exitErr(mess);
		};

		//loop thro' ME st devs multiples
		do{
			readMEFile >> variableName >> stDevMultiple;
	
			if(preVariableName != variableName)
			{
				variableMEStDevMultiples[variableName] = stDevMultiple;
			};

			preVariableName = variableName;
		} while(!readMEFile.eof());

		readMEFile.close();

		for(ocd = origCtsData.begin(); ocd != origCtsData.end(); ++ocd)
		{
			vsd = variableMEStDevMultiples.find((*ocd)->name);
			if(vsd != variableMEStDevMultiples.end())
			{
				stDevMultiple = vsd->second;
			}
			else
			{
				stDevMultiple = measurementErrorMultiple;
			};

			if(stDevMultiple != 0) useMeasurementError = true;

			stDev = stDevMultiple*(*ocd)->getStDev();
			measurementErrors.push_back(stDev);
		};
	}
	else
	{
		for(ocd = origCtsData.begin(); ocd != origCtsData.end(); ++ocd)
		{
			stDev = measurementErrorMultiple * (*ocd)->getStDev();
			measurementErrors.push_back(stDev);
		};
		if(measurementErrorMultiple != 0) useMeasurementError = true;
	};

};

//! Sets up measurement errors for each variable
void BootstrapNetworkTaskHelp::setMeasurementErrors(const map<unsigned int, double> & mes, const double & me, const map<unsigned int, double> & mesm, const double & mem, Network * network)
{
	map<unsigned int, CtsNode *> origNetworkCtsNodes = network->getCtsNodes();
	double measureErr;
	double measureErrMultiple;
	map<unsigned int, double>::const_iterator varME;
	map<unsigned int, double>::const_iterator varMEMulti;

	for(map<unsigned int, CtsNode *>::const_iterator oncn = origNetworkCtsNodes.begin(); oncn != origNetworkCtsNodes.end(); ++oncn)
	{
		varMEMulti = mesm.find(oncn->first);
		varME = mes.find(oncn->first);

		if(varMEMulti != mesm.end())
		{
			measureErrMultiple = varMEMulti->second;
			measureErr = measureErrMultiple * (oncn->second->getCtsData()->getStDev());
		}
		else if(varME != mes.end())
		{
			measureErr = varME->second;
		}
		else if(mem != 0)
		{
			measureErr = mem * (oncn->second->getCtsData()->getStDev());
		}
		else
		{
			measureErr = me;
		};

		if(measureErr != 0) useMeasurementError = true;

		measurementErrors.push_back(measureErr);
	};

};

//! Sets up family IDs for sampling from
void BootstrapNetworkTaskHelp::setupFamilyIDs(const string & famFile, TaskCentral * taskCentral, Network * network)
{
	famFilename =  famFile;

	ifstream readFamilyFile;
	readFamilyFile.open(famFilename.c_str());
	if(!readFamilyFile.is_open())
	{
		string mess = "Cannot read family file: " + famFilename + "!";
		exitErr(mess);
	};

	string famID, indivID, FatherId, MotherID, sexID, famIndivID;
	string prevFamIndivID = "";
	double phenoType;
	unsigned int noCases = 0;
	unsigned int subjectNo = 0;
	list<string> dataIDs;
	set<string> uniqueFamIDs;

	NetworkMissingData * networkMissing = network->getNetworkMissingData();

	unsigned int count = 1;
	unsigned int dataCount = 1;

	//do not want to choose missing data when bootstrapping
	list<bool>::const_iterator m = networkMissing->missing.begin();

	//loop thro' subjects 
	do{		
		readFamilyFile >> famID >> indivID >> FatherId >> MotherID >> sexID >> phenoType;
		famIndivID = famID + " " + indivID;
	
		//do not duplicate the last row
		if(famIndivID != prevFamIndivID) 
		{
			dataIDs.push_back(famID);
			dataIDs.push_back(indivID);		
			subjectNo++;			
		};

		familyIDs[subjectNo] = famID;

		if(m != networkMissing->missing.end() && !*m) uniqueFamIDs.insert(famID);

		prevFamIndivID = famIndivID;
		++count;
		++m;
	}while(!readFamilyFile.eof());

	readFamilyFile.close();

	taskCentral->checkIDs(famFilename, dataIDs);

	noFamiliesNotMiss = uniqueFamIDs.size();

};

//! Draws a random value from random dist, equivalent to R command: rnorm(1, 0, measurmentError)
double getRandomMeasurmentError(const double & measurmentError)
{
	int which = 2; //set for drawing random number
	double p = (double)rand() / (double)RAND_MAX;
	double q = 1 - p;
	int status = 0;
	double bound = 0;
	double value;
	double mean = 0;
	double sd = measurmentError;

	cdfnor(&which, &p, &q, &value, &mean, &sd, &status, &bound);

	return value;
};

//! Updates bootstrap data for network.
void BootstrapNetworkTaskHelp::updateBootstrapData(Network * bootstrapNetwork, Network * network)
{
	//family stuff
	bool useFamily = famFilename != "";
	set<unsigned int> usedIndivIDs; //indiv ID
	set<string> usedFamIDs; //fam ID
	map<unsigned int, string>::const_iterator fi;

	unsigned int noIndivstoPick;
	if(useFamily) noIndivstoPick = noFamiliesNotMiss;
	else noIndivstoPick = totalNonMissingData;

	//reset boostrap counts
	unsigned int dataNo = 1;
	map<unsigned int, unsigned int>::iterator nb;
	for(nb = noTimesDataForBootstrap.begin(); nb != noTimesDataForBootstrap.end(); ++nb, ++dataNo)
	{
		if(useBootstraps || dataNo > noIndivstoPick)
		{
			nb->second = 0;
		}
		else
		{
			nb->second = 1;
		};
	};

	unsigned int randomDataNo;
	map<unsigned int, unsigned int>::const_iterator dfb; // dataForBootstrap;

	bool dataOK = true;
	
	//choose random bootstraps
	if(useBootstraps)
	{
		for(dataNo = 1; dataNo <= noIndivstoPick; )
		{
			randomDataNo = rand() % totalNonMissingData + 1;

			//ensure diff indivs from same family are not chosen
			if(useFamily)
			{
				dataOK = false;

				if(usedIndivIDs.find(randomDataNo) != usedIndivIDs.end()) dataOK = true; //can reuse for bootstrapping
				else
				{
					fi = familyIDs.find(randomDataNo);
					if(fi == familyIDs.end())
					{
						string mess = "Problem finding family ID for individual: " + toString(randomDataNo) + "!";
						exitErr(mess);
					};

					if(usedFamIDs.find(fi->second) == usedFamIDs.end()) //do not pick if fam already used
					{
						dataOK = true;
						usedIndivIDs.insert(randomDataNo);
						usedFamIDs.insert(fi->second);
					};

				};

			};

			if(dataOK)
			{
				dfb = dataForBootstrap.find(randomDataNo);
				nb = noTimesDataForBootstrap.find(dfb->second);
				nb->second++;
				++dataNo;
			};

		};
	}
	else
	{

	};

	//update discrete node data
	list<DiscreteData *>::const_iterator od = origDiscreteData.begin();
	list<DiscreteData *>::iterator bd = bootDiscreteData.begin();
	list<unsigned int>::iterator bdv;
	
	for( ; od != origDiscreteData.end(); ++od, ++bd)
	{
		bdv = (*bd)->values.begin();
		nb = noTimesDataForBootstrap.begin();

		//loop thro' data
		for(list<unsigned int>::const_iterator odv = (*od)->values.begin(); odv != (*od)->values.end(); ++odv, ++nb)
		{	
			//change bootstrap data (missing data pts will never be chosen for bootstrapping)
			for(unsigned int i = 1; i <= nb->second; ++i, ++bdv)
			{
				*bdv = *odv;
			};
		};

	};

	//update cts node data
	list<CtsData *>::const_iterator oc = origCtsData.begin();
	list<CtsData *>::iterator bc = bootCtsData.begin();
	list<double>::iterator bcv;
	unsigned int test = 1;
	list<double>::const_iterator me = measurementErrors.begin();

	for( ; oc != origCtsData.end(); ++oc, ++bc)
	{
		bcv = (*bc)->values.begin();		
		nb = noTimesDataForBootstrap.begin();

		//loop thro' data
		for(list<double>::const_iterator ocv = (*oc)->values.begin(); ocv != (*oc)->values.end(); ++ocv, ++nb)
		{	
			//change bootstrap data
			for(unsigned int i = 1; i <= nb->second; ++i, ++bcv)
			{
				*bcv = *ocv;
				if(useMeasurementError)
				{
					if(!useEffectSize) *bcv += getRandomMeasurmentError(*me); //noise
					else *bcv *= (1 + getRandomMeasurmentError(*me)); //effect size noise
				}
			};	
		};
		
		if(useMeasurementError) ++me;
	};

	//clear caches for network evaluations
	bootstrapNetwork->clearCache();
	network->clearCache();
};

//! Updates bootstrap data for network by taking a percentage subset of data without replacement - technically not a bootstrap
void BootstrapNetworkTaskHelp::updateBootstrapDataSubset(Network * bootstrapNetwork, Network * network, const bool & useRandomData, const unsigned int & percent)
{
	//reset boostrap counts
	map<unsigned int, unsigned int>::iterator nb;
	for(nb = noTimesDataForBootstrap.begin(); nb != noTimesDataForBootstrap.end(); ++nb)
	{
		nb->second = 0;
	};

	unsigned int randomDataNo;
	map<unsigned int, unsigned int>::const_iterator dfb; // dataForBootstrap;

	//bool dataOK = true;
	vector<unsigned int> availableNums;
	unsigned int pickedNum;

	unsigned int noIndivstoPick;
	if(useRandomData)
	{
		unsigned int totalData = totalNonMissingData + totalMissingData;
		noIndivstoPick = (unsigned int)((double)(totalData)*((double)(percent)/100.0) + 0.5);
		for(unsigned int i = 1; i <= totalData; ++i) availableNums.push_back(i);
	} 
	else
	{
		noIndivstoPick = (unsigned int)((double)(totalNonMissingData)*((double)(percent)/100.0) + 0.5);
		for(unsigned int i = 1; i <= totalNonMissingData; ++i) availableNums.push_back(i);
	};

	//choose random indivs
	unsigned int dataNo = 1;
	while(dataNo <= noIndivstoPick)
	{
		pickedNum = rand() % availableNums.size();
	
		randomDataNo = availableNums[pickedNum];
		availableNums.erase(availableNums.begin() + pickedNum);
			
		dfb = dataForBootstrap.find(randomDataNo);
		nb = noTimesDataForBootstrap.find(dfb->second);
			
		if(nb->second == 0) //taking a subsample
		{
			nb->second++;
			++dataNo;			
		};
			
	};

	//update discrete node data
	list<DiscreteData *>::const_iterator od = origDiscreteData.begin();
	list<DiscreteData *>::iterator bd = bootDiscreteData.begin();
	list<unsigned int>::iterator bdv;	
	list<bool>::const_iterator modv;

	for( ; od != origDiscreteData.end(); ++od, ++bd)
	{
		bdv = (*bd)->values.begin();
		nb = noTimesDataForBootstrap.begin();

		modv = (*od)->missingValues.begin();
		//loop thro' data
		for(list<unsigned int>::const_iterator odv = (*od)->values.begin(); odv != (*od)->values.end(); ++odv, ++nb, ++modv)
		{	
			//change bootstrap data (missing data pts will never be chosen for subsample - unless doing random replace)
			for(unsigned int i = 1; i <= nb->second; ++i, ++bdv)
			{
				if(!*modv) *bdv = *odv;
				else
				{
					*bdv = (*od)->getRandomLevel(); //draws from complete data
				};
			};
		};

	};

	//update cts node data
	list<CtsData *>::const_iterator oc = origCtsData.begin();
	list<CtsData *>::iterator bc = bootCtsData.begin();
	list<double>::iterator bcv;
	list<bool>::const_iterator mocv;
	unsigned int test = 1;

	for( ; oc != origCtsData.end(); ++oc, ++bc)
	{
		bcv = (*bc)->values.begin();		
		nb = noTimesDataForBootstrap.begin();

		mocv = (*oc)->missingValues.begin();
		//loop thro' data
		for(list<double>::const_iterator ocv = (*oc)->values.begin(); ocv != (*oc)->values.end(); ++ocv, ++nb, ++mocv)
		{	

			//change bootstrap data
			for(unsigned int i = 1; i <= nb->second; ++i, ++bcv)
			{				
				if(!*mocv) *bcv = *ocv;
				else					
				{
					*bcv = (*oc)->getRandomCompleteValue(); //samples from complete data
				};
			};	
		};
	};

	//clear caches for network evaluations
	bootstrapNetwork->clearCache();
	network->clearCache();
};

//! Averages network.
void AverageNetworksTask::doTask()
{
	//use the network score method if told to do so
	if(useNetworkScoreMethod)
	{
		doNetworkScoreMethod();
		return;
	};

	//choose network for analyses
	if(networkName != "") network = taskCentral->getNetwork(networkName);
	else network = taskCentral->getLatestNetwork(networkName);
	
	bootstrapNetworkTaskHelp.setUseBootstraps(useBootstraps);
	if(measurementErrorMultipleFile != "")
	{
		bootstrapNetworkTaskHelp.setMeasurementErrorFile(measurementErrorMultipleFile);
		bootstrapNetworkTaskHelp.setupMEStDevMultiples(measurementErrorMultiple);
	}
	else if(measurementErrorFile != "")
	{
		bootstrapNetworkTaskHelp.setMeasurementErrorFile(measurementErrorFile);
		bootstrapNetworkTaskHelp.setupMEStDevs(measurementError);
	}
	else bootstrapNetworkTaskHelp.setMeasurementErrors(measurementErrors, measurementError, measurementErrorMultiples, measurementErrorMultiple, network);

	bootstrapNetworkTaskHelp.setUseEffectSize(useEffectSize);

	//check all of the Nodes have the same amount of data
	taskCentral->checkData();
	taskCentral->addNetworkMissingData(network->updateNetworkMissingData());

	if(famFilename != "") bootstrapNetworkTaskHelp.setupFamilyIDs(famFilename, taskCentral, network);

	//get copy of the initial network
	string bootNetworkName = name;
	string preNodeName = name + "-";

	bootstrapNetwork = bootstrapNetworkTaskHelp.createBootstrapNetwork(taskCentral, network, bootNetworkName);

	//set up list of data for bootstrapping
	bootstrapNetworkTaskHelp.setupBootstrap(taskCentral, bootstrapNetwork, network);


	//setup search object to do the searches for the bootstrapping
	SearchNetworkModelsTask searchNetworkModelsTask;
	searchNetworkModelsTask.setNetworkName(bootNetworkName);
	searchNetworkModelsTask.setRandomRestarts(randomRestarts);
	searchNetworkModelsTask.setJitterRestarts(jitterRestarts);
	searchNetworkModelsTask.setTaskCentral(taskCentral);
	searchNetworkModelsTask.initialiseTask();
	searchNetworkModelsTask.setPreNodeName(preNodeName);


	string PDAG;
	double bootScore = 0;


	ofstream outLikes;
	unsigned int randomSeed = taskCentral->getRandomSeed();

	if(likelihoodFilename != "")
	{
		outLikes.open(likelihoodFilename.c_str());
	};

	for(unsigned int iter = 1; iter <= noBootstraps; ++iter)
	{
		srand(randomSeed++);	//ensure the same bootstrap data for every run when the same seed is set between separate runs 
		taskCentral->setRandomSeed(randomSeed);

		bootstrapNetworkTaskHelp.updateBootstrapData(bootstrapNetwork, network);

		bootstrapNetwork->removeAllEdges(); // with data changed the current
		bootstrapNetwork->updateBlackWhiteDifferentData(network, preNodeName);

		searchNetworkModelsTask.doTask();

		if(likelihoodFilename != "") outLikes << bootstrapNetwork->getScoreFromCache() << "\n";

		//object will have changed after search to load new one with same name, the best fit network
		PDAG = bootstrapNetwork->getPDAGNetworkString();
		updateArcStrengths(bootstrapNetwork);

	};

	if(likelihoodFilename != "") outLikes.close();

	outputAverageNetworkFile();
};

//! Outputs average network to a text file.
void AverageNetworksTask::outputAverageNetworkFile()
{
	map<string, pair<double, double> >::iterator tas;
	map<string, double>::iterator tac;
	double strength;
	double direction;

	//calc strength and direction for network score method so the edges are ordered properly
	if(useNetworkScoreMethod)
	{
		for(tas = arcStrengths.begin(); tas != arcStrengths.end(); ++tas)
		{		
			strength = (tas->second.first + tas->second.second) / totalProb; 
			direction = tas->second.first / (tas->second.first + tas->second.second);
		
			//update records for setting the final network
			tas->second.first = strength;
			tas->second.second = direction;
		};

	};

	//order the arcs/edges in order of strength
	multimap<double, string> orderedArcs; //strength, arc name
	
	for(map<string, pair<double, double> >::iterator tas = arcStrengths.begin(); tas != arcStrengths.end(); ++tas)
	{
		orderedArcs.insert(make_pair(tas->second.first, tas->first));
	};

	double noIts;
	if(!useNetworkScoreMethod) noIts = (double)noBootstraps;
	else noIts = (double)noNetworksEval;

	string node1, node2;
	unsigned int spacePos;
	
	Node * nodeOb1;
	Node * nodeOb2;
	Node * nodeObTemp;
	unsigned int nodeNumber1;
	unsigned int nodeNumber2;
	
	AllNodeData * allNodeData = taskCentral->getAllNodeData();

	string outFilename = filename;
	if(useMeasurementError) outFilename = filename + "-ave.dat";
	ofstream outArcs(outFilename.c_str());

	outArcs << "from\ttype1\tto\ttype2\tstrength\tdirection\n";

	multimap<double, string>::const_reverse_iterator oa = orderedArcs.rbegin();

	double maxStrength = 0;
	if(oa != orderedArcs.rend()) maxStrength = oa->first;

	//scale to between 0 and 1 for this PDAG
	for(; oa != orderedArcs.rend(); ++oa)
	{
		tas = arcStrengths.find(oa->second);
		tac = arcCounts.find(oa->second);

		if(useNetworkWeightMethod)
		{
			strength = tas->second.first / noIts;// / maxStrength; 
			direction = tas->second.second / tac->second;

			//update records for setting the final network
			tas->second.first = strength;
			tas->second.second = direction;
		}
		else if(!useNetworkScoreMethod)
		{
			strength = tas->second.first / noIts; 
			direction = tas->second.second / tac->second;

			//update records for setting the final network
			tas->second.first = strength;
			tas->second.second = direction;
		}
		else
		{
			strength = tas->second.first; 
			direction = tas->second.second;
		};
		

		spacePos = (unsigned int)tas->first.find_first_of('\t');
		node1 = tas->first.substr(0, spacePos);
		node2 = tas->first.substr(spacePos+1);
		nodeNumber1 = allNodeData->getNodeDataNumber(node1);
		nodeNumber2 = allNodeData->getNodeDataNumber(node2);

		nodeOb1 = network->getNetworkNode(nodeNumber1);
		nodeOb2 = network->getNetworkNode(nodeNumber2);
		
		if(direction < 0.5)
		{
			nodeObTemp = nodeOb1;
			nodeOb1 = nodeOb2;
			nodeOb2 = nodeObTemp;
			direction = (1 - direction);
		};
	
		outArcs << nodeOb1->getDisplayName() << "\t";

		if(nodeOb1->getIsDiscreteNode()) outArcs << "d\t";
		else if(nodeOb1->getIsFactorNode()) outArcs << "f\t";
		else outArcs << "c\t";

		outArcs << nodeOb2->getDisplayName() << "\t";
			
		if(nodeOb2->getIsDiscreteNode()) outArcs << "d\t";
		else if(nodeOb2->getIsFactorNode()) outArcs << "f\t";
		else outArcs << "c\t";
			
		outArcs << strength << "\t" << direction << "\n";

		if(!arcThresholdSet) orderedArcStrengths.push_back(strength);
		if(useMeasurementError)
		{
			fromToProb[make_pair(nodeNumber1, nodeNumber2)] = strength * direction;
		};
	};
	
	outArcs.close();

	if(useMeasurementError)
	{
		//output measurment error summaries against the "true" network befor it is overwritten
		outputNodeRobustnessSummary();
		outputEdgeRobustnessSummary();
	};

	setFinalNetwork();

	if(igraphPrefix != "") outputRGraph();

	//output file with estimated threshold in it
	if(thresholdFilename != "")
	{		
		ofstream outThres;
		outThres.open(thresholdFilename.c_str());
		outThres << arcThreshold << "\n";
		outThres.close();
	};
};

//! Get probs of edges connected to node.
void AverageNetworksTask::getNodeAveEdgeProbs(string & nodeName, double & aveParentProb, double & aveChildProb, double & aveProb)
{
	aveParentProb = 0;
	aveChildProb = 0;
	aveProb = 0;
	double parTot = 0;
	double childTot = 0;

	for(map<pair<unsigned int, unsigned int>, double>::const_iterator ftp = fromToProb.begin(); ftp != fromToProb.end(); ++ftp)
	{
		if(network->edgeExists((*ftp).first.first, (*ftp).first.second)) //check is in true network
		{
			if(network->getNetworkNode((*ftp).first.second)->getDisplayName() == nodeName)
			{
				aveParentProb += ftp->second;
				parTot++;
			} else if(network->getNetworkNode((*ftp).first.first)->getDisplayName() == nodeName)
			{
				aveChildProb += ftp->second;
				childTot++;
			}

		};
	};

	if(parTot != 0) aveParentProb /= parTot;
	else
	{
		if(network->getNumberOfParents(nodeName) > 0) aveParentProb = 0; else aveParentProb = -1;
	};

	if(childTot != 0) aveChildProb /= childTot;
	else
	{
		if(network->getNumberOfChildren(nodeName) > 0) aveChildProb = 0; else aveChildProb = -1;	
	};

	if(childTot != 0 && parTot != 0)
	{
		aveProb = (aveChildProb + aveParentProb) / (childTot + parTot);
	}
	else if(childTot != 0)
	{
		aveProb = aveChildProb;
	}
	else if(parTot != 0)
	{
		aveProb = aveParentProb;
	}
	else
	{
		if(aveParentProb == 0 || aveChildProb == 0) aveProb = 0; else aveProb = -1;
	};
};

//! Outputs node robustness summary
void AverageNetworksTask::outputNodeRobustnessSummary()
{
	list<double> measurementErrors = bootstrapNetworkTaskHelp.getMeasurementErrors();
	list<CtsData *> origCtsData = bootstrapNetworkTaskHelp.getOrigCtsData();
	string nodeName;
	//unsigned int nodeNumber;
	double stDev;
	double stDevME;
	double ratio;
	double aveParentProb;
	double aveChildProb;
	double aveProb;

	string outFilename = filename;
	if(useMeasurementError) outFilename = filename + "-nodes.dat";
	ofstream outNodes(outFilename.c_str());

	outNodes << "node\tst_dev\tME_st_dev\tratio\tave_parent_edge_prob\tave_child_edge_prob\tave_edge_prob\n";

	list<double>::const_iterator me = measurementErrors.begin();
	for(list<CtsData *>::const_iterator ocd = origCtsData.begin(); ocd != origCtsData.end(); ++ocd, ++me)
	{
		nodeName = (*ocd)->name;		
		stDev = (*ocd)->getStDev();
		stDevME = *me;
		ratio = stDevME / stDev;
		getNodeAveEdgeProbs(nodeName, aveParentProb, aveChildProb, aveProb);
		
		outNodes << nodeName << "\t" << stDev << "\t" << stDevME << "\t" << ratio << "\t";
		if(aveParentProb != -1) outNodes << aveParentProb << "\t"; else outNodes << "NA\t";
		if(aveChildProb != -1) outNodes << aveChildProb << "\t"; else outNodes << "NA\t";
		if(aveProb != -1) outNodes << aveProb << "\t"; else outNodes << "NA\t";
		outNodes << "\n";
	};

	outNodes.close();
};

//! Outputs edge robustness summary
void AverageNetworksTask::outputEdgeRobustnessSummary()
{
	list<double> measurementErrors = bootstrapNetworkTaskHelp.getMeasurementErrors();
	list<CtsData *> origCtsData = bootstrapNetworkTaskHelp.getOrigCtsData();
	map<string, double> nodeME;
	map<string, double>::const_iterator nME;
	list<double>::const_iterator me = measurementErrors.begin();
	
	double stDev = 0;
	double stDevME = 0;
	double ratio = 0;

	set<pair<unsigned int, unsigned int> > outPuttedEdges;

	string outFilename = filename;
	if(useMeasurementError) outFilename = filename + "-edges.dat";
	ofstream outEdges(outFilename.c_str());

	for(list<CtsData *>::const_iterator ocd = origCtsData.begin(); ocd != origCtsData.end(); ++ocd, ++me)
	{
		nodeME[(*ocd)->name] = *me;		
	};

	outEdges << "from\tto\tst_dev_from\tME_st_dev_from\tratio_from\tst_dev_to\tME_st_dev_to\tratio_to\tprob\n";

	for(map<pair<unsigned int, unsigned int>, double>::const_iterator ftp = fromToProb.begin(); ftp != fromToProb.end(); ++ftp)
	{
		if(network->edgeExists((*ftp).first.first, (*ftp).first.second)) //check is in true network
		{
			outPuttedEdges.insert(make_pair((*ftp).first.first, (*ftp).first.second));
			outEdges << network->getNetworkNode((*ftp).first.first)->getDisplayName() << "\t" << network->getNetworkNode((*ftp).first.second)->getDisplayName() << "\t";
			
			//from node
			if(network->getNetworkNode((*ftp).first.first)->getIsDiscreteNode())
			{
				stDev = 0;
				stDevME = 0; 
				ratio = 0;
			}
			else {
				stDev = network->getNetworkNode((*ftp).first.first)->getCtsData()->getStDev();
				nME = nodeME.find(network->getNetworkNode((*ftp).first.first)->getCtsData()->name);
				if(nME != nodeME.end())
				{
					stDevME = nME->second;
				}
				else
				{
					stDevME = 0;
				};
				ratio = stDevME / stDev;
			};

			outEdges << stDev << "\t" << stDevME << "\t" << ratio << "\t";

			//to node
			if(network->getNetworkNode((*ftp).first.second)->getIsDiscreteNode())
			{
				stDev = 0;
				stDevME = 0;
				ratio = 0;
			}
			else {
				stDev = network->getNetworkNode((*ftp).first.second)->getCtsData()->getStDev();
				nME = nodeME.find(network->getNetworkNode((*ftp).first.second)->getCtsData()->name);
				if(nME != nodeME.end())
				{
					stDevME = nME->second;
				}
				else
				{
					stDevME = 0;
				};
			};
			ratio = stDevME / stDev;
			outEdges << stDev << "\t" << stDevME << "\t" << ratio << "\t" << ftp->second <<"\n";
		};
	};

	//add missing edges...
	map<unsigned int, Node *> networkNodes = network->getNetworkNodes();
	map<unsigned int, Node *> parents;
	set<pair<unsigned int, unsigned int> >::const_iterator ope;

	for(map<unsigned int, Node *>::const_iterator nn = networkNodes.begin(); nn != networkNodes.end(); ++nn)
	{
		parents = nn->second->getParents();
		for(map<unsigned int, Node *>::const_iterator p = parents.begin(); p != parents.end(); ++p)
		{
			if(outPuttedEdges.find(make_pair(p->first, nn->first)) == outPuttedEdges.end())
			{
				outEdges << p->second->getDisplayName() << "\t" << nn->second->getDisplayName() << "\t";
				
				//from node
				if(p->second->getIsDiscreteNode())
				{
					stDev = 0;
					stDevME = 0;
					ratio = 0;
				}
				else {
					stDev = p->second->getCtsData()->getStDev();
					nME = nodeME.find(p->second->getCtsData()->name);
					if(nME != nodeME.end())
					{
						stDevME = nME->second;
					}
					else
					{
						stDevME = 0;
					};
				};
				ratio = stDevME / stDev;			
				outEdges << stDev << "\t" << stDevME << "\t" << ratio << "\t";

				//to node
				if(nn->second->getIsDiscreteNode())
				{
					stDev = 0;
					stDevME = 0;
					ratio = 0;
				}
				else {
					stDev = nn->second->getCtsData()->getStDev();
					nME = nodeME.find(nn->second->getCtsData()->name);
					if(nME != nodeME.end())
					{
						stDevME = nME->second;
					}
					else
					{
						stDevME = 0;
					};
				};
				ratio = stDevME / stDev;
				outEdges << stDev << "\t" << stDevME << "\t" << ratio << "\t0\n";
			};
		};
	};

	outEdges.close();
};


//! Outputs R code for plotting the average network.
void AverageNetworksTask::outputRGraph()
{
	string rcodeFileName = igraphPrefix + ".R";

	//output R code for igraph
	ofstream rcode(rcodeFileName.c_str());
	
	rcode << "#threshold, an arc must be greater than the threshold to be plotted\n";
	rcode << "threshold<-" << arcThreshold <<"\n";
	rcode << "plotThresholdEst<-TRUE\n\n";

	rcode << "#load igraph library, http://igraph.org/r/\n";
	rcode << "library(igraph)\n\n";
 
	rcode << "#load average network graph\n";
	rcode << "aveGraph<-read.table(\""<< filename <<"\", header=TRUE, stringsAsFactors=FALSE)\n\n";

	rcode << "#plot arc strength versus cumulative number of arcs with strength <= arc strength\n";
	rcode << "if(plotThresholdEst) {\n";
	rcode << "png(filename=\"" << igraphPrefix << "-thresholdEst.png\", width=600, height=600)\n";
	rcode << "y<-c()\n";
	rcode << "for(stren in aveGraph$strength) y<-append(y, sum(aveGraph$strength <= stren))\n";	
	rcode << "plot.stepfun(aveGraph$strength, xlab=\"arc strength\", ylab=\"cumulative distribution function\", verticals=FALSE, xlim=c(0,1), pch=19, main=\"\")\n";
	rcode << "abline(v=threshold, lty=2)\n";
	rcode << "dev.off()\n";
	rcode << "}\n\n";

	rcode << "#create node and edge tables for igraph\n";
	rcode << "#map node names to numbers\n";
	rcode << "nodeList<-as.numeric(as.factor(c(aveGraph$from, aveGraph$to)))\n";
	rcode << "noArcs<-length(aveGraph$from)\n";
	rcode << "fromNum<-nodeList[1:noArcs]\n";
	rcode << "toNum<-nodeList[(noArcs+1):(2*noArcs)]\n";
	rcode << "nodes1<-as.data.frame(cbind(fromNum, aveGraph$from, aveGraph$type1))\n";
	rcode << "colnames(nodes1)<-c(\"id\", \"name\", \"type\")\n";
	rcode << "nodes2<-as.data.frame(cbind(toNum, aveGraph$to, aveGraph$type2))\n";
	rcode << "colnames(nodes2)<-c(\"id\", \"name\", \"type\")\n";
	rcode << "nodes<-unique(rbind(nodes1, nodes2))\n";
	rcode << "edges<-as.data.frame(cbind(fromNum, toNum, aveGraph$strength, aveGraph$direction))\n";
	rcode << "colnames(edges)<-c(\"from\", \"to\", \"strength\", \"direction\")\n\n";

	rcode << "#apply threshold for plotting arc/edge\n";
	rcode << "edges<-edges[edges$strength > threshold,]\n\n";

	rcode << "#create graph\n";
	rcode << "graph<-graph_from_data_frame(edges, directed = TRUE, vertices = nodes)\n\n";

	rcode << "#plot the network and output png file, edit style as required\n\n";
	
	rcode << "#style for continuous nodes\n";
	rcode << "shape<-rep(\"circle\", length(nodes$type))\n";
	rcode << "vcolor<-rep(\"#eeeeee\", length(nodes$type))\n";
	rcode << "vsize<-rep(25, length(nodes$type))\n";
	rcode << "color<-rep(\"black\", length(nodes$type))\n\n";

	rcode << "#style for discrete nodes\n";
	rcode << "shape[nodes$type==\"d\"]<-\"rectangle\"\n";
	rcode << "vcolor[nodes$type==\"d\"]<-\"#111111\"\n";
	rcode << "vsize[nodes$type==\"d\"]<-20\n";
	rcode << "color[nodes$type==\"d\"]<-\"white\"\n\n";

	rcode << "#style for factor nodes\n";
	rcode << "shape[nodes$type==\"f\"]<-\"rectangle\"\n";
	rcode << "vcolor[nodes$type==\"f\"]<-\"#eeeeee\"\n";
	rcode << "vsize[nodes$type==\"f\"]<-20\n";
	rcode << "color[nodes$type==\"f\"]<-\"black\"\n\n";

	rcode << "#edge widths for significances\n";
	rcode << "minWidth<-0.3\n";
	rcode << "maxWidth<-10\n";
	rcode << "edgeMax<-max(edges$strength)\n";
	rcode << "edgeMin<-min(edges$strength)\n";
	rcode << "widths<-((edges$strength-edgeMin)/(edgeMax-edgeMin))*(maxWidth - minWidth) + minWidth\n";
	rcode << "styles<-rep(1, length(widths))\n\n";

	rcode << "edge.labels<-ifelse(edges$direction==1, round(edges$strength,2), paste(round(edges$strength,2),\"(\",round(edges$direction,2),\")\",sep=\"\"))\n\n";

	rcode << "#plot to a png file\n";
	rcode << "png(filename=\"" << igraphPrefix << ".png\", width=800, height=800)\n\n";

	rcode << "plot(graph, vertex.shape=shape, vertex.size=vsize, vertex.color=vcolor, vertex.label.color=color, ";
	rcode << "edge.width=widths, edge.lty=styles, edge.color=\"black\", edge.arrow.size=1.5, edge.label = edge.labels, edge.label.cex=1.5, edge.label.color=\"red\")\n\n";

	rcode << "#finish png file\n";
	rcode << "dev.off()\n\n";

	rcode.close();

};

//! Estimates the arc/edge threshold for the average network.
void AverageNetworksTask::estimateArcThreshold()
{
	map<double, unsigned int> cumulativeArcStrengths;
	unsigned int cumCount = 0;
	double arcStrength = 0;

	for(list<double>::const_reverse_iterator oas = orderedArcStrengths.rbegin(); oas != orderedArcStrengths.rend(); )
	{
	
		arcStrength = *oas;
		
		do
		{			
			++oas;
			++cumCount;

		}while(*oas <= arcStrength && oas != orderedArcStrengths.rend());

		cumulativeArcStrengths[arcStrength] = cumCount;

	};

	//page 152, Bayes network book, Scutari and Denis
	unsigned int noArcsInFinalNet = 1;
	double bestThreshold = 0; //the threshold that gives minimum L1
	double L1score;
	double bestL1score = (double)cumCount; //upper bound for L1
	double k = (double)cumCount; //no of arcs
	double prevStrength, prevCum;
	double propArcs; //, cum;
	map<double, unsigned int>::const_iterator cas;

	for(map<double, unsigned int>::const_iterator casTh = cumulativeArcStrengths.begin(); casTh != cumulativeArcStrengths.end(); ++casTh) 
	{
		cas = cumulativeArcStrengths.begin();
		propArcs = (double)(casTh->second)/k;

		L1score = propArcs * cas->first;

		prevCum = (double)(cas->second)/k;
		prevStrength = cas->first;
		++cas;

		while(cas != cumulativeArcStrengths.end()) 
		{
			if(prevCum > propArcs) L1score += (prevCum - propArcs) * (cas->first - prevStrength);
			else L1score += (propArcs - prevCum) * (cas->first - prevStrength);
			
			prevCum = (double)(cas->second)/k;
			prevStrength = cas->first;
			++cas;
		};

		if(L1score < bestL1score)
		{
			bestL1score = L1score;
			bestThreshold = casTh->first;
		};
	
	};

	arcThreshold = bestThreshold;
};

//! Sets the final average network using the arc strength threshold.
void AverageNetworksTask::setFinalNetwork()
{

	if(!arcThresholdSet) estimateArcThreshold();

	//remove all existing edges firstly
	network->removeAllEdges();

	string node1, node2;
	unsigned int spacePos;

	for(map<string, pair<double, double> >::iterator tas = arcStrengths.begin(); tas != arcStrengths.end(); ++tas)
	{
		if(tas->second.first > arcThreshold)
		{
			//add edge to the network
			spacePos = (unsigned int)tas->first.find_first_of('\t');
			node1 = tas->first.substr(0, spacePos);
			node2 = tas->first.substr(spacePos + 1);

			if(tas->second.second < 0.5)  network->addEdgeInit(node2, node1);
			else network->addEdgeInit(node1, node2);
		};
	};

	//get rid of bootstrap network and node data which cannot be used by other tasks
	if(!useNetworkScoreMethod)
	{
		taskCentral->removeNetwork(name);
		bootstrapNetwork->deleteAllNetworkNodeData();
		delete bootstrapNetwork;
	};
};

//! Updates arc strengths for average network.
void AverageNetworksTask::updateArcStrengths(Network * net)
{
	list<string> eqNetworks;
	if(freeClearMemory) map<string, pair<double, double> >().swap(networkArcStrengths); else networkArcStrengths.clear();

	string preNodeName = name + "-";

	if(!useEquivNets)
	{
		updateArcStrengthsNetwork(net->getNetworkString(0), preNodeName);
	}
	else
	{
		//get list of equivalent networks for the current network
		unsigned int noEqivNets = net->outputEquivalentNetworks(true, eqNetworks, preNodeName);

		//loop thro' equiv networks and calc and record arc strengths
		for(list<string>::const_iterator en = eqNetworks.begin(); en != eqNetworks.end(); ++en)
		{
			updateArcStrengthsNetwork(*en, preNodeName);
		};

		double noNets = (double)eqNetworks.size();

		//scale to between 0 and 1 for this PDAG
		for(map<string, pair<double, double> >::iterator nas = networkArcStrengths.begin(); nas != networkArcStrengths.end(); ++nas)
		{
			nas->second.first /= noNets;
			nas->second.second /= noNets;
		};
	};

	updateArcStrengthsTotal();
};

//! Updates arc strengths for given network.
void AverageNetworksTask::updateArcStrengthsNetwork(const string & networkStr, const string & preNodeName)
{
	if(useNetworkWeightMethod)
	{
		bootstrapNetwork->removeAllEdges();
		bootstrapNetwork->setNetwork(networkStr, preNodeName); 
		bootstrapNetwork->calcEdgeSignifs();
	};

	string networkString = networkStr;
	map<string, pair<double, double> >::iterator nas;

	bool foundNode;
	string nodeParentsStr;	
	string nodeStr;
	string parentsStr, aParent;
	string arcName;
	unsigned int pos1, pos2, posBar, posColon;
	bool node1ToNode2;
	double strength, direction, prob1to2 = 0, prob2to1 = 0;

	do{
		foundNode = false;
		pos1 = (unsigned int)networkString.find_first_of('[');
		pos2 = (unsigned int)networkString.find_first_of(']');

		if(pos1 < networkString.length() && pos2 < networkString.length() && pos1 + 1 < pos2)
		{
			foundNode = true;
			nodeParentsStr = networkString.substr(pos1+1, pos2-pos1-1);

			//extract node name and parents, so A|B:C gives A with B and C as parents
			posBar = (unsigned int)nodeParentsStr.find_first_of('|');

			//add node
			if(posBar < nodeParentsStr.length()) nodeStr = nodeParentsStr.substr(0, posBar);
			else nodeStr = nodeParentsStr;
				
			//check if parents exist
			if(posBar < nodeParentsStr.length())
			{
				parentsStr = nodeParentsStr.substr(posBar+1);

				do{
					if(useNetworkScoreMethod)
					{
						prob1to2 = 0;
						prob2to1 = 0;
					};

					posColon = (unsigned int)parentsStr.find_first_of(':');

					if(posColon < parentsStr.length()) aParent = parentsStr.substr(0, posColon);
					else aParent = parentsStr;

					//add edge					
					if(nodeStr < aParent)
					{
						arcName = nodeStr + "\t" + aParent;
						node1ToNode2 = false;
					}
					else
					{
						arcName = aParent + "\t" + nodeStr;	
						node1ToNode2 = true;
					};

					nas = networkArcStrengths.find(arcName);
					if(nas != networkArcStrengths.end())
					{
						if(useNetworkScoreMethod)
						{
							if(node1ToNode2) nas->second.first += networkProb; //prob(nd1 --> nd2)
							else nas->second.second += networkProb; //prob(nd1 <-- nd2)
						}
						else if(useNetworkWeightMethod)
						{ 		
							nas->second.first += bootstrapNetwork->getEdgeProb(preNodeName+aParent, preNodeName+nodeStr); //strength
							if(node1ToNode2) nas->second.second++; //direction
						}
						else 
						{
							nas->second.first++; //strength
							if(node1ToNode2) nas->second.second++; //direction
						};				

						
					}
					else
					{
						if(useNetworkScoreMethod)							
						{
							if(node1ToNode2) prob1to2 = networkProb; //prob(nd1 --> nd2)
							else prob2to1 = networkProb; //prob(nd1 <-- nd2)

							networkArcStrengths[arcName] = make_pair(prob1to2, prob2to1);
						}
						else if(useNetworkWeightMethod)
						{ 						
							strength = bootstrapNetwork->getEdgeProb(preNodeName+aParent, preNodeName+nodeStr);
							if(node1ToNode2) direction = 1;
							else direction = 0;
							networkArcStrengths[arcName] = make_pair(strength, direction);
						}
						else
						{
							strength = 1;
							if(node1ToNode2) direction = 1;
							else direction = 0;
							networkArcStrengths[arcName] = make_pair(strength, direction);
						};
								
					};

					if(posColon < parentsStr.length())
					{
						//chop off previous found parent node
						parentsStr = parentsStr.substr(posColon+1);
					}
					else 
						parentsStr = "";

				}while(parentsStr != "");

			};

			//chop off previous found node
			networkString = networkString.substr(pos2+1);
		};

	}while(foundNode);

};

//! Updates total count of arc strengths
void AverageNetworksTask::updateArcStrengthsTotal()
{
	double strength, direction;
	map<string, pair<double, double> >::iterator tas; //total arc strengths
	map<string, double>::iterator tac; //total arc counts

	for(map<string, pair<double, double> >::const_iterator nas = networkArcStrengths.begin(); nas != networkArcStrengths.end(); ++nas)
	{
		tas = arcStrengths.find(nas->first);
		if(tas != arcStrengths.end())
		{
			tas->second.first += nas->second.first; 
			tas->second.second += nas->second.second;	
		}
		else
		{
			strength = nas->second.first; 
			direction = nas->second.second;
			arcStrengths[nas->first] = make_pair(strength, direction);	
		};

		tac = arcCounts.find(nas->first);

		if(tac != arcCounts.end())
		{
			tac->second++;
		}
		else
		{
			arcCounts[nas->first] = 1;
		};
	};

};

//! Updates arc strengths for average network using network score method.
void AverageNetworksTask::updateArcStrengthsScoreMethod()
{
	if(freeClearMemory) map<string, pair<double, double> >().swap(networkArcStrengths); else networkArcStrengths.clear();

	if(!offSetSet)
	{
		offSet = -network->calcScore();
		offSetSet = true;
	};

	double netScore = network->calcScore();

	if(netScore*0 == 0)
	{
		networkProb = exp(network->calcScore() + offSet);
		totalProb += networkProb;
	}
	else
	{
		networkProb = 0;
	};

	updateArcStrengthsNetwork(network->getNetworkString(0));		

	updateArcStrengthsTotal();
};

//! Set the offset for the Score method - do a search to find best network
void AverageNetworksTask::setOffSetScoreMethod()
{

	//setup search object to do the searches for the bootstrapping
	SearchNetworkModelsTask searchNetworkModelsTask;
	searchNetworkModelsTask.setNetworkName(networkName);
	searchNetworkModelsTask.setRandomRestarts(randomRestarts);
	searchNetworkModelsTask.setJitterRestarts(jitterRestarts);
	searchNetworkModelsTask.setTaskCentral(taskCentral); 
	searchNetworkModelsTask.initialiseTask();
	searchNetworkModelsTask.doTask();

	offSet = -network->calcScore();
	offSetSet = true;

};

//! Do the network score method as preposed by So-Youn
void AverageNetworksTask::doNetworkScoreMethod()
{
	//choose network for analyses
	if(networkName != "") network =	taskCentral->getNetwork(networkName);
	else network = taskCentral->getLatestNetwork(networkName);
	
	unsigned int oldScoreType = network->getScoreType();
	network->setScoreType(1); //BIC

	offSetSet = false;

	setOffSetScoreMethod();

	network->calculateAllScoresSetup();
	network->cacheAllNodes();
	
	bool updated = true;	
	double score;	
    noNetworksEval = 0;
	string PDAG;

	while(updated)
	{	
		//evaluate network
		if(!network->hasLoop())
		{
			noNetworksEval++;

			score = network->getScoreFromCache();

			updateArcStrengthsScoreMethod();

		};

		network->updateNextNetwork(updated);
	};

	outputAverageNetworkFile();

	//return network to previous score type
	network->setScoreType(oldScoreType);
};

//! Initialises the compare networks task.
void CompareNetworksTask::initialiseTask()
{
		
		
};

//! Outputs task header for the compare networks task.
void CompareNetworksTask::outputTaskHeader()
{
	outputTaskName();
	out("Comparing networks\n");	
};

//! Outputs task details for the compare networks task.
void CompareNetworksTask::outputTaskDetails()
{
		
		
};

//! Compares networks.
void CompareNetworksTask::doTask()
{
	//check all of the Nodes have the same amount of data
	//taskCentral->checkData(); //do not update missing data here
	//taskCentral->addNetworkMissingData(network->updateNetworkMissingData());	
};

//! Initialises the simulate network data task.
void SimulateNetworkDataTask::initialiseTask()
{
	
		
};

//! Outputs task header for the simulate network data task.
void SimulateNetworkDataTask::outputTaskHeader()
{
	outputTaskName();
	out("Simulating network data\n");	
};

//! Outputs task details for the simulate network data task.
void SimulateNetworkDataTask::outputTaskDetails()
{
	
	if(networkName!="") {out("Data simulation network given by network: "); out(networkName); out("\n");};
	out("Number of simulations: "); out(noSims); out("\n");
	if(paraFilename!="") {out("Parameter file name: "); out(paraFilename); out("\n");};
	
	out("Network: "); out(name); out("\n");
	out("Network type: "); out(simNetwork->getNetworkType()); out("\n");
	if(simNetwork->getNetworkType() != "deal") {out("Network score type: "); out(simNetwork->getScoreTypeName()); out("\n"); };
	if(simNetwork->getScoreFixName() != "none") { out("Network score fix: "); out(simNetwork->getScoreFixName()); out("\n"); };
	out("Total number of nodes: "); out(totalNoDataNodes + totalDisNodes + totalFactorNodes + totalCtsNodes); out(" (Discrete: ");  out(totalDisNodes); out(" | Factor: ");  out(totalFactorNodes); out(" | Continuous: "); out(totalCtsNodes);
			if(totalNoDataNodes != 0) {out(" | No data: "); out(totalNoDataNodes);} out(")"); out("\n");
	out("Total number of edges: "); out(totalEdges); out("\n");
	out("Network Structure: "); out(simNetwork->getNetworkString()); out("\n");
	if(simNetwork->hasLoop()) {out("The network contains a loop!\n");};
	if(simNetwork->getScoreTypeName() == "bayes" || simNetwork->getNetworkType() == "deal") { out("Imaginary sample size: "); out(simNetwork->getImaginarySampleSize()); out("\n"); };
	if(simNetwork->hasNodesWithNoData()) out("The network has nodes with no data\n");
	else
	{
		NetworkMissingData * networkMissingData = simNetwork->updateNetworkMissingData();
		pair<unsigned int, unsigned int> missingNotMissing = simNetwork->getNoMissingNotMissing();
		out("Total data at each node: "); out(missingNotMissing.second); out("\n");
		out("Missing data at each node: "); out(missingNotMissing.first); out("\n");
		delete networkMissingData;
	};
};
	
//! Simulates network data.
void SimulateNetworkDataTask::doTask()
{
	if(networkName != "" && paraFilename != "")
	{
		exitErr("Choose either a network or a network parameter file to simulate network node data!");
	};

	//create a network with empty data 
	if(paraFilename != "") createSimNetworkFromParaFile();
	else createSimNetworkFromNetwork();

	set<pair<unsigned int, unsigned int> > allowedEdgeTypesFileNos; // fileNo1, fileNo2
	set<pair<unsigned int, unsigned int> > notAllowedEdgeTypesFileNos;
	
	simNetwork->setUpAllowedEdgeTypes(allowedEdgeTypesFileNos, notAllowedEdgeTypesFileNos);

	//simulate the data
	simNetwork->simulateData(noSims, snpDataIntegers);

	//set some missing data
	simNetwork->setSomeMissingData(missingFirst, missingLast, missingProb);

};

//! Simulates network data using current network for parameter estimates.
void SimulateNetworkDataTask::createSimNetworkFromNetwork()
{
	
	//choose network to perform analysis
	if(networkName != "") network =	taskCentral->getNetwork(networkName);
	else network = taskCentral->getLatestNetwork(networkName);

	if(network->getNetworkType() == "deal")
	{
		exitErr("A deal network may not be used to simulete data!");
	};

	simNetwork = (NetworkBNLearn*)taskCentral->createSimNetworkFromNetwork(network, name, createDifferentNodeData);

	simNetwork->getNumberNodesAndEdges(totalDisNodes, totalFactorNodes, totalCtsNodes, totalNoDataNodes, totalEdges);
	simNetwork->setImaginarySampleSize(imaginarySampleSize);
};

//! Simulates network data using given network for parameter estimates and returns network with different data.
Network * TaskCentral::createSimNetworkFromNetwork(Network * network, const string & taskName, const bool & createDifferentNodeData)
{
	Network * simNetwork;

	bool emptyDataNetwork  = network->hasNodesWithNoData();

	//addNetworkMissingData(network->updateNetworkMissingData());

	//check all of the Nodes have the same amount of data
	checkData();

	if(!emptyDataNetwork) addNetworkMissingData(network->updateNetworkMissingData());
	
	//calculate the posterior in order to ensure everything is set up in nodes for simulating the data
	network->calculatePriorsAndPosteriorsBNL();

	network->calcScore();

	if(network->getNetworkType() == "deal")
	{
		simNetwork = new NetworkDeal(getAllNodeData(), network->getScoreType());
	}
	else
	{
		simNetwork = new NetworkBNLearn(getAllNodeData(), network->getScoreType(), network->getScoreFix());
	};

	//set the same network score type
	simNetwork->setScoreType(network->getScoreType());
	
	//add nodes
	unsigned int nd = 1;
	unsigned int noCurNetworkNodes = 0;
	unsigned int noNodes = network->getNoNodes();
	AllNodeData * allNodeData = getAllNodeData();
	unsigned int noAllDataNodes = allNodeData->getNoNodeData();
	Node * aNode;
	Node * simNode;
	Data * simData;
	string dataName;
	string preNodeName = taskName + "-";
	string displayName;

	unsigned int nodeNumber;

	unsigned int fileNo;
	unsigned int noSimLevels;

	while(nd <= noAllDataNodes && noCurNetworkNodes < noNodes)
	{
		
		if(network->nodeExists(nd)) 
		{
			aNode = network->getNetworkNode(nd);
			fileNo = aNode->getFileNo(); //use same file no so that black/white lists types can be copied

			if(createDifferentNodeData)
			{
				dataName = preNodeName + aNode->getName();

				//add new empty data
				if(aNode->getNodeType() == "c")
				{
					simData = new CtsData(dataName, fileNo);
				}
				else if(aNode->getNodeType()=="d")
				{
					simData = new DiscreteData(dataName, fileNo);					
				}
				else if(aNode->getSimNodeType() == 0) simData = new CtsData(dataName, fileNo);
				else
				{
						simData = new DiscreteData(dataName, fileNo);						
				};				

				allNodeData->addNodeData(simData);
			}
			else
			{
				dataName = aNode->getName();

				if(aNode->getNodeType() == "c")
				{
					simData = aNode->getCtsData();
				}
				else if(aNode->getNodeType()=="d")
				{
					simData = aNode->getDiscreteData();					
				}
				else
				{
					if(aNode->getSimNodeType() == 0) simData = new CtsData(dataName, fileNo); //for empty data nodes
					else
					{
						simData = new DiscreteData(dataName, fileNo);			
					};

					allNodeData->replaceNodeData(simData);
				};

				simData->clearData(); //remove any previous data
			};

			simData->isSNPData = aNode->getIsSNPNode();			
			
			//add node to the network
			simNetwork->addNodeInit(dataName);

			//get the added node and set some things
			nodeNumber = allNodeData->getNodeDataNumber(dataName);
			simNode = simNetwork->getNetworkNode(nodeNumber);

			displayName = aNode->getName();
			simNode->setDisplayName(displayName);	
			
			simNode->setIsFactorNode(aNode->getIsFactorNode());
			simNode->setIsChildFactorNode(aNode->getIsFactorChildNode());

			if(!emptyDataNetwork)
			{
				simNode->simDataCopyDisLevels(aNode);
			}
			else
			{
				noSimLevels = aNode->getSimNoLevels();
				simNode->setSimNoLevels(noSimLevels);
			};

			++noCurNetworkNodes;
		};
		++nd;
	};

	//add the same edges as the prev network		
	nd = 1;
	noCurNetworkNodes = 0;
	Node * origNode;
	Node * copyingNode;
	string nodeName;

	while(nd <= noAllDataNodes && noCurNetworkNodes < noNodes)
	{
		if(network->nodeExists(nd))
		{
			origNode = network->getNetworkNode(nd);
			if(createDifferentNodeData) nodeName = preNodeName + origNode->getName(); //get corresponding node in sim network
			else nodeName = origNode->getName();
			
			nodeNumber = allNodeData->getNodeDataNumber(nodeName);
			copyingNode = simNetwork->getNetworkNode(nodeNumber);

			origNode->simDataCopyParents(taskName, copyingNode, simNetwork, createDifferentNodeData);
			if(!emptyDataNetwork) origNode->simDataCopyParas(taskName, copyingNode, network, simNetwork, createDifferentNodeData);

			++noCurNetworkNodes;
		};
		++nd;
	};

	//set default parameters
	if(emptyDataNetwork)
	{
		simNetwork->setDefaultParameters();
	};

	//create new missing data object for this sim network
	addNetworkMissingData(simNetwork->updateNetworkMissingData());

	if(createDifferentNodeData) simNetwork->updateBlackWhiteDifferentData(network, preNodeName); 
	else simNetwork->updateBlackWhite(network); 

	//copy edge costs
	map<pair<unsigned int, unsigned int>, double> costEdges = network->getCostEdges();
	map<pair<unsigned int, unsigned int>, double> costEdgeTypes = network->getCostEdgeTypes();;
	simNetwork->setCostEdges(costEdges, costEdgeTypes);

	addNetwork(taskName, simNetwork); //add network with name of this task

	return simNetwork;
};

//! Creates network to simulate data given parameter file.
void SimulateNetworkDataTask::createSimNetworkFromParaFile()
{

	//create network for simulating data
	simNetwork = new NetworkBNLearn(taskCentral->getAllNodeData(), getNetworkType(score), getNetworkScoreFix(scoreFix));

	//add nodes
	addNodesToSimNetwork();

	//add parents for each node
	addParentsToSimNetwork();

	//find the number of discrete levels for discrete nodes
	findNoLevelsFromParasFile();

	//add parents and parameters for each node
	addParasToSimNetwork();

	//create new missing data object for this sim network	
	taskCentral->addNetworkMissingData(simNetwork->updateNetworkMissingDataSimNet(noSims));

	simNetwork->getNumberNodesAndEdges(totalDisNodes, totalFactorNodes, totalCtsNodes, totalNoDataNodes, totalEdges);

	simNetwork->setImaginarySampleSize(imaginarySampleSize);

	taskCentral->addNetwork(name, simNetwork); //add network with name of this task
};

//! Adds nodes to the simulate data network.
void SimulateNetworkDataTask::addNodesToSimNetwork()
{
	ifstream readParaFile;
	readParaFile.open(paraFilename.c_str());
	if(!readParaFile.is_open())
	{
		string mess = "Cannot read parameter file: " + paraFilename + "!\n";
		exitErr(mess);
	};
	
	string word1 = "", word2 = "", word3 = "", word4 = "";
	string dataName = "";
	bool isDiscrete = false;
	bool isSNPNode = false;
	Data * simData;
	Node * simNode;
	bool foundNode;
	AllNodeData * allNodeData = taskCentral->getAllNodeData();
	unsigned int ctsFileNo = taskCentral->getFileNo();//need to add different file nos in order to exclude dis-->cts edges (even though data is not from a file)
	taskCentral->nextFileNo();
	unsigned int disFileNo = taskCentral->getFileNo();
	unsigned int nodeNumber;

	do{

		isDiscrete = false;
		isSNPNode = false;
		foundNode = true;

		//read in next 4 words
		word1 = word2;
		word2 = word3;
		word3 = word4;
		readParaFile >> word4;
		
		if(word2 == "CONTINUOUS" && word3 == "NODE:")
		{
			dataName = word4;			
		}
		else if(word1 == "CONTINUOUS" && word2 == "SNP" && word3 == "NODE:")
		{
			dataName = word4;			
			isSNPNode = true;
		}
		else if(word2 == "DISCRETE" && word3 == "NODE:")
		{
			dataName = word4;
			isDiscrete = true;
		}
		else if(word1 == "DISCRETE" && word2 == "SNP" && word3 == "NODE:")
		{
			dataName = word4;
			isDiscrete = true;
			isSNPNode = true;
		}
		else
		{
			foundNode = false;
		};
		
		if(foundNode)
		{
			//add new empty data
			if(!isDiscrete)
			{
				simData = new CtsData(dataName, ctsFileNo);
			}
			else 
			{
				simData = new DiscreteData(dataName, disFileNo);			
			};

			simData->isSNPData = isSNPNode;

			allNodeData->addNodeData(simData);

			//add node to the network
			simNetwork->addNodeInit(dataName);

			//get the added node and set some things
			nodeNumber = allNodeData->getNodeDataNumber(dataName);
			simNode = simNetwork->getNetworkNode(nodeNumber);
			simNode->setDisplayName(dataName);				
		};

	}while(!readParaFile.eof());

	readParaFile.close();
};

//! Removes any dodgey returns that may be from windows.
void removeReturns(string & s)
{
	if(!s.empty() && s.substr(s.size() - 1) == "\r") s = s.substr(0, s.size() - 1);
};

//! Trims string.
void trim(string & str, const string & whitespace)
{
    size_t strBegin = str.find_first_not_of(whitespace);

    if(strBegin == string::npos) return; // no content

    size_t strEnd = str.find_last_not_of(whitespace);
    size_t strRange = strEnd - strBegin + 1;

    str = str.substr(strBegin, strRange);
}

//! Get parents names from string parent1:parent2:parent3.
list<string> splitString(const string & strParents, string delimiter)
{
	list<string> parents;
	string s = strParents;	
	size_t pos = 0;
	string token;

	while((pos = s.find(delimiter)) != string::npos)
	{
		token = s.substr(0, pos);
		parents.push_back(token);
		s.erase(0, pos + delimiter.length());
	};

	if(s.length() > 0) parents.push_back(s);

	//remove any empty strings
	list<string>::iterator pa = parents.begin();
	while(pa != parents.end())
	{
		if(*pa == "" || *pa == "\r" || *pa == "\n" || *pa == "\t")
		{
			parents.erase(pa);
			pa = parents.begin();
		}
		else
			++pa;
	};

	//remove any spurious returns and any white space	
	for(pa = parents.begin(); pa != parents.end(); ++pa)
	{		
		removeReturns(*pa);
		trim(*pa);
	};	

	return parents;
};

//! Adds the parents to the simulation network.
void SimulateNetworkDataTask::addParentsToSimNetwork()
{
	//add discrete parents firstly
	ifstream readParaFile;
	readParaFile.open(paraFilename.c_str());
	if(!readParaFile.is_open())
	{
		string mess = "Cannot read parameter file: " + paraFilename + "!\n";
		exitErr(mess);
	};
	
	string word1 = "", word2 = "", word3 = "", word4 = "";
	Node * simNode;
	Node * parentSimNode;

	string nodeName;		
	list<string> parents;
	bool foundNode;
	AllNodeData * allNodeData = taskCentral->getAllNodeData();
	unsigned int nodeNumber;

	do{
	
		foundNode = true;

		//read in next 4 words
		word1 = word2;
		word2 = word3;
		word3 = word4;
		readParaFile >> word4;
		
		if(word2 == "DISCRETE" && word3 == "PARENTS:" )
		{
			nodeName = word1;
			parents = splitString(word4);
			nodesWithDiscreteParents.insert(nodeName);
		}
		else
		{
			foundNode = false;
		};
		
		if(foundNode)
		{
			
			//get the added node and set some things
			nodeNumber = allNodeData->getNodeDataNumber(nodeName);
			simNode = simNetwork->getNetworkNode(nodeNumber);
			
			//add parents
			for(list<string>::const_iterator pa = parents.begin(); pa != parents.end(); ++pa)
			{
				nodeNumber = allNodeData->getNodeDataNumber(*pa);
				parentSimNode = simNetwork->getNetworkNode(nodeNumber);
				simNetwork->addEdge(parentSimNode, simNode);
			};
		};

	}while(!readParaFile.eof());

	readParaFile.close();

	//find the cts parents next
	readParaFile.open(paraFilename.c_str());
	if(!readParaFile.is_open())
	{
		string mess = "Cannot read parameter file: " + paraFilename + "!\n";
		exitErr(mess);
	};

	string line;
	bool ctsSNPFound;
	bool coeffsFound;
	list<string> coeffs;
	list<string>::const_iterator co;
	string parentName;

	do{

		getline(readParaFile, line);

		ctsSNPFound = false;

		if(line.length() > 17 && line.substr(0, 17) == "CONTINUOUS NODE: ") 		
		{
			nodeName = line.substr(17);
			removeReturns(nodeName);

			ctsSNPFound = true;
		}
		else if((line.length() > 21 && line.substr(0, 21) == "CONTINUOUS SNP NODE: "))
		{
			nodeName = line.substr(21);
			removeReturns(nodeName);
			
			ctsSNPFound = true;
		};

		if(ctsSNPFound)
		{
			coeffsFound = false;
			nodeNumber = allNodeData->getNodeDataNumber(nodeName);
			simNode = simNetwork->getNetworkNode(nodeNumber);

			while(!coeffsFound)
			{
				getline(readParaFile, line);

				coeffs = splitString(line, " ");

				if(!coeffs.empty() && *coeffs.begin() == "Coefficients:") coeffsFound = true;
				
				if(readParaFile.eof())
				{
					outErr("Coefficients not found for continuous node "); outErr(nodeName); outErr("!");
					exitErr("");
				};
			};

			co = coeffs.begin();
			co++;

			while(co != coeffs.end())
			{
				parentName = *co;
				if(parentName.substr(parentName.length()-1,1) == ":")
				{
					parentName = parentName.substr(0, parentName.length()-1);
					nodeNumber = allNodeData->getNodeDataNumber(parentName);
					parentSimNode = simNetwork->getNetworkNode(nodeNumber);
					simNetwork->addEdge(parentSimNode, simNode);
				};

				++co;
			};
		};
	}while(!readParaFile.eof());

	readParaFile.close();

};

//! Returns whether a node had discrete parents or not.
bool SimulateNetworkDataTask::getHasDiscreteParents(const string & nd)
{
	set<string>::const_iterator nwdp = nodesWithDiscreteParents.find(nd);

	if(nwdp != nodesWithDiscreteParents.end()) return true; 

	return false;
};

//! Finds the number of levels for a discrete node.
void SimulateNetworkDataTask::findNoLevelsFromParasFile()
{
	ifstream readParaFile;
	readParaFile.open(paraFilename.c_str());
	if(!readParaFile.is_open())
	{
		string mess = "Cannot read parameter file: " + paraFilename + "!\n";
		exitErr(mess);
	};

	string line;
	string nodeName;
	bool snpFound;
	unsigned int nodeNoLevels = 0;
	list<string> lineWords;
	bool foundLevel;
	string levelName;
	unsigned int nodeNumber;

	Node * simNode;
	AllNodeData * allNodeData = taskCentral->getAllNodeData();
	
	getline(readParaFile, line);

	do
	{
		snpFound = false;

		if(line.length() > 15 && line.substr(0, 15) == "DISCRETE NODE: ") 		
		{
			nodeName = line.substr(15);	
			removeReturns(nodeName);
			snpFound = true;
		}
		else if((line.length() > 19 && line.substr(0, 19) == "DISCRETE SNP NODE: "))
		{
			nodeName = line.substr(19);	
			removeReturns(nodeName);
			snpFound = true;
		};

		if(snpFound)
		{
			nodeNoLevels = 0;
			nodeNumber = allNodeData->getNodeDataNumber(nodeName);
			simNode = simNetwork->getNetworkNode(nodeNumber);
			
			if(simNode->getNoParents() > 0)
			{
				if(!getline(readParaFile, line)) break; //PARENTS:
				if(!getline(readParaFile, line)) break; //parent levels
			};

			foundLevel = true;

			while(foundLevel && getline(readParaFile, line)) 
			{
				lineWords = splitString(line, ": ");
				if(lineWords.size() == 2)
				{
					nodeNoLevels++;
					foundLevel = true;
					//add name of level to discrete node
					levelName = *lineWords.begin();
					simNode->getLevelNo(levelName);
				}
				else
					foundLevel = false;
			};

			noLevels[nodeName] = nodeNoLevels;
		}
		else
		{
			getline(readParaFile, line);
		};

	}while(!readParaFile.eof());

	readParaFile.close();
};

//! Adds parameters for simulating data to a discrete node.
void SimulateNetworkDataTask::addDiscreteSimParas(Node * simNode, ifstream & readParaFile)
{
	unsigned int parentGroupNo = simNode->getParentsGroupNo();
	string nodeName = simNode->getName();
	map<unsigned int, double> probs; //level, prob
	map<unsigned int, double> cumProbs; //level, cum prob
	string levelName;
	string probStr;
	double prob;
	double probTot = 0;

	unsigned int paraLevel;

	map<string, unsigned int>::const_iterator nlv = noLevels.find(nodeName);
	if(nlv == noLevels.end())
	{
		outErr("Problem with node "); outErr(nodeName); outErr(" when parsing parameter file for simulating data!");
		exitErr("");
	};

	for(unsigned int lv = 1; lv <= nlv->second; ++lv)
	{
		readParaFile >> levelName;
		readParaFile >> probStr;

		//remove colon
		if(levelName.substr(levelName.length()-1, 1)==":") levelName = levelName.substr(0, levelName.length()-1);

		if(readParaFile.eof())
		{
			outErr("Problem with node "); outErr(nodeName); outErr(" when parsing parameter file for simulating data!");
			exitErr("");
		};

		prob = atof(probStr.c_str());
		
		paraLevel = simNode->getLevelNo(levelName);
		probs[paraLevel] = prob;
		probTot += prob;
	};

	//set all to probs to 1, if all 0
	if(probTot == 0)
	{
		for(map<unsigned int, double>::iterator pr = probs.begin(); pr != probs.end(); ++pr)
		{
			pr->second = 1;			
		};
		probTot = probs.size();
	};

	//normalise probs
	if(probTot != 1)
	{
		for(map<unsigned int, double>::iterator pr = probs.begin(); pr != probs.end(); ++pr)
		{
			pr->second /= probTot;			
		};
	};

	double cumProb = 0;
	unsigned int lv;

	//set cumalative probs
	for(map<unsigned int, double>::iterator pr = probs.begin(); pr != probs.end(); )
	{
		cumProb += pr->second;
		lv = pr->first;
		++pr;
		if(pr == probs.end()) break; //do not add last one as should be 1 any way
		cumProbs[lv] = cumProb;	
	};
	
	simNode->setSimDataCumLevelProb(parentGroupNo, cumProbs);
};

//! Adds parameters for simulating data to a discrete node.
void SimulateNetworkDataTask::addCtsSimParas(Node * simNode, ifstream & readParaFile, string & word)
{
	unsigned int nodeParentGroupNo = simNode->getNodeAndParentsGroupNo();
	CtsPriLocalDist * ctsParas = new CtsPriLocalDist();
	bool interceptDone = false;
	bool coeffsDone = false;
	bool varianceDone = false;	
	double var;
	bool readNextWord = false;
	string ctsParentName;
	Node * ctsParentNode;
	AllNodeData * allNodeData = taskCentral->getAllNodeData();
	unsigned int nodeNumber;

	readParaFile >> word;

	while((!interceptDone || !coeffsDone || !varianceDone) && !readParaFile.eof())
	{
		
		if(word == "Intercept:")
		{
			readParaFile >> word;			
			ctsParas->intercept = atof(word.c_str());
			interceptDone = true;
			readParaFile >> word;
		}
		else if(word == "Variance:")
		{
			readParaFile >> word;
			var = atof(word.c_str());
			ctsParas->variance = var;
			ctsParas->stDev = sqrt(var);
			varianceDone = true;
			readParaFile >> word;
		}
		else if(word == "Coefficients:")
		{
			readParaFile >> word;			

			while(!readParaFile.eof() && word != "Variance:" && word != "Mean:")
			{
				ctsParentName = word;

				//remove colon
				if(ctsParentName.substr(ctsParentName.length()-1, 1)==":") ctsParentName = ctsParentName.substr(0, ctsParentName.length()-1);

				nodeNumber = allNodeData->getNodeDataNumber(ctsParentName);
				ctsParentNode = simNetwork->getNetworkNode(nodeNumber);

				//add coeff with value
				readParaFile >> word;
				ctsParas->coeffs[ctsParentNode->getNodeID()] = atof(word.c_str()); 

				//see if the next word is another coeff
				readParaFile >> word;
			};

			coeffsDone = true;
		}
		else
		{
			readParaFile >> word;
		};

		
	};

	if(readParaFile.eof() && (!interceptDone || !coeffsDone || !varianceDone))
	{
		outErr("Problem setting the parameters for continuous node "); outErr(simNode->getName()); outErr("!");
		exitErr("");
	};

	simNode->setCtsParas(nodeParentGroupNo, ctsParas);
};

//! Adds the parameters to the simulation network.
void SimulateNetworkDataTask::addParasToSimNetwork()
{
	//add parents firstly
	ifstream readParaFile;
	readParaFile.open(paraFilename.c_str());
	if(!readParaFile.is_open())
	{
		string mess = "Cannot read parameter file: " + paraFilename + "!\n";
		exitErr(mess);
	};
	
	Node * simNode;
	Node * parentSimNode;

	string nodeName;
	string parentString, strParents;
	list<string> discreteParents;
	list<string> discreteParentLevels;
	string word;
	unsigned int nodeNumber;

	bool foundNode;
	bool isDiscrete;
	bool hasDiscreteParents;

	AllNodeData * allNodeData = taskCentral->getAllNodeData();
	
	do{
	
		foundNode = false;
		isDiscrete = false;
		hasDiscreteParents = false;
		discreteParents.clear();
		discreteParentLevels.clear();

		//read in word		
		readParaFile >> word;
		
		if(word == "NODE:")
		{			
			readParaFile >> nodeName;
			foundNode = true;
		};
		
		if(foundNode)
		{
			
			//get the added node and set some things
			nodeNumber = allNodeData->getNodeDataNumber(nodeName);
			simNode = simNetwork->getNetworkNode(nodeNumber);
			
			if(getHasDiscreteParents(nodeName))
			{
				hasDiscreteParents = true;
				readParaFile >> word >> parentString >> strParents;
				if(parentString != "PARENTS:")
				{
					outErr("Problem with node "); outErr(nodeName); outErr(" and setting parameters.\n");
					exitErr("Problem parsing parameter file for simulating data!");
				};

				discreteParents = splitString(strParents);

				//read in levels
				readParaFile >> word;
			};
		
			//read in parameters for each set of parent levels
			do{

				
				if(hasDiscreteParents)
				{
					discreteParentLevels = splitString(word);

					if(discreteParentLevels.size() != discreteParents.size())
					{
						outErr("Number of discrete parent levels do not match the number of discrete parents for node "); outErr(simNode->getName()); outErr("!");
						exitErr("");
					};

					list<string>::const_iterator palv = discreteParentLevels.begin();					

					//set levels of discrete parents
					for(list<string>::const_iterator pa = discreteParents.begin(); pa != discreteParents.end(); ++pa, ++palv)
					{
						nodeNumber = allNodeData->getNodeDataNumber(*pa);
						parentSimNode = simNetwork->getNetworkNode(nodeNumber);
						parentSimNode->setLevel(*palv);
					};
				};

				if(simNode->getIsDiscreteNode())
				{
					addDiscreteSimParas(simNode, readParaFile);
					readParaFile >> word;
				}
				else
				{
					addCtsSimParas(simNode, readParaFile, word);
				};
				
			}while(!readParaFile.eof() && hasDiscreteParents && word != "DISCRETE" && word != "CONTINUOUS");
		};

	}while(!readParaFile.eof());

	readParaFile.close();
};

//! Initialises the task for outputting the network.
void OutputNetworkTask::initialiseTask()
{
	if(jobNo > 0) taskCentral->setStartAndEndIndivsBasedOnJobNo(jobNo, jobTotal, startIndivNo, endIndivNo);		
};

//! Outputs the task header for outputting the network.
void OutputNetworkTask::outputTaskHeader()
{
	outputTaskName();
	out("Outputting network\n");
};

//! Outputs the task details for outputting the network.
void OutputNetworkTask::outputTaskDetails()
{
	out("Network: "); out(networkName); out("\n");
	out("Network Structure: "); out(network->getNetworkString()); out("\n");
	if(filenamePrefix != "")
	{
		out("Network output to igraph files:\n");
		out(" "); out(filenamePrefix); out("-nodes.dat\n");
		out(" "); out(filenamePrefix); out("-edges.dat\n");	
		out("R code to plot network using igraph package: "); out(filenamePrefix); out("-plot.R\n");	
	};

	if(filename != "")
	{
		out("Network output to file: "); out(filename); out("\n");
	};

	if(filename2 != "")
	{
		out("Network output in bnlearn format to file: "); out(filename2); out("\n");
	};

	if(equivFilename != "")
	{
		if(!loopInNetwork)
		{
			if(noEquivNets > 1) {out("A list of "); out(noEquivNets); out(" equivalent networks are output to file: ");} 
			else {out("A list of "); out(noEquivNets); out(" equivalent network is output to file: ");}; 	
		}
		else
		{
			out("The network has a loop in it so a list of equivalent networks could not be output to file: ");
		};

		out(equivFilename); out("\n");
	};

	if(nodeDataFilePrefix!="" && (wasDisData || wasCtsData || wasSNPData || wasSNPDisData || wasSNPCtsData))
	{
		out("Node data output to files:\n"); 
		if(wasDisData) out(" "+nodeDataFilePrefix+"-discrete.dat\n");
		if(wasSNPDisData) out(" "+nodeDataFilePrefix+"-SNPs-discrete.dat\n");
		if(wasCtsData) out(" "+nodeDataFilePrefix+"-cts.dat\n");	
		if(wasSNPCtsData) out(" "+nodeDataFilePrefix+"-SNPs-cts.dat\n");
		if(wasSNPData) out(" "+nodeDataFilePrefix+".bed/.bim/.fam\n");	
		if(jobNo > 0) {out("Job "); out(jobNo); out(" of "); out(jobTotal); out("\n");};
		if(startIndivNo > 0) {out("Individuals "); out(startIndivNo); out(" to "); out(endIndivNo); out("\n");};
	};

};
	
//! Outputs the network.
void OutputNetworkTask::doTask()
{
	//choose network to perform analysis
	if(networkName != "") network =	taskCentral->getNetwork(networkName);
	else network = taskCentral->getLatestNetwork(networkName);

	if(filename == "" && filenamePrefix == "") filename = "network.dat";

	//output in igraph
	if(filenamePrefix != "") outputNetworkIGraph();
	
	//output network
	if(filename != "") outputNetwork();

	//output network in bnlearn format
	if(filename2 != "") outputNetwork2();

	//output network node data
	if(nodeDataFilePrefix != "") outputNodeData();

	//output equivalent networks
	if(equivFilename != "")
	{
		loopInNetwork = network->hasLoop();
		if(!loopInNetwork) outputEquivalentNetworks();
	};
};

//! Outputs network node data to cts, discrete and possibly .bed files
void OutputNetworkTask::outputNodeData()
{
	network->outputNodeData(nodeDataFilePrefix, outputBedFile, wasDisData, wasCtsData, wasSNPDisData, wasSNPCtsData, wasSNPData, startIndivNo, endIndivNo, jobNo);
};

//! Outputs the network in BayesNetty format.
void OutputNetworkTask::outputNetwork()
{
	network->outputNetwork(filename);
};

//! Outputs the network in [A][C][B|A:C] format.
void OutputNetworkTask::outputNetwork2()
{
	ofstream outNet(filename2.c_str());

	outNet << network->getNetworkString(0) << "\n";
	
	outNet.close();
};

//! Outputs equivalent networks to file.
void OutputNetworkTask::outputEquivalentNetworks()
{
	list<string> dummyList; //not used but needed as method has dual purpose

	noEquivNets = network->outputEquivalentNetworks(false, dummyList, equivFilename);
	
};

//! Outputs the network in format for use with igraph R package use.
void OutputNetworkTask::outputNetworkIGraph()
{
	//check all of the Nodes have the same amount of data - needed for edge weights
	taskCentral->checkData();
	taskCentral->addNetworkMissingData(network->updateNetworkMissingData());

	string edgeFileName = filenamePrefix + "-edges.dat";
	string nodesFileName = filenamePrefix + "-nodes.dat";
	string rcodeFileName = filenamePrefix + "-plot.R";

	network->outputNetworkIGraphEdges(edgeFileName);
	network->outputNetworkIGraphNodes(nodesFileName);

	//output R code for igraph
	ofstream rcode(rcodeFileName.c_str());
	
	rcode << "#load igraph library, http://igraph.org/r/\n";
	rcode << "library(igraph)\n\n";
 
	rcode << "#load network graph\n";
	rcode << "nodes<-read.table(\"" << nodesFileName <<"\", header=TRUE)\n";
	rcode << "edges<-read.table(\""<< edgeFileName <<"\", header=TRUE)\n\n";

	rcode << "#create graph\n";
	rcode << "graph<-graph_from_data_frame(edges, directed = TRUE, vertices = nodes)\n\n";

	rcode << "#plot the network and output png file, edit style as required\n\n";
	
	rcode << "#style for continuous nodes\n";
	rcode << "shape<-rep(\"circle\", length(nodes$type))\n";
	rcode << "vcolor<-rep(\"#eeeeee\", length(nodes$type))\n";
	rcode << "vsize<-rep(25, length(nodes$type))\n";
	rcode << "color<-rep(\"black\", length(nodes$type))\n\n";

	rcode << "#style for discrete nodes\n";
	rcode << "shape[nodes$type==\"d\"]<-\"rectangle\"\n";
	rcode << "vcolor[nodes$type==\"d\"]<-\"#111111\"\n";
	rcode << "vsize[nodes$type==\"d\"]<-20\n";
	rcode << "color[nodes$type==\"d\"]<-\"white\"\n\n";

	rcode << "#style for factor nodes\n";
	rcode << "shape[nodes$type==\"f\"]<-\"rectangle\"\n";
	rcode << "vcolor[nodes$type==\"f\"]<-\"#eeeeee\"\n";
	rcode << "vsize[nodes$type==\"f\"]<-20\n";
	rcode << "color[nodes$type==\"f\"]<-\"black\"\n\n";

	rcode << "#edge widths for significances\n";
	rcode << "minWidth<-0.3\n";
	rcode << "maxWidth<-10\n";
	rcode << "edgeMax<-max(edges$chisq)\n";
	rcode << "edgeMin<-min(edges$chisq)\n";
	rcode << "widths<-((edges$chisq-edgeMin)/(edgeMax-edgeMin))*(maxWidth - minWidth) + minWidth\n";
	rcode << "styles<-rep(1, length(widths))\n\n";

	/*
	rcode << "#if the data is imputed draw some edges dashed to suggest possible false positive edges (most likely when the intermediate variable is heavily imputed)\n";
	rcode << "isImputedData<-";
	if(network->getIsImputedData()) rcode << "TRUE"; else rcode << "FALSE";
	rcode<< "\n\n";

	rcode << "existsStrongerIndirectRoute<-function(edgeNo, edges)\n";
	rcode << "{\n";
	rcode << "isOtherRoute<-FALSE\n";
	rcode << "allFroms<-edges[edges$from==edges$from[edgeNo],]\n\n";
 
	rcode << "for(i in 1:length(allFroms$to))\n";
	rcode << "{\n";
	rcode << "nextEdges<-edges[allFroms$to[i]==edges$from,]\n";
	rcode << "num<-which(edges$to[edgeNo] == nextEdges$to)\n"; 
	rcode << "if(length(num)==1 && edges[edgeNo,3]<allFroms[i, 3] && edges[edgeNo,3]<nextEdges[num, 3]) isOtherRoute<-TRUE\n"; 
	rcode << "}\n\n";

	rcode << "isOtherRoute\n";
	rcode << "}\n\n";

	rcode << "getAllDashed<-function(edges)\n";
	rcode << "{\n";
	rcode << "results<-rep(FALSE, length(edges$from))\n";
	rcode << "for(i in 1:length(edges$from)) results[i]<-existsStrongerIndirectRoute(i, edges)\n";  
	rcode << "results\n";
	rcode << "}\n\n";

	rcode << "#set some edges to be dashed using the two above functions\n";
	rcode << "if(isImputedData) styles[getAllDashed(edges)]<-2\n\n";
	*/

	rcode << "#plot to a png file\n";
	rcode << "png(filename=\"" << filenamePrefix << ".png\", width=800, height=800)\n\n";

	rcode << "plot(graph, vertex.shape=shape, vertex.size=vsize, vertex.color=vcolor, vertex.label.color=color, ";
	rcode << "edge.width=widths, edge.lty=styles, edge.color=\"black\", edge.arrow.size=1.5)\n\n";

	rcode << "#finish png file\n";
	rcode << "dev.off()\n\n";

	rcode.close();
};

//! Initialises the task for outputting the network priors.
void OutputPriorsTask::initialiseTask()
{
	
		
};

//! Outputs the task header for outputting the network priors.
void OutputPriorsTask::outputTaskHeader()
{
	outputTaskName();
	out("Outputting priors\n");
};

//! Outputs the task details for outputting the network priors.
void OutputPriorsTask::outputTaskDetails()
{
	
	out("Network: "); out(networkName); out("\n");
	out("Network Structure: "); out(network->getNetworkString()); out("\n");
	out("Output priors to file: ");
	out(filename); out("\n");
		
};
	
//! Outputs the network priors.
void OutputPriorsTask::doTask()
{
	//choose network to perform analysis
	if(networkName != "") network =	taskCentral->getNetwork(networkName);
	else network = taskCentral->getLatestNetwork(networkName);

	//check all of the Nodes have the same amount of data
	taskCentral->checkData(); //do not update missing data here
	taskCentral->addNetworkMissingData(network->updateNetworkMissingData());

	//output in priors
	network->calculatePriorsAndPosteriors();
	network->outputPriors(filename);
	
};

//! Initialises the task for outputting the network posteriors.
void OutputPosteriorsTask::initialiseTask()
{
	
		
};

//! Outputs the task header for outputting the network posteriors.
void OutputPosteriorsTask::outputTaskHeader()
{
	outputTaskName();
	out("Outputting posteriors\n");
};

//! Outputs the task details for outputting the network posteriors.
void OutputPosteriorsTask::outputTaskDetails()
{
	out("Network: "); out(networkName); out("\n");
	out("Network Structure: "); out(network->getNetworkString()); out("\n");
	out("Output posteriors to file: ");
	out(filename); out("\n");
		
};
	
//! Outputs the network posteriors.
void OutputPosteriorsTask::doTask()
{
	//choose network to perform analysis
	if(networkName != "") network =	taskCentral->getNetwork(networkName);
	else network = taskCentral->getLatestNetwork(networkName);

	//check all of the Nodes have the same amount of data
	taskCentral->checkData(); //do not update missing data here
	taskCentral->addNetworkMissingData(network->updateNetworkMissingData());

	//output in posteriors
	network->calculatePriorsAndPosteriors();
	network->calcScore(); //need to do this for discrete posteriors
	network->outputPosteriors(filename);
};

//! Initialises the task for imputing network data.
void ImputeNetworkDataTask::initialiseTask()
{
	if(jobNo > 0) taskCentral->setStartAndEndIndivsBasedOnJobNo(jobNo, jobTotal, impStartIndivNo, impEndIndivNo);		
};

//! Outputs the task header for for imputing network data.
void ImputeNetworkDataTask::outputTaskHeader()
{
	outputTaskName();
	out("Imputing network data\n");
};

//! Outputs the task details for imputing network data.
void ImputeNetworkDataTask::outputTaskDetails()
{
	out("Network: "); out(networkName); out("\n");
	out("Network Structure: "); out(network->getNetworkString()); out("\n");
	if(jobNo > 0) {out("Job "); out(jobNo); out(" of "); out(jobTotal); out("\n");};
	if(impStartIndivNo > 0) {out("Imputing individuals "); out(impStartIndivNo); out(" to "); out(impEndIndivNo); out("\n");};
	out("Number of individuals with missing data: "); out(noIndivsWithMissing); out("\n");
	out("Number of individuals imputed: "); out(noIndivsImputed); out("\n");
	double percentImputed1 = 100 * (1 - ((double)noBitsOfDataNotImputed1)/((double)noBitsOfDataToImpute)); //always 100% now
	double percentImputed2 = 100 * (1 - ((double)noBitsOfDataNotImputed2)/((double)noBitsOfDataToImpute)) - percentImputed1; 
	double percentImputed3 = 100 * (1 - ((double)noBitsOfDataNotImputed3)/((double)noBitsOfDataToImpute)) - percentImputed1 - percentImputed2;
	double percent4 = 100 - percentImputed1 - percentImputed2 - percentImputed3;
	out("Percentage of data imputed (at each step): "); 
	if(noIndivsWithMissing == 0) out("N/A"); else
	{
		out(percentImputed1);
		if(percentImputed2 > 0 || percentImputed3 > 0 || percent4 > 0) {out("..."); out(percentImputed2);};
		if(percentImputed3 > 0 || percent4 > 0) {out("..."); out(percentImputed3);};
		if(percent4 > 0) {out("..."); out(percent4);};		
	};
	out("\n");
	if(useRandomData) out("Training data: missing data replaced with randomly sampled values\n");
	else out("Training data: individuals with complete data\n");
	out("Minimum percentage of non-missing edges (or singleton nodes) required to impute individual: "); out(100*minNonMissingEdges); out("\n");
	out("Random restarts: "); out(randomRestarts); out("\n");
	out("Random jitter restarts: "); out(jitterRestarts); out("\n");
	
	out("\n");
};
	
//! Starts the task of imputing the network data.
void ImputeNetworkDataTask::doTask()
{
	//choose network to perform analysis
	if(networkName != "") network =	taskCentral->getNetwork(networkName);
	else network = taskCentral->getLatestNetwork(networkName);

	//check all of the Nodes have the same amount of data
	taskCentral->checkData(); //do not update missing data here

	taskCentral->addNetworkMissingData(network->updateNetworkMissingData());

	//impute the network data
	doNetworkDataImputation();

	noBitsOfDataNotImputed1 = noBitsOfDataNotImputed;
	unsigned int noBitsOfDataToImpute1 = noBitsOfDataToImpute;
	noBitsOfDataNotImputed2 = 0;
	noBitsOfDataNotImputed3 = 0;

	if(noBitsOfDataNotImputed1 > 0)
	{
		//try again
		doNetworkDataImputation(false, true, 2);
		noBitsOfDataNotImputed2 = noBitsOfDataNotImputed;
		if(noBitsOfDataNotImputed2 > 0)
		{
			//and again without adjusting
			doNetworkDataImputation(true, true, 3);
			noBitsOfDataNotImputed3 = noBitsOfDataNotImputed;
			if(noBitsOfDataNotImputed3 > 0)
			{
				//resort to random replacement
				network->randomlyFillMissingValues();
			};
		};
	};

	noBitsOfDataToImpute = noBitsOfDataToImpute1;

	taskCentral->removeNetwork(name);
	bootstrapNetwork->deleteAllNetworkNodeData();
	delete bootstrapNetwork;
};

//! Imputes the network data.
void ImputeNetworkDataTask::doNetworkDataImputation(const bool & doNotDoAdjust, const bool & updateImpImmed, const unsigned int & tryNo)
{	
	//loop thro' data pts, i.e 1 to N individuals
	unsigned int amountData = taskCentral->getAmountOfData();
	list<unsigned int> nodesWithMissingData;
	//set<unsigned int> missingNodesGrouped;
	map<unsigned int, set<unsigned int> > groupNodesWithMissingData; //groupId, nodeID
	map<unsigned int, set<unsigned int> > groupNodesWithCompleteData; //groupID as above, nodeID
	map<unsigned int, set<unsigned int> > completeDataPartialNodes; //node ID, node IDs to 'take away'
	map<unsigned int, set<unsigned int> >::const_iterator gwmd2;
	map<unsigned int, set<unsigned int> >::iterator gwmd3;
	set<unsigned int> groupsToDelete;
	unsigned int count, ranGroup, groupToAddTo;
	set<unsigned int>::const_iterator md;
	set<unsigned int> possibleGroupsToAddNodeTo;
	set<unsigned int>::const_iterator pgtat;
	set<unsigned int> missingDataNodesSet;
	set<unsigned int>::const_iterator mdns, mdns2;
	stack<unsigned int> missingDataNodes;
	stack<unsigned int> connectedMissingDataNodes;
	unsigned int missNodeID;
	map<unsigned int, Node *> bootNetworkNodes;
	//unsigned int subsetPercent = 90; //default 90, but changable

	list<unsigned int> connectedBootNodes;
	set<unsigned int>::const_iterator nd;
	string prevNetworkStr = network->getNetworkString(); // used as starting pt for bootstrap data search if req'd

	//decide if using randomly sampled data to replace missing values or complete individuals as training data
	if(!useCompleteData && !useRandomData)
	{
		pair<unsigned int, unsigned int> missingNotMissing = network->getNoMissingNotMissing();	
		if(missingNotMissing.first < minIndivsForCompleteData) useRandomData = true; else useCompleteData = true;
	};

	//record if data was imputed then update after all indivs are imputed
	network->setupWasImputed();

	//record that data was imputed
	network->setIsImputedData();

	//setup iterators for data
	network->setInitialDataIterators2();

	//setup standard devs for cts nodes
	network->setStandardDevs();

	//setup bootstrap network
	BootstrapNetworkTaskHelp bootstrapNetworkTaskHelp;

	//get copy of the initial network
	string bootNetworkName = name; 
	string preNodeName = name + "-";
	
	//setup so if subsetPercent == 0 then the "bootstrap" network is only fitted once to speed things up
	unsigned int subsetPercent2 = subsetPercent;
	bool bootstrapFitted = false;
	if(subsetPercent == 0) subsetPercent2 = 100;
	
	if(tryNo==1)
	{
		bootstrapNetwork = bootstrapNetworkTaskHelp.createBootstrapNetwork(taskCentral, network, bootNetworkName);

		//set up list of data for bootstrapping (ok, we are not really bootstraping, we are taking a randomly sampled subset here)
		bootstrapNetworkTaskHelp.setupBootstrap(taskCentral, bootstrapNetwork, network, useRandomData, subsetPercent2);
	};

	//setup search object to do the searches for the bootstrapping
	SearchNetworkModelsTask searchNetworkModelsTask;
	searchNetworkModelsTask.setNetworkName(bootNetworkName);
	searchNetworkModelsTask.setRandomRestarts(randomRestarts);
	searchNetworkModelsTask.setJitterRestarts(jitterRestarts);
	searchNetworkModelsTask.setTaskCentral(taskCentral); 
	searchNetworkModelsTask.initialiseTask();
	searchNetworkModelsTask.setPreNodeName(preNodeName);
	searchNetworkModelsTask.setCopyCoeffs(true);

	//set up node maps between bootstrap network and the original network
	network->setNodeMap(preNodeName, bootstrapNetwork);

	if(tryNo==1)
	{
		noIndivsWithMissing = 0;
		noIndivsImputed = 0;
	};

	noBitsOfDataToImpute = 0;
	noBitsOfDataNotImputed = 0;
//	bool nodeAdded;
	unsigned int groupID = 1;
//	unsigned int addedGroupID;
	set<unsigned int> someNodes;
	//unsigned int nodeIDBootstrap;
	//unsigned int nodeIDOrig;
	//unsigned int OLD = 0;
	//bool noJoin = false;//true; //if missing data is joined in network, do not impute or if A->B<-C, A and C missing
	bool skipImp = false;
	unsigned int noEdgesNotMissing;
	unsigned int totalEdges;
	unsigned int noNonMissingSingltonNodes;
	unsigned int totalSingltonNodes;
	Node * aNode;
	list<unsigned int> otherFactors;
	
	//bootstrapNetwork->setNetwork("[sA1][sA2][sB1][sB2][sC1][sCD][sD2][se1][se2][sf1][sf2][sg1][sg2][sh1][sh2][sh3][sN1][sN2][sM1][sM2][exA|sA1:sA2][exB|sB1:sB2][exC|sC1:sCD][exD|sCD:sD2][exN|sN1:sN2][exM|sM1:sM2][trait|exA:exB:exC:exD][exe|trait:se1:se2][exf|exe:sf1:sf2][exg|trait:sg1:sg2][exh|trait:sh1:sh2:sh3]", preNodeName);

	for(unsigned int indiv = 1; indiv <= amountData; ++indiv)
	{
		if((impStartIndivNo == 0 || indiv >= impStartIndivNo) && (impEndIndivNo == 0 || indiv <= impEndIndivNo) && (impStartIndivNo <= impEndIndivNo))
		{
			nodesWithMissingData = network->getNodesWithMissingData(); //for data iterator 2, corresponding to indiv no

			
			if(nodesWithMissingData.size() > 0)
			{
				if(tryNo==1) noIndivsWithMissing++;
				
				//do imputation
				if(maxMissing == 0 || nodesWithMissingData.size() <= maxMissing)
				{
					if(freeClearMemory)
					{
						map<unsigned int, set<unsigned int> >().swap(groupNodesWithMissingData); 
						map<unsigned int, set<unsigned int> >().swap(groupNodesWithCompleteData);
					}
					else
					{
						groupNodesWithMissingData.clear();
						groupNodesWithCompleteData.clear();
						//missingNodesGrouped.clear();
					};

					if(subsetPercent != 0 || !bootstrapFitted)
					{
						//bootstrap the complete data
				
						bootstrapNetworkTaskHelp.updateBootstrapDataSubset(bootstrapNetwork, network, useRandomData, subsetPercent2);
				
						bootstrapNetwork->removeAllEdges(); // with data changed the current
						if(usePrevNetwork) bootstrapNetwork->setNetwork(prevNetworkStr, preNodeName);
					
						bootstrapNetwork->updateBlackWhiteDifferentData(network, preNodeName);
					
						//fit network with bootstrap data
						searchNetworkModelsTask.doTask();
						bootstrapFitted = true;
					};

	
					//set groups of missing node data
					for(list<unsigned int>::const_iterator nwmd = nodesWithMissingData.begin(); nwmd != nodesWithMissingData.end(); ++nwmd)
					{
						//missingDataNodes.push(*nwmd);
						missingDataNodesSet.insert(*nwmd);
					};
					
					groupID = 0;

					while(!missingDataNodesSet.empty())
					{
						missNodeID = *missingDataNodesSet.begin();						
						if(freeClearMemory) set<unsigned int>().swap(someNodes); else someNodes.clear();
						someNodes.insert(missNodeID);
						missingDataNodesSet.erase(missNodeID);
						
						//put factor nodes in the same group to be imputed
						if(bootstrapNetwork->getNetworkNode(bootstrapNetwork->convertID(missNodeID))->getIsFactorNode())
						{
							//add other factors variables if this node is a head factor node
							otherFactors = taskCentral->getAllNodeData()->getFactorDataGroup(missNodeID);
						
							for(list<unsigned int>::const_iterator of = otherFactors.begin(); of != otherFactors.end(); ++of)
							{
								mdns2 = missingDataNodesSet.find(*of);
								if(mdns2 != missingDataNodesSet.end())
								{
									someNodes.insert(*of);
									missingDataNodesSet.erase(*of);
								};
							};
						};					

						groupID++;
						groupNodesWithMissingData[groupID] = someNodes;
					};
				

					//set up connected nodes with complete data for each group of missing data
					for(map<unsigned int, set<unsigned int> >::const_iterator gwmd = groupNodesWithMissingData.begin(); gwmd != groupNodesWithMissingData.end(); ++gwmd)
					{
						if(freeClearMemory) set<unsigned int>().swap(someNodes); else someNodes.clear();
						for(set<unsigned int>::const_iterator md = gwmd->second.begin(); md != gwmd->second.end(); ++md)
						{
							if(!useAllNN)
							{
								connectedBootNodes = bootstrapNetwork->getConnectedNodes(bootstrapNetwork->convertID(*md));

								for(list<unsigned int>::iterator cn3 = connectedBootNodes.begin(); cn3 != connectedBootNodes.end(); ++cn3)
								{
									//check the node has complete data for the indiv that is being imputed
									aNode = network->getNetworkNode(network->convertID(*cn3));
									if(!aNode->nodeDataIsMissing2()) someNodes.insert(aNode->getNodeID());
								};
							}
							else
							{
								//use all other nodes with complete data - just here for evaluation purposes not to be used for proper analysis of real data
								bootNetworkNodes = bootstrapNetwork->getNetworkNodes();

								for(map<unsigned int, Node *>::iterator bnn = bootNetworkNodes.begin(); bnn != bootNetworkNodes.end(); ++bnn)
								{
									//check the node has complete data for the indiv that is being imputed
									aNode = network->getNetworkNode(network->convertID(bnn->first));
									if(!aNode->nodeDataIsMissing2()) someNodes.insert(aNode->getNodeID());
								};
							};
						};

						groupNodesWithCompleteData[gwmd->first] = someNodes;
					};
				
					//find single missing nodes which have no connected complete data nodes, then if connected to other missing node(s) then add to this one of these groups randomly
					if(false)
					{

					if(freeClearMemory) set<unsigned int>().swap(groupsToDelete); else groupsToDelete.clear();
					gwmd2 = groupNodesWithMissingData.begin();
					for(map<unsigned int, set<unsigned int> >::const_iterator gwcd2 = groupNodesWithCompleteData.begin(); gwcd2 != groupNodesWithCompleteData.end(); ++gwcd2, ++gwmd2)			
					{
						if(gwcd2->second.size() == 0)
						{
							connectedBootNodes = bootstrapNetwork->getConnectedNodes(bootstrapNetwork->convertID(*gwmd2->second.begin()));
						
							if(freeClearMemory) set<unsigned int>().swap(possibleGroupsToAddNodeTo); else possibleGroupsToAddNodeTo.clear();

							for(list<unsigned int>::iterator cn3 = connectedBootNodes.begin(); cn3 != connectedBootNodes.end(); ++cn3)
							{
								//check the node has complete data for the indiv that is being imputed
								aNode = network->getNetworkNode(network->convertID(*cn3));
								if(aNode->nodeDataIsMissing2())
								{
									//find group this node is in
									gwmd3 = groupNodesWithMissingData.begin();
									while(gwmd3 != groupNodesWithMissingData.end())
									{
										md = gwmd3->second.find(aNode->getNodeID());
										if(md != gwmd3->second.end())
										{
											possibleGroupsToAddNodeTo.insert(gwmd3->first);
											groupsToDelete.insert(gwmd2->first);
										};

										++gwmd3;
									};
								};
							};

							//add to random connected group when the above conditions are met, then delete original group later
							if(possibleGroupsToAddNodeTo.size() > 0)
							{
								if(possibleGroupsToAddNodeTo.size() == 1)
								{
									groupToAddTo = *possibleGroupsToAddNodeTo.begin();															
								}
								else 
								{
									ranGroup = rand() % possibleGroupsToAddNodeTo.size() + 1;
									pgtat = possibleGroupsToAddNodeTo.begin();
									count = 1;
									while(pgtat != possibleGroupsToAddNodeTo.end() && count < ranGroup)
									{
										if(count == ranGroup)
										{
											groupToAddTo = *pgtat;
										};

										++pgtat;
										++count;
									};
								};

								gwmd3 = groupNodesWithMissingData.find(groupToAddTo);
								someNodes = gwmd3->second;
								someNodes.insert(*gwmd2->second.begin());
								groupNodesWithMissingData[*possibleGroupsToAddNodeTo.begin()] = someNodes;	
							};
						};
					};

					//now delete groups that were merged with other missing data that had connected nodes with data
					for(set<unsigned int>::const_iterator gtd = groupsToDelete.begin(); gtd != groupsToDelete.end(); ++gtd)
					{
						gwmd3 = groupNodesWithMissingData.find(*gtd);
						if(gwmd3 != groupNodesWithMissingData.end()) groupNodesWithMissingData.erase(gwmd3);
					};
					};//end of merging groups

					//try no joining... not great
					skipImp = false;
					/*if(noJoin)
					{
						for(map<unsigned int, set<unsigned int> >::iterator gwmdJ = groupNodesWithMissingData.begin(); gwmdJ != groupNodesWithMissingData.end(); ++gwmdJ)
						{
							if(gwmdJ->second.size() > 1) skipImp = true;
						};

					};*/
				
					//try require minimum no of edges not connected to a missing edge
					if(useMinNonMissingEdges)
					{
						bootstrapNetwork->getNoEdgesNotMissing(groupNodesWithMissingData, network, noEdgesNotMissing, totalEdges, noNonMissingSingltonNodes, totalSingltonNodes);
						if((double)((double)noEdgesNotMissing+(double)noNonMissingSingltonNodes)/(double)((double)totalEdges+(double)totalSingltonNodes) < minNonMissingEdges) skipImp = true; 			
					};
				
					//loop thro' data and find nearest neighbour
					if(!skipImp)
					{
						if(tryNo==1) noIndivsImputed++;
						if(useMean) noBitsOfDataNotImputed += network->imputeDataMean(groupNodesWithMissingData, bootstrapNetwork, updateImpImmed);
						else noBitsOfDataNotImputed += network->imputeDataNN(groupNodesWithMissingData, groupNodesWithCompleteData, bootstrapNetwork, doNotDoAdjust, updateImpImmed);
					};
				
					noBitsOfDataToImpute += nodesWithMissingData.size();
				};

			};//end skipping indivs
		};//end indiv with missing data

		network->advanceDataIterators2();
	};

	//actually set imputed data as non missing now, values have already been set
	network->updateImputedDataAsNonMissing();

	//remove the bootstrap network ... later
	//if(noBitsOfDataNotImputed==0) taskCentral->removeNetwork(bootNetworkName); 

	//remove bootstrap data nodes...
	//

	//update the list of missing data, now that data has been imputed for missing data
	taskCentral->addNetworkMissingData(network->updateNetworkMissingData());
};



//! Initialises the task for estimating the recall and precision for imputing network data.
void ImputeEstimateRecallPrecisionTask::initialiseTask()
{
			
};

//! Outputs the task header for estimating the recall and precision for imputing network data.
void ImputeEstimateRecallPrecisionTask::outputTaskHeader()
{
	outputTaskName();
	out("Estimating the recall and precision when imputing network data\n");
};

//! Outputs the task details estimating the recall and precision for imputing network data.
void ImputeEstimateRecallPrecisionTask::outputTaskDetails()
{
	out("Network: "); out(networkName); out("\n");
	//out("Network Structure: "); out(network->getNetworkString()); out("\n");
	out("Number of iterations: "); out(iterations); out("\n");
	out("Random restarts: "); out(randomRestarts); out("\n");
	out("Random jitter restarts: "); out(jitterRestarts); out("\n");
	out("Minimum percentage of non-missing edges (or singleton nodes) required to impute individual: "); out(100*minNonMissingEdges); out("\n");
	out("Individuals with data: "); out(missingNotMissing.second); out("\n");
	out("Individuals with missing data: "); out(missingNotMissing.first); out("\n");
	if(simDataFilename != "") {out("Simulation data output to file: "); out(simDataFilename); out("\n");};
	if(simNetFilename != "") {out("Simulation network output to file: "); out(simNetFilename); out("\n");};

	if(missingNotMissing.second > 1)
	{
		out("\n");
		out("Recall: the percentage of edges found from the original true network\n");
		out("Precision: the percentage of edges in the network that are also in the original true network\n");
		out("\n");
		out("                                   Recall     Precision\n");		
		out("No imputation                      "); out(toString2DP(recallNoImp*100)); out("     "); out(toString2DP(preNoImp*100)); out("\n");
		if(doImps)
		{
			out("Imputation                         "); out(toString2DP(recallImpRT*100)); out("     "); out(toString2DP(preImpRT*100)); out("\n");
			out("Imputation (complete training)     "); out(toString2DP(recallImp*100)); out("     "); out(toString2DP(preImp*100)); out("\n");
		};
		out("Full data (upper limit)            "); out(toString2DP(recallFull*100)); out("     "); out(toString2DP(preFull*100)); out("\n");
		out("\n");
	}
	else {
		out("\n");
		out("Not enough data to fit an initial network to base estimations from!\n");
	}
};
	
//! Starts the task of estimating the recall and precision for imputing network data.
void ImputeEstimateRecallPrecisionTask::doTask()
{
	//choose network to perform analysis
	if(networkName != "") network =	taskCentral->getNetwork(networkName);
	else network = taskCentral->getLatestNetwork(networkName);

	string simDataTaskName = "s d"; //sim data
	string simDataPreNodeName = simDataTaskName + "-";

	//check all of the Nodes have the same amount of data
	taskCentral->checkData(); //do not update missing data here
	NetworkMissingData * networkMissingData = network->updateNetworkMissingData();
	taskCentral->addNetworkMissingData(networkMissingData);
	
	missingNotMissing = network->getNoMissingNotMissing();
	if(missingNotMissing.second <= 1) return;

	//NEED to copy node level missing data for later then write fn to replicate the missing data pattern
	SearchNetworkModelsTask searchNetworkModelsTask;
	SimulateNetworkDataTask simulateNetworkDataTask;
	unsigned int noSims = taskCentral->getAmountOfData();
	
	ImputeNetworkDataTask imputeNetworkDataTask; 

	string dummyNetName = "s d dummy";
	string dummyPreNodeName = dummyNetName + "-";
	Network * dummyNet;

	//sim data just to replicate missing data pattern
	simulateNetworkDataTask.setNetworkName(networkName);
	simulateNetworkDataTask.setCreateDifferentNodeData(true);
	simulateNetworkDataTask.setTaskName(dummyNetName);
	simulateNetworkDataTask.setTaskCentral(taskCentral);
	simulateNetworkDataTask.setNoSims(noSims);
	simulateNetworkDataTask.initialiseTask();
	simulateNetworkDataTask.doTask();
	
	dummyNet = taskCentral->getNetwork(dummyNetName);

	//add missingness to dummy data as in orig data
	network->copyMissingness(dummyNet, dummyPreNodeName);

	for(unsigned int it = 1; it <= iterations; ++it)
	{
		
		///////////////////////////////////////////////////
		////fit the data to a network, call it network-orig
			
		//impute data to get initial estimate of network to use for further sims to test methods against
		imputeNetworkDataTask.setTaskName("i ");
		imputeNetworkDataTask.setNetworkName(networkName);
		imputeNetworkDataTask.setRandomRestarts(randomRestarts);
		imputeNetworkDataTask.setJitterRestarts(jitterRestarts);
		imputeNetworkDataTask.setTaskCentral(taskCentral); 
		imputeNetworkDataTask.setMinNonMissingEdges(minNonMissingEdges);
		imputeNetworkDataTask.initialiseTask();	
		imputeNetworkDataTask.setUseRandomData();

		imputeNetworkDataTask.doTask();

		//fit network to missing sim data, call it network-miss
		searchNetworkModelsTask.setNetworkName(networkName);
		searchNetworkModelsTask.setRandomRestarts(randomRestarts);
		searchNetworkModelsTask.setJitterRestarts(jitterRestarts);
		searchNetworkModelsTask.setTaskCentral(taskCentral); 
		searchNetworkModelsTask.initialiseTask();
		//searchNetworkModelsTask.setPreNodeName(preNodeName);
		searchNetworkModelsTask.setCopyCoeffs(true);

		searchNetworkModelsTask.initialiseTask();
		searchNetworkModelsTask.doTask();

		//add back missingness to orig data as in copied sim data
		dummyNet->copyMissingness(network, dummyPreNodeName, true);
	
		if(iterations == 1) dummyNet->deleteAllNetworkNodeData();
	
		/////////////////////////////////////////////////
		//simulate data for the fitted network	
		
		//simulateNetworkDataTask.set
		simulateNetworkDataTask.setNetworkName(networkName);
		simulateNetworkDataTask.setCreateDifferentNodeData(true);
		simulateNetworkDataTask.setTaskName(simDataTaskName);
		simulateNetworkDataTask.setTaskCentral(taskCentral);
		simulateNetworkDataTask.setNoSims(noSims);
		simulateNetworkDataTask.initialiseTask();
		simulateNetworkDataTask.doTask();
	
		simNetwork = taskCentral->getNetwork(simDataTaskName);
	
		simNetwork->removeAllEdges();

		/////////////////////////////////////////////////
		//fit a network to the full data, call it network-full
		searchNetworkModelsTask.setNetworkName(simDataTaskName);		
		searchNetworkModelsTask.doTask();
		
		/////////////////////////////////////////////////
		//compare the full net to the orig net to get recall and precision
	
		network->calcRecallPrecision(simNetwork, recallFull, simDataPreNodeName, true);
		simNetwork->calcRecallPrecision(network, preFull, "", false);
	
		/////////////////////////////////////////////////
		//add missingness to sim data as in orig data
		network->copyMissingness(simNetwork, simDataPreNodeName);
	
		simNetwork->removeAllEdges();
		/////////////////////////////////////////////////
		//fit network to missing sim data, call it network-miss
		searchNetworkModelsTask.setNetworkName(simDataTaskName);
		searchNetworkModelsTask.setRandomRestarts(randomRestarts);
		searchNetworkModelsTask.setJitterRestarts(jitterRestarts);
		searchNetworkModelsTask.setTaskCentral(taskCentral); 
		searchNetworkModelsTask.initialiseTask();	
		searchNetworkModelsTask.setCopyCoeffs(true);
	
		searchNetworkModelsTask.doTask();
		
		/////////////////////////////////////////////////
		//compare the miss net to the orig net to get recall and precision
		network->calcRecallPrecision(simNetwork, recallNoImp, simDataPreNodeName, true);
		simNetwork->calcRecallPrecision(network, preNoImp, "", false);

		if(doImps)
		{
			/////////////////////////////////////////////////
			//impute data using complete training, call it network-imp-com
		
			imputeNetworkDataTask.setTaskName("i ");
			imputeNetworkDataTask.setNetworkName(simDataTaskName);
			imputeNetworkDataTask.setRandomRestarts(randomRestarts);
			imputeNetworkDataTask.setJitterRestarts(jitterRestarts);
			imputeNetworkDataTask.setTaskCentral(taskCentral); 
			imputeNetworkDataTask.setMinNonMissingEdges(minNonMissingEdges);
			imputeNetworkDataTask.initialiseTask();	
			imputeNetworkDataTask.setUseCompleteData();
		
			imputeNetworkDataTask.doTask();
		
			simNetwork->removeAllEdges();
			//fit network to missing sim data, call it network-miss
			searchNetworkModelsTask.initialiseTask();
			searchNetworkModelsTask.doTask();
		
			/////////////////////////////////////////////////
			//compare the imp complete net to the orig net to get recall and precision
			network->calcRecallPrecision(simNetwork, recallImp, simDataPreNodeName, true);
			simNetwork->calcRecallPrecision(network, preImp, "", false);
		
			/////////////////////////////////////////////////
			//impute data using random training, call it network-imp-rt
	
			//add missingness to sim data as in orig data
			network->copyMissingness(simNetwork, simDataPreNodeName);
		
			imputeNetworkDataTask.setUseRandomData();	
			imputeNetworkDataTask.initialiseTask();	
			imputeNetworkDataTask.doTask();
	
			simNetwork->removeAllEdges();
			//fit network to missing sim data, call it network-miss	
			searchNetworkModelsTask.initialiseTask();	
			searchNetworkModelsTask.doTask();
		
			/////////////////////////////////////////////////
			//compare the imp RT net to the orig net to get recall and precision
			network->calcRecallPrecision(simNetwork, recallImpRT, simDataPreNodeName, true);
			simNetwork->calcRecallPrecision(network, preImpRT, "", false);
		};

		if(simNetFilename != "") outputSimNetwork(it);
		if(simDataFilename != "")
		{
			//add missingness to sim data as in orig data
			network->copyMissingness(simNetwork, simDataPreNodeName);
	
			outputSimData(it);
		};

		taskCentral->removeNetwork(simDataTaskName);
		simNetwork->deleteAllNetworkNodeData();
		delete simNetwork;
	};

	recallNoImp /= (double)iterations; preNoImp /= (double)iterations;
	recallImp /= (double)iterations; preImp /= (double)iterations;
	recallImpRT /= (double)iterations; preImpRT /= (double)iterations;
	recallFull /= (double)iterations; preFull /= (double)iterations;
};

//! Outputs sim data to file for use elsewhere. 
void ImputeEstimateRecallPrecisionTask::outputSimData(const unsigned int & it)
{
	bool outputBedFile = false;
	bool wasDisData = false;
	bool wasCtsData = false;
	bool wasSNPDisData = false;
	bool wasSNPCtsData = false;
	bool wasSNPData = false;
	unsigned int jobNo = 0;
	unsigned int jobTotal = 0;
	unsigned int startIndivNo = 0;
	unsigned int endIndivNo = 0;

	string filename = simDataFilename;
	if(iterations > 1) filename += "-" + toString(it);

	simNetwork->outputNodeData(filename, outputBedFile, wasDisData, wasCtsData, wasSNPDisData, wasSNPCtsData, wasSNPData, startIndivNo, endIndivNo, jobNo);
};

//! Outputs sim network.
void ImputeEstimateRecallPrecisionTask::outputSimNetwork(const unsigned int & it)
{
	string filename = simNetFilename;
	if(iterations > 1) filename += "-" + toString(it) + ".dat"; else filename += ".dat";

	ofstream outNet(filename.c_str());

	outNet << simNetwork->getNetworkString(0) << "\n";
	
	outNet.close();
};

void CalculateRecallPrecisionTask::initialiseTask()
{


};

void CalculateRecallPrecisionTask::outputTaskHeader()
{
	outputTaskName();
	out("Calculating the recall and precision\n");
};

void CalculateRecallPrecisionTask::outputTaskDetails()
{
	out("Network: "); out(networkName); out("\n");
	out("Network Structure: "); out(network->getNetworkString()); out("\n");
	out("True Network: "); out(trueNetworkName); out("\n");
	out("True Network Structure: "); out(trueNetwork->getNetworkString()); out("\n");
	if(filename != "") { out("Recall and precision written to file: "); out(filename); out("\n"); };
	out("\n");
	out("Recall: the percentage of edges found from the original true network\n");
	out("Precision: the percentage of edges in the network that are also in the original true network\n");
	out("\n");
	out("Recall: "); out(toString2DP(recall*100)); out("\n");
	out("Precision: "); out(toString2DP(precision*100)); out("\n");
	if(filename != "")
	{
		out("\n");
		out("Recall and precision written to file: "); out(filename); out("\n");
	};
};


void CalculateRecallPrecisionTask::doTask()
{
	//choose network to perform analysis
	if(networkName != "") network =	taskCentral->getNetwork(networkName);
	else network = taskCentral->getLatestNetwork(networkName);

	if(trueNetworkName == "")
	{
		out("You must set the true network structure to calculate the recall and precision!\n");
		return;
	};

	trueNetwork = taskCentral->getNetwork(trueNetworkName);

	//compare the net to the true net to get recall and precision
	trueNetwork->calcRecallPrecision(network, recall, "", true);
	network->calcRecallPrecision(trueNetwork, precision, "", false);

	//network->calcRecallPrecision(simNetwork, recallImpRT, simDataPreNodeName, true);
	//simNetwork->calcRecallPrecision(network, preImpRT, "", false);

	//write results to file if req'd
	if(filename != "")
	{
		ofstream outResults(filename.c_str());
		outResults << recall * 100 << " " << precision * 100 << "\n";
		outResults.close();
	};
};


//! Initialises the task for analysis of the robustness.
void MeasurementErrorRobustnessTask::initialiseTask()
{


};

//! Outputs task header for analysis of the robustness.
void MeasurementErrorRobustnessTask::outputTaskHeader()
{
	outputTaskName();
	out("Analysing robustness of the fitted network when accounting for measurement error\n");

};

//!  Outputs task details for analysis of the robustness.
void MeasurementErrorRobustnessTask::outputTaskDetails()
{
	out("Network: "); out(networkName); out("\n");
	
	if(measurementErrorFile != "")
	{
		out("Standard deviation of measurement errors given in file: "); out(measurementErrorFile); out("\n");		
	}
	
	if(measurementErrorMultipleFile != "")
	{
		out("Multiple of standard deviations defining measurement errors given in file: "); out(measurementErrorFile); out("\n");
	}

	if(measurementErrors.empty() && measurementErrorMultiples.empty())
	{
		if(measurementErrorMultiple != 0) { out("Multiple of standard deviations defining measurement error for all variables: "); out(measurementError); out("\n"); }
		else { out("Standard deviation of measurement error for all variables: "); out(measurementError); out("\n"); }
	}
	else
	{
		if(!measurementErrorMultiples.empty())
		{
			out("Multiple of standard deviations defining measurement error for each variable is set to:\n");
			for(map<unsigned int, double>::const_iterator me = measurementErrorMultiples.begin(); me != measurementErrorMultiples.end(); ++me)
			{
				out(network->getNetworkNode(me->first)->getDisplayName()); out(" ");  out(me->second); out("; ");
			};			
		}
		else if(!measurementErrors.empty())
		{
			out("Standard deviation of measurement error for each variable is set to:\n");
			for(map<unsigned int, double>::const_iterator me = measurementErrors.begin(); me != measurementErrors.end(); ++me)
			{
				out(network->getNetworkNode(me->first)->getDisplayName()); out(" ");  out(me->second); out("; ");
			};		
		};

		if(measurementErrorMultiple != 0) { out("\nFor other variables the multiple of standard deviations defining measurement error is set to: "); out(measurementErrorMultiple); out("\n"); }
		else { out("\nFor other variables the standard deviation of measurement error is set to: "); out(measurementError); out("\n"); }
	};
	
	if(useBootstraps) { out("Number of bootstrap iterations: "); out(noIterations); out("\n"); }
	else { out("Number of iterations: "); out(noIterations); out("\n"); };
	out("Random restarts: "); out(randomRestarts); out("\n");
	out("Random jitter restarts: "); out(jitterRestarts); out("\n");
	

	if(filename != "")
	{
		out("Average network output to file: "); out(filename); out("-ave.dat"); out("\n");
		out("Nodes summary output to file: "); out(filename); out("-nodes.dat"); out("\n");
		out("Edges summary output to file: "); out(filename); out("-edges.dat"); out("\n");
	};	
	if(igraphPrefix != "") { out("R code to plot average network: "); out(igraphPrefix + ".R"); out("\n"); };

	if(arcThresholdSet) { out("Set edge threshold: "); out(arcThreshold); out("\n"); }
	else { out("Estimated edge threshold: "); out(arcThreshold); out("\n"); };
	if(thresholdFilename != "") { out("Estimated edge threshold output to file: "); out(thresholdFilename); out("\n"); };
	out("Network structure (after above threshold): "); out(network->getNetworkString()); out("\n");
	out("Network score type: "); out(network->getScoreTypeName()); out("\n");
	if(network->getScoreFixName() != "none") { out("Network score fix: "); out(network->getScoreFixName()); out("\n"); };
	if(network->hasLoop()) { out("The network contains a loop!\n"); }
	if(!network->checkNetworkIsValid()) { out("The network is not valid!\n"); }
	else { out("Network score = "); out(network->calcScore()); out("\n"); };

};

//! Do task for analysing robustness of the fitted network when accounting for measurement error
void MeasurementErrorRobustnessTask::doTask()
{
	//choose network for analyses
	if(networkName != "") network = taskCentral->getNetwork(networkName);
	else network = taskCentral->getLatestNetwork(networkName);

	//set up average network task
	averageNetworksTask.setTaskCentral(taskCentral);
	averageNetworksTask.setNetworkName(networkName);	
	averageNetworksTask.setFilename(filename);
	averageNetworksTask.setThresholdFilename(thresholdFilename);
	averageNetworksTask.setFamilyFile(famFilename);
	if(arcThresholdSet) averageNetworksTask.setArcThreshold(arcThreshold);
	averageNetworksTask.setNoBootstraps(noIterations);
	averageNetworksTask.setUseBootstraps(useBootstraps);
	averageNetworksTask.setRandomRestarts(randomRestarts);
	averageNetworksTask.setJitterRestarts(jitterRestarts);
	averageNetworksTask.setRFilenamePrefix(igraphPrefix);
	averageNetworksTask.setUseEquivNets(useEquivNets);

	unsigned int nodeID;
	for(map<string, double>::const_iterator mse = measurementStrErrors.begin(); mse != measurementStrErrors.end(); ++mse)
	{
		nodeID = taskCentral->getNodeDataNumberInit(mse->first);
		if(nodeID != 0) measurementErrors[nodeID] = mse->second;
		else
		{
			string mess = "Attempt to add measurement error " + toString2DP(mse->second) + " to node \"" + mse->first + "\" but no such node exists!\n";
			exitErr(mess);
		};
	};

	for(map<string, double>::const_iterator msem = measurementStrErrorMultiples.begin(); msem != measurementStrErrorMultiples.end(); ++msem)
	{
		nodeID = taskCentral->getNodeDataNumberInit(msem->first);
		if(nodeID != 0) measurementErrorMultiples[nodeID] = msem->second;
		else
		{
			string mess = "Attempt to add measurement error multiple " + toString2DP(msem->second) + " to node \"" + msem->first + "\" but no such node exists!\n";
			exitErr(mess);
		};
	};

	averageNetworksTask.setMeasurementError(measurementError);	
	averageNetworksTask.setMeasurementErrors(measurementErrors);
	averageNetworksTask.setMeasurementErrorFile(measurementErrorFile);
	
	averageNetworksTask.setMeasurementErrorMultiple(measurementErrorMultiple);
	averageNetworksTask.setMeasurementErrorMultiples(measurementErrorMultiples);
	averageNetworksTask.setMeasurementErrorMultipleFile(measurementErrorMultipleFile);

	averageNetworksTask.setUseEffectSize(useEffectSize);

	//averageNetworksTask.setUseMeasurementError(true);

	averageNetworksTask.doTask();

	if(!arcThresholdSet) arcThreshold = averageNetworksTask.getArcThreshold();
};

void MeasurementErrorRobustnessTask::setMeasurementError(const string & nd, const double & me)
{
	measurementStrErrors[nd] = me;
};

void MeasurementErrorRobustnessTask::setMeasurementErrorMultiple(const string & nd, const double & me)
{
	measurementStrErrorMultiples[nd] = me;
};


