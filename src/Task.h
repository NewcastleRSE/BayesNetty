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


/*! \file Task.h
    \brief This file organises the tasks.
    
*/

#ifndef __TASK
#define __TASK

#include <list>
#include <string>
#include <map>
#include <set>
#include "main.h"
#include "Network.h"
#include "Search.h"

using namespace std; // initiates the "std" or "standard" namespace

void removeReturns(string & s);
list<string> splitString(const string & strParents, string delimiter = ":");
void trim(string & str, const string & whitespace = " \t");

//! Class to act as central storage and linking place for all Tasks.
class TaskCentral
{
private:

	AllNodeData * allNodeData; //potential data(nodes) for any network
	unsigned int fileNoCount;
	map<string, unsigned int> dataNameFileNo; //name of data input task and corresponding file no

	map<string, Network *> allNetworks;
	list<string> orderedNetworks;
	list<NetworkMissingData *> allNetworkMissingData; //keep separate as may be temporary networks for searches with nodes ref'ing these also
	unsigned int defaultScoreType; //0 = loglik, 1 = BIC, 2 = AIC, 3 = bayes with iss - not all available depending on network type

	unsigned int randomSeed;

public:

	TaskCentral() : randomSeed(0)
	{ 
		allNodeData = new AllNodeData();
		fileNoCount = 1;	
		defaultScoreType = 1;
	};

	//! Delete Task things
	~TaskCentral()
	{
		delete allNodeData;
		for(map<string, Network *>::iterator nt = allNetworks.begin(); nt != allNetworks.end(); ++nt) delete nt->second;
		allNetworks.clear();
		for(list<NetworkMissingData *>::iterator mnt = allNetworkMissingData.begin(); mnt != allNetworkMissingData.end(); ++mnt) delete *mnt;
		allNetworkMissingData.clear();
	};

	void addNetwork(const string & name, Network * network);
	void removeNetwork(const string & name);
	Network * getNetwork(const string & name);
	Network * getLatestNetwork(string & networkName);
	Network * getNoEdgeNetwork();
	void setRandomSeed(const unsigned int & seed) {randomSeed = seed;};
	unsigned int getRandomSeed() {return randomSeed;};
	void setDefaultScoreType(const unsigned int & dst) {defaultScoreType = dst;};
	void setupRestartNetwork(Network * restartNetwork, Network * network, const string & preNodeName = "");
	void updateRestartNetwork(Network * restartNetwork, Network * network, const unsigned int & minEdges, const unsigned int & maxEdges, const bool & doJitter);
	void copyNetworkEdges(Network * copyNetwork, Network * network);
	Network * createSimNetworkFromNetwork(Network * network, const string & taskName, const bool & createDifferentNodeData);
	void addNetworkMissingData(NetworkMissingData * nmd) {allNetworkMissingData.push_back(nmd);};
	AllNodeData * getAllNodeData() const {return allNodeData;};
	void addNodeData(Data * data) {allNodeData->addNodeData(data);};
	unsigned int getFileNo() {return fileNoCount;};
	void nextFileNo() {fileNoCount++;};
	void setFileNoFromDataName(const string & fileDataName, const unsigned int & fileNo) {dataNameFileNo[fileDataName] = fileNo;};
	unsigned int getFileNoFromDataName(const string & fileDataName);
	void checkData() {allNodeData->checkData();};
	void updateDataNames() {allNodeData->updateDataNames();};
	unsigned int getAmountOfData() const {return allNodeData->getAmountOfData();};
	void checkIDs(string & filename, list<string> & dataIDsToCheck) {allNodeData->checkIDs(filename, dataIDsToCheck);};
	void setNameOfIDs(list<string> & nids) {allNodeData->setNameOfIDs(nids);};
	list<string> getNameOfIDs() {return allNodeData->getNameOfIDs();};
	list<string> getDataIDs() {return allNodeData->getDataIDs();};
	unsigned int getNodeDataNumberInit(const string & nodeName) const {return allNodeData->getNodeDataNumberInit(nodeName);};
	void setStartAndEndIndivsBasedOnJobNo(unsigned int & jobNo, unsigned int & jobTotal, unsigned int & startIndivNo, unsigned int & endIndivNo);
};

//! Organises the correct analysis to perform.
class Task
{
private:
	

protected:
	string name;	
	TaskCentral * taskCentral;

public:

	Task() : name(""), taskCentral(0)
	{ 
		
	};

	//! Delete Task things
	virtual ~Task()
	{
		
	};

	virtual void initialiseTask() {};
	virtual void outputTaskHeader() {};
	virtual void outputTaskDetails() {};
	virtual void doTask() {};
	void setTaskName(const string & nm) {name= nm;};
	string getTaskName() {return name;};
	void outputTaskName() {out("Task name: "); out(name); out("\n");};
	void setTaskCentral(TaskCentral * tc) {taskCentral = tc;};
};

//! Calculate the Markov Blanket of a node.
class CalculateMarkovBlanketTask : public Task
{
private:

	string networkName;
	string nodeName;
	Network * network;
	set<Node *> markovBlanket;
	Network * blanketNetwork;
	
public:

	CalculateMarkovBlanketTask() : nodeName(""), network(0), blanketNetwork(0) {};

	~CalculateMarkovBlanketTask() {};

	//general task methods
	void initialiseTask();
	void outputTaskHeader();
	void outputTaskDetails();
	void doTask();

	//blanket methods
	void setNetworkName(const string & nm) {networkName = nm;};
	void setNodeName(const string & nm) {nodeName = nm;};
	
};

//! Stores info about a SNP.
struct SNPInfo
{
	unsigned int chromosome;
	string name;
	unsigned int bp;

	SNPInfo(unsigned int c, string n, unsigned int b) : chromosome(c), name(n), bp(b) {};

	~SNPInfo() {};
};

//! Names of SNPs.
class DescriptionOfSNPs
{
private:
	map<unsigned int, SNPInfo *> snpInfos; //snpNo, SNPInfo, base-pair positions of all SNPs in order with SNP names
	map<unsigned int, SNPInfo *>::const_iterator snpInfo;

public:

	DescriptionOfSNPs()
	{
		
	};

	~DescriptionOfSNPs()
	{
		for(map<unsigned int, SNPInfo *>::iterator si = snpInfos.begin(); si != snpInfos.end(); ++si) delete si->second;
	};

	//methods regarding the SNP descriptions
	void setUpSNPDesciptionData(string & filename);
	void advanceIterator() {++snpInfo;};
	SNPInfo * getSNPInfo() {return snpInfo->second;};
	void outputSNPInfo();

	unsigned int getNoSNPs() const
	{
		return snpInfos.size();
	};	
};

//! Input data, stores a collection of data, discrete cts and SNP data.
class InputDataTask : public Task
{
private:

	string filename, variableFilename;
	unsigned int fileno;
	set<string> includeVariables;
	set<string> excludeVariables;

	unsigned int taskDataType;
	//1 = cts data in a table, 2 = discrete data in a table, 3 = SNP data in .bed file
	//4 = SNP data in .bed to be treated as cts
	//5 = discrete data in a table, treated as a factor
	//6 = SNP data in .bed file, treated as a factor
	//7 = cts SNP data in a table, 8 = discrete SNP data in a table, 9 = factor SNP data in a table
	bool isCSV;
	char comma;		
	string line;
	string aNumber;
	string prefixName;
	istringstream lineStream;

	unsigned int totalNodeDatas;
	string ctsMissingValueStr;
	double ctsMissingValue;
	bool ctsMissingValueSet;
	string discreteMissingValue;
	unsigned int numberOfIDs;

	//for loading .bed files
	DescriptionOfSNPs descSNPs;

	unsigned int totalNoSubjects;
	unsigned int bitCount;
	int one;
	int aBit;	
	char buffer[1];
	ifstream readSNPData;

public:

	InputDataTask() : filename(""), variableFilename(""), fileno(0), includeVariables(), excludeVariables(), taskDataType(0), isCSV(false), comma(','), line(""), aNumber(""), prefixName(""), lineStream(""), descSNPs(), totalNodeDatas(0),
		ctsMissingValueStr(""), ctsMissingValue(-9), ctsMissingValueSet(false), discreteMissingValue("NA"), numberOfIDs(2), totalNoSubjects(0), bitCount(9), one('\1'), aBit(0) { buffer[0] = ' '; };

	~InputDataTask() {};

	//general task methods
	void initialiseTask();
	void outputTaskHeader();
	void outputTaskDetails();
	void doTask();

	//input task methods
	void addDataFile(const string & file) {filename = file;};
	void setDataType(const unsigned int & ft) {taskDataType = ft;};
	void setIsCSV(const bool & csv) {isCSV = csv;};
	void getStringInputData(ifstream & readData, string & aString);
	void getNumberInputData(ifstream & readData, double & aNum, bool & missing, string & badNumber);
	void setCtsMissingValue(const string & cmv);
	void setDiscreteMissingValue(const string & dmv) {discreteMissingValue = dmv;};
	void setIncludeVariables(const string & filename);
	void setExcludeVariables(const string & filename);
	void setPrefixName(const string& pfnm) {prefixName = pfnm;};
	void setNumberOfIDs(const unsigned int & no) {numberOfIDs = no;};
	void checkIncludeExcludeSet();
	bool isIncluded(const string & variable);
	void getVariableNames(list<string> & variableNames, const bool & secondLine = false);
	void loadSNPData();
	void loadSNPDataAsFactors();
	void loadCtsData(const bool & isSNPData = false);
	void loadDiscreteData(const bool & isSNPData = false);
	void loadDiscreteDataAsFactors(const bool & isSNPData = false);
	void handleFamFile();
	void handleBimFile();

	void openBinaryFileFirst();
	unsigned char getNextNoOfMinorAlleles();
};

//! Inputs network.
class InputNetworkTask : public Task
{
private:

	string filename;
	string filename2;
	string filenamePrefix; //for igraph
	string networkType;
	string whitelistFilename;
	string blacklistFilename;
	string localProbFilename;
	string score;
	string scoreFix;
	set<pair<string, string> > allowedEdgeTypes;
	set<pair<string, string> > notAllowedEdgeTypes;
	set<string> noParentsNodes;
	set<string> noChildrenNodes;
	map<pair<string, string>, double> costEdges;
	map<pair<string, string>, double> costEdgeTypes;

	Network * network;
	unsigned int totalDisNodes;
	unsigned int totalFactorNodes;
	unsigned int totalCtsNodes;
	unsigned int totalNoDataNodes;
	unsigned int totalEdges;
	double imaginarySampleSize;
	bool emptyNetwork;
	bool usingImputedData;
	
public:

	InputNetworkTask() : filename(""), filename2(""), filenamePrefix(""), networkType("bnlearn"), whitelistFilename(""), blacklistFilename(""), localProbFilename(""), score("BIC"), scoreFix(""), network(0), totalDisNodes(0), totalFactorNodes(0), totalCtsNodes(0), totalNoDataNodes(0), totalEdges(0),
		imaginarySampleSize(10), emptyNetwork(false), usingImputedData(false) {};

	~InputNetworkTask() {};

	//general task methods
	void initialiseTask();
	void outputTaskHeader();
	void outputTaskDetails();
	void doTask();
	
	//input network methods
	void setNetworkType(const string & nm);
	void addNetworkFile(const string & file) {filename = file;};
	void addNetworkFile2(const string & file) {filename2 = file;};
	void addNetworkFileIGraph(const string & file) {filenamePrefix = file;};
	void addWhitelistDataFile(const string & file) {whitelistFilename = file;};
	void addBlacklistDataFile(const string & file) {blacklistFilename = file;};
	void setLocalProbFilename(const string & file) {localProbFilename = file;};
	void setImaginarySampleSize(const double & i) {imaginarySampleSize = i;};
	void setEmptyNetwork(const bool & et) {emptyNetwork = et;};
	void setUsingImputedData() {usingImputedData = true;};
	void setScore(const string & sc) {score = sc;};
	void setScoreFix(const string & scf) {scoreFix = scf;};
	void loadNetworkData(const unsigned int & type); //type = 0, network; type = 1, whitelist; type = 2, blacklist
	void loadNetworkDataFormat1(const unsigned int & type); //BayesNetty format
	void loadNetworkDataFormat2(); //bnlearn format
	void loadNetworkDataFormat3(); //igraph format
	void setAllowedEdgeType(const string & dataName1, const string & dataName2) {allowedEdgeTypes.insert(make_pair(dataName1, dataName2));};
	void setNotAllowedEdgeType(const string & dataName1, const string & dataName2) {notAllowedEdgeTypes.insert(make_pair(dataName1, dataName2));};
	void setCostEdgeType(const string & dataName1, const string & dataName2, const double & cost) {costEdgeTypes[make_pair(dataName1, dataName2)] = cost;};
	void setCostEdge(const string & node1, const string & node2, const double & cost) {costEdges[make_pair(node1, node2)] = cost;};
	void addNoParentsNode(const string & npn) {noParentsNodes.insert(npn);};
	void addNoChildrenNode(const string & npn) {noChildrenNodes.insert(npn);};
	
};

//! Calculate the posterior of a network.
class CalculatePosteriorTask : public Task
{
private:

	string networkName;
	Network * network;
	bool hasLoop;

public:

	CalculatePosteriorTask() : networkName(""), network(0), hasLoop(false) {};

	~CalculatePosteriorTask() {};

	//general task methods
	void initialiseTask();
	void outputTaskHeader();
	void outputTaskDetails();
	void doTask();
	
	//calc posterior methods
	void setNetworkName(const string & nm) {networkName = nm;};
};

//! Calculate the score of a network.
class CalculateNetworkScoreTask : public Task
{
private:

	string networkName;
	Network * network;
	double networkScore;
	string bestNetworkStructure;
	bool hasLoop;
	string filename;
	string allScoresFilename;
	unsigned int noNetworksEval; //no of networks evaluated

public:

	CalculateNetworkScoreTask() : filename(""), allScoresFilename(""), noNetworksEval(0), network(0), networkScore(0), hasLoop(false) {};

	~CalculateNetworkScoreTask() {};

	//general task methods
	void initialiseTask();
	void outputTaskHeader();
	void outputTaskDetails();
	void doTask();
	
	//calc posterior methods
	void setNetworkName(const string & nm) {networkName = nm;};
	void setFilename(const string & fn) {filename = fn;};
	void setAllScoresFilename(const string & fn) {allScoresFilename = fn;};
	unsigned int calculateAllScores();
};

//! Search networks for best fitting model.
class SearchNetworkModelsTask : public Task
{
private:

	string networkName;
	SearchNetworks * searchNetworks;
	string filename;
	string searchType;
	double bestNetworkScore;
	Network * bestNetwork;
	Network * restartNetwork;
	unsigned int randomRestarts;
	unsigned int jitterRestarts;
	double minEdgeNodeRatio;
	double maxEdgeNodeRatio;
	unsigned int minEdges;
	unsigned int maxEdges;
	double minEdgeNodeRatioJitter;
	double maxEdgeNodeRatioJitter;
	unsigned int minEdgesJitter;
	unsigned int maxEdgesJitter;
	string preNodeName;
	bool copyCoeffs;

public:

	SearchNetworkModelsTask() : filename(""), searchType("Greedy"), randomRestarts(0), jitterRestarts(0), minEdgeNodeRatio(0.2), maxEdgeNodeRatio(3), minEdgeNodeRatioJitter(0.3), maxEdgeNodeRatioJitter(1), preNodeName(""), copyCoeffs(false),
	  bestNetwork(0), bestNetworkScore(0), maxEdges(0), maxEdgesJitter(0), minEdges(0), minEdgesJitter(0), restartNetwork(0), searchNetworks(0) {};

	~SearchNetworkModelsTask()
	{
		
	};

	//general task methods
	void initialiseTask();
	void outputTaskHeader();
	void outputTaskDetails();
	void doTask();
	
	//search network models methods
	void setNetworkName(const string & nm) {networkName = nm;};
	void setFilename(const string & fn) {filename = fn;};
	void setRandomRestarts(const unsigned int & rr) {randomRestarts = rr;};
	void setJitterRestarts(const unsigned int & jr) {jitterRestarts = jr;};
	void setPreNodeName(const string & pnm) {preNodeName = pnm;};
	void setCopyCoeffs(const bool & cpco) {copyCoeffs = cpco;};
};

//! Class to handle a bootstrap network.
class BootstrapNetworkTaskHelp
{
private:

	list<DiscreteData *> origDiscreteData;
	list<DiscreteData *> bootDiscreteData;
	list<CtsData *> origCtsData;
	list<CtsData *> bootCtsData;

	map<unsigned int, unsigned int> dataForBootstrap; //ordinal number, data no in list
	unsigned int totalNonMissingData;
	unsigned int totalMissingData;
	map<unsigned int, unsigned int> noTimesDataForBootstrap; // data no in list, no times randomly chosen for bootstrap

	string famFilename;
	unsigned int noFamiliesNotMiss;
	map<unsigned int, string> familyIDs;
	
	bool useMeasurementError;
	bool useEffectSize;
	bool useBootstraps;
	//double measurementError;
	string measurementErrorFile;
	list<double> measurementErrors;

public:

	BootstrapNetworkTaskHelp() : noFamiliesNotMiss(0), totalMissingData(0), totalNonMissingData(0), useMeasurementError(false), useEffectSize(false), useBootstraps(true), measurementErrors(), measurementErrorFile("") {};

	~BootstrapNetworkTaskHelp() {};

	Network * createBootstrapNetwork(TaskCentral * taskCentral, Network * network, const string & name);
	void setupBootstrap(TaskCentral * taskCentral, Network * bootstrapNetwork,  Network * network, const bool & useRandomData = false, const unsigned int & percent = 0);
	void updateBootstrapData(Network * bootstrapNetwork, Network * network);
	void updateBootstrapDataSubset(Network * bootstrapNetwork, Network * network, const bool & useRandomData, const unsigned int & percent);
	void setupFamilyIDs(const string & famFile, TaskCentral * taskCentral, Network * network);
	void setMeasurementErrors(const map<unsigned int, double> & mes, const double & me, const map<unsigned int, double> & mesm, const double & mem, Network * network);
	void setMeasurementErrorFile(const string & mef) { measurementErrorFile = mef; useMeasurementError = true; };
	void setUseBootstraps(const bool & ub) {useBootstraps = ub;};
	void setUseEffectSize(const bool & ues) { useEffectSize = ues; };
	void setupMEStDevs(const double & measurementError);
	void setupMEStDevMultiples(const double & measurementErrorMultiple);
	list<double> getMeasurementErrors() {return measurementErrors;};
	list<CtsData *> getOrigCtsData() {return origCtsData;};
};

//! Calculate average networks.
class AverageNetworksTask : public Task
{
private:

	string networkName;
	string filename;
	string igraphPrefix;
	Network * network;
	Network * bootstrapNetwork;
	string famFilename;
	string likelihoodFilename;

	double arcThreshold;
	bool arcThresholdSet;
	string thresholdFilename;
	bool useEquivNets;
	unsigned int noBootstraps;
	unsigned int randomRestarts;
	unsigned int jitterRestarts;
	
	BootstrapNetworkTaskHelp bootstrapNetworkTaskHelp;
	map<string, pair<double, double> > arcStrengths; //"node1 node2", arc strength, arc direction
	map<string, double> arcCounts; //number of times an arc apeears in any direction
	map<string, pair<double, double> > networkArcStrengths;
	list<double> orderedArcStrengths;

	bool useNetworkScoreMethod;
	bool useNetworkWeightMethod;
	double networkProb;
	double totalProb;
	unsigned int noNetworksEval;
	double offSet;
	bool offSetSet;
	bool useMeasurementError;
	bool useBootstraps;
	bool useEffectSize;
	double measurementError;
	string measurementErrorFile;
	map<unsigned int, double> measurementErrors;	
	double measurementErrorMultiple;
	string measurementErrorMultipleFile;
	map<unsigned int, double> measurementErrorMultiples;
	map<pair<unsigned int, unsigned int>, double> fromToProb; //from node, to node, prob. for robustness output

public:

	AverageNetworksTask() : filename("averageNetwork.dat"), thresholdFilename(""), igraphPrefix(""), arcThreshold(0), arcThresholdSet(false), useEquivNets(false), noBootstraps(100), randomRestarts(0), jitterRestarts(0), bootstrapNetworkTaskHelp(), useNetworkScoreMethod(false), useNetworkWeightMethod(false), famFilename(""), likelihoodFilename(""),
	   bootstrapNetwork(0), network(0), networkProb(0), noNetworksEval(0), offSet(0), offSetSet(false), totalProb(0), useMeasurementError(false), useBootstraps(true), useEffectSize(false), measurementError(0), measurementErrorFile(""), measurementErrors(), measurementErrorMultiple(), measurementErrorMultipleFile(""), measurementErrorMultiples() {};

	~AverageNetworksTask()
	{
		
	};

	//general task methods
	void initialiseTask();
	void outputTaskHeader();
	void outputTaskDetails();
	void doTask();
	
	//average network methods
	void setNetworkName(const string & nm) {networkName = nm;};
	void setFilename(const string & fn) {filename = fn;};
	void setThresholdFilename(const string & fn) {thresholdFilename = fn;};
	void setFamilyFile(const string & ff) {famFilename = ff;};
	void setArcThreshold(const double & th) {arcThreshold = th; arcThresholdSet = true;};
	double getArcThreshold() {return arcThreshold;};
	void setNoBootstraps(const unsigned int & nb) {noBootstraps = nb;};
	void setUseEquivNets(const bool & dqc) { useEquivNets = dqc; };
	void setRandomRestarts(const unsigned int & rr) {randomRestarts = rr;};
	void setJitterRestarts(const unsigned int & jr) {jitterRestarts = jr;};
	void setRFilenamePrefix(const string & rfn) {igraphPrefix = rfn;};
	void setLikelihoodFilename(const string & fn) {likelihoodFilename = fn;};
	void setMeasurementError(const double & me) {measurementError = me; useMeasurementError = true;};
	void setMeasurementErrors(const map<unsigned int, double> & mes)
	{ 
		for(map<unsigned int, double>::const_iterator m = mes.begin(); m != mes.end(); ++m)
		{
			measurementErrors[m->first] = m->second;
		};
		useMeasurementError = true;
	};
	void setMeasurementErrorFile(const string & mef) { measurementErrorFile = mef; useMeasurementError = true; };
	void setMeasurementErrorMultiple(const double & me) { measurementErrorMultiple = me; useMeasurementError = true; };
	void setMeasurementErrorMultiples(const map<unsigned int, double> & mes)
	{
		for(map<unsigned int, double>::const_iterator m = mes.begin(); m != mes.end(); ++m)
		{
			measurementErrorMultiples[m->first] = m->second;
		};
		useMeasurementError = true;
	};
	void setMeasurementErrorMultipleFile(const string & mef) { measurementErrorMultipleFile = mef; useMeasurementError = true; };
	void setUseMeasurementError(const bool & um) { useMeasurementError = um; };
	void setUseEffectSize(const bool & ues) { useEffectSize = ues; };
	void outputNodeRobustnessSummary();
	void outputEdgeRobustnessSummary();
	void getNodeAveEdgeProbs(string & nodeName, double & aveParentProb, double & aveChildProb, double & aveProb);
	void setUseBootstraps(const bool & ub) {useBootstraps = ub;};
	void updateArcStrengths(Network * net);
	void updateArcStrengthsNetwork(const string & networkStr, const string & preNodeName = "");
	void updateArcStrengthsTotal();
	void estimateArcThreshold();
	void setFinalNetwork();
	void outputRGraph();
	void outputAverageNetworkFile();

	//average methods for So-Youn approach using network scores
	void setUseNetworkScoreMethod(const bool & u) {useNetworkScoreMethod = u;};
	void doNetworkScoreMethod();
	void setOffSetScoreMethod();
	void updateArcStrengthsScoreMethod();

	//average methods using weighted bootstrap
	void setUseNetworkWeightMethod(const bool & u) {useNetworkWeightMethod = u;};
};


//! Compare some networks.
class CompareNetworksTask : public Task
{
private:

	
	
public:

	CompareNetworksTask() {};

	~CompareNetworksTask() {};

	//general task methods
	void initialiseTask();
	void outputTaskHeader();
	void outputTaskDetails();
	void doTask();
	
};

//! Simulate a network.
class SimulateNetworkDataTask : public Task
{
private:
	
	string networkName;
	Network * network;
	NetworkBNLearn * simNetwork;
	string paraFilename;
	unsigned int noSims;
	string whitelistFilename;
	string blacklistFilename;
	string score;
	string scoreFix;
	map<string, unsigned int> noLevels; //node name, no levels; no of levels for discrete node
	set<string> nodesWithDiscreteParents;
	bool createDifferentNodeData; //if creating network from existing network with data then if true creates different data with different name
	bool snpDataIntegers;
	map<string, unsigned int> missingFirst;
	map<string, unsigned int> missingLast;
	map<string, double> missingProb;

	unsigned int totalDisNodes;
	unsigned int totalFactorNodes;
	unsigned int totalCtsNodes;
	unsigned int totalNoDataNodes;
	unsigned int totalEdges;
	double imaginarySampleSize;

public:

	SimulateNetworkDataTask() : networkName(""), network(0), simNetwork(0), paraFilename(""), noSims(100), whitelistFilename(""), blacklistFilename(""), score("BIC"), scoreFix(""), createDifferentNodeData(false),
		totalDisNodes(0), totalFactorNodes(0), totalCtsNodes(0), totalEdges(0), imaginarySampleSize(10), totalNoDataNodes(0), snpDataIntegers(false) {};

	~SimulateNetworkDataTask() {};

	//general task methods
	void initialiseTask();
	void outputTaskHeader();
	void outputTaskDetails();
	void doTask();
	
	// simulate data task methods
	void setNetworkName(const string & nm) {networkName = nm;};	
	void setParameterFile(const string & pf) {paraFilename = pf;};
	void setNoSims(const unsigned int ns) {noSims = ns;};
	void setSNPIntegers(const bool & si) {snpDataIntegers = si;};
	void setCreateDifferentNodeData(const bool & cdnd) {createDifferentNodeData = cdnd;};
	void createSimNetworkFromNetwork();
	void createSimNetworkFromParaFile();	
	void addNodesToSimNetwork();
	void addParentsToSimNetwork();
	void addParasToSimNetwork();
	bool getHasDiscreteParents(const string & nd);
	void findNoLevelsFromParasFile();
	void addDiscreteSimParas(Node * simNode, ifstream & readParaFile);
	void addCtsSimParas(Node * simNode, ifstream & readParaFile, string & word);
	void setImaginarySampleSize(const double & i) {imaginarySampleSize = i;};	
	void setScore(const string & sc) {score = sc;};
	void setMissingFirst(const string & v, const unsigned int & i) {missingFirst[v] = i;};
	void setMissingLast(const string & v, const unsigned int & i) {missingLast[v] = i;};
	void setMissingProb(const string & v, const double & pr) {missingProb[v] = pr;};
	void setScoreFix(const string & scf) {score = scf;};
	void addWhitelistDataFile(const string & file) {whitelistFilename = file;};
	void addBlacklistDataFile(const string & file) {blacklistFilename = file;};
};

//! Outputs the structure of a network.
class OutputNetworkTask : public Task
{
private:

	string networkName;
	Network * network;
	string filename;  //standard BayesNetty format
	string filename2; //bnlearn style, [a][b][b|a:b]
	string filenamePrefix; //for igraph
	string nodeDataFilePrefix;
	string equivFilename;
	unsigned int noEquivNets;
	bool outputBedFile;

	bool wasDisData;
	bool wasCtsData;
	bool wasSNPDisData;
	bool wasSNPCtsData;
	bool wasSNPData;
	bool loopInNetwork;
	unsigned int jobNo;
	unsigned int jobTotal;
	unsigned int startIndivNo;
	unsigned int endIndivNo;

public:

	OutputNetworkTask() : network(0), filename(""), filename2(""), filenamePrefix(""), nodeDataFilePrefix(""), equivFilename(""), outputBedFile(false), wasDisData(false), wasCtsData(false), wasSNPDisData(false), wasSNPCtsData(false), wasSNPData(false),
		noEquivNets(0), loopInNetwork(0), jobNo(0), jobTotal(0), startIndivNo(0), endIndivNo(0) {};

	~OutputNetworkTask() {};

	//general task methods
	void initialiseTask();
	void outputTaskHeader();
	void outputTaskDetails();
	void doTask();
	
	//output network methods
	void setNetworkName(const string & nm) {networkName = nm;};
	void setFilenamePrefix(const string & fn) {filenamePrefix = fn;};
	void setFilename(const string & fn) {filename = fn;};
	void setFilename2(const string & fn) {filename2 = fn;};
	void setNodeDataFilePrefix(const string & fp) {nodeDataFilePrefix = fp;};
	void setOutputBedFile(const bool & ob) {outputBedFile = ob;};
	void setEquivNetsFilename(const string & eqfn) {equivFilename = eqfn;};
	void outputNetworkIGraph();
	void outputNetwork();
	void outputNetwork2();
	void outputNodeData();
	void outputEquivalentNetworks();
	void setStartIndivNo(const unsigned int & sn) {startIndivNo = sn;};
	void setEndIndivNo(const unsigned int & en) {endIndivNo = en;};
	void setJobNo(const unsigned int & jn) {jobNo = jn;};
	void setJobTotal(const unsigned int & jt) {jobTotal = jt;};
	
};


//! Outputs the priors of a network.
class OutputPriorsTask : public Task
{
private:

	string networkName;
	Network * network;
	string filename;

public:

	OutputPriorsTask() : network(0), filename("priors.dat") {};

	~OutputPriorsTask() {};

	//general task methods
	void initialiseTask();
	void outputTaskHeader();
	void outputTaskDetails();
	void doTask();
	
	//output network methods
	void setNetworkName(const string & nm) {networkName = nm;};
	void setFilename(const string & fn) {filename = fn;};
	
};

//! Outputs the posteriors of a network.
class OutputPosteriorsTask : public Task
{
private:

	string networkName;
	Network * network;
	Network * bootstrapNetwork;
	string filename;

public:

	OutputPosteriorsTask() : network(0), bootstrapNetwork(0), filename("posteriors.dat") {};

	~OutputPosteriorsTask() {};

	//general task methods
	void initialiseTask();
	void outputTaskHeader();
	void outputTaskDetails();
	void doTask();
	
	//output network methods
	void setNetworkName(const string & nm) {networkName = nm;};
	void setFilename(const string & fn) {filename = fn;};
	
};

//! Imputes data for a network.
class ImputeNetworkDataTask : public Task
{
private:

	string networkName;
	Network * network;
	Network * bootstrapNetwork;	
	bool usePrevNetwork; //use previous network as starting point in bootstrap search
	bool useCompleteData;
	bool useRandomData;
	bool useMean; //just substitute mean for cts variables - just for evalution purposes
	bool useAllNN; //use all variables for NN - just for evalution purposes
	unsigned int minIndivsForCompleteData; //if training data method not set then this is the min amount to auto set to complete data
	unsigned int randomRestarts;
	unsigned int jitterRestarts;
	unsigned int noBitsOfDataToImpute;
	unsigned int noBitsOfDataNotImputed;
	unsigned int noBitsOfDataNotImputed1, noBitsOfDataNotImputed2, noBitsOfDataNotImputed3;
	unsigned int impStartIndivNo;
	unsigned int impEndIndivNo;
	unsigned int jobNo;
	unsigned int jobTotal;
	unsigned int maxMissing;
	unsigned int noIndivsWithMissing;
	unsigned int noIndivsImputed;
	bool useMinNonMissingEdges;
	double minNonMissingEdges;
	unsigned int subsetPercent;

public:

	ImputeNetworkDataTask() : usePrevNetwork(false), useCompleteData(false), useRandomData(false), useMean(false), useAllNN(false), minIndivsForCompleteData(40), randomRestarts(0), jitterRestarts(0), impStartIndivNo(0), impEndIndivNo(0), jobNo(0), jobTotal(0),
		network(0), bootstrapNetwork(0), maxMissing(0), noIndivsWithMissing(0), noIndivsImputed(0), useMinNonMissingEdges(true), minNonMissingEdges(0), noBitsOfDataNotImputed(0), noBitsOfDataNotImputed1(0), noBitsOfDataNotImputed2(0), noBitsOfDataNotImputed3(0), noBitsOfDataToImpute(0), subsetPercent(90) {};

	~ImputeNetworkDataTask() {};

	//general task methods
	void initialiseTask();
	void outputTaskHeader();
	void outputTaskDetails();
	void doTask();
	
	//impute methods
	void setNetworkName(const string & nm) {networkName = nm;};
	void doNetworkDataImputation(const bool & doNotDoAdjust = false, const bool & updateImpImmed = false, const unsigned int & tryNo = 1);
	void setRandomRestarts(const unsigned int & rr) {randomRestarts = rr;};
	void setJitterRestarts(const unsigned int & jr) {jitterRestarts = jr;};
	void setStartIndivNo(const unsigned int & sn) {impStartIndivNo = sn;};
	void setEndIndivNo(const unsigned int & en) {impEndIndivNo = en;};
	void setJobNo(const unsigned int & jn) {jobNo = jn;};
	void setJobTotal(const unsigned int & jt) {jobTotal = jt;};
	void setMaxMissing(const unsigned int & mm) {maxMissing = mm;};
	void setMinNonMissingEdges(const double & mnme) {minNonMissingEdges = mnme/100.0; if(minNonMissingEdges == 0) useMinNonMissingEdges = false;};
	void setSubsetPercent(const unsigned int & sp) {subsetPercent = sp;};
	void setUseCompleteData() {useCompleteData = true; useRandomData = false;};
	void setUseRandomData() {useRandomData = true; useCompleteData = false;};
	void setUseMean() {useMean = true; useCompleteData = true; useRandomData = false;};
	void setUseAllNN() {useAllNN = true;};
	void setUsePrevNetwork() {usePrevNetwork = true;};
};

//! Estimates the recall and precision for a network and data.
class ImputeEstimateRecallPrecisionTask : public Task
{
private:
	string networkName;
	Network * network;
	Network * simNetwork;
	unsigned int randomRestarts;
	unsigned int jitterRestarts;
	double recallNoImp;
	double preNoImp;
	double recallImp;
	double preImp;
	double recallImpRT;
	double preImpRT;
	double recallFull;
	double preFull;
	pair<unsigned int, unsigned int> missingNotMissing;
	double minNonMissingEdges;
	unsigned int iterations;
	bool doImps;
	string simDataFilename;
	string simNetFilename;

public:

	ImputeEstimateRecallPrecisionTask() : network(0), simNetwork(0), randomRestarts(0), jitterRestarts(0), recallNoImp(0),	preNoImp(0), recallImp(0), preImp(0),
		recallImpRT(0), preImpRT(0), recallFull(0), preFull(0), minNonMissingEdges(0), iterations(1), doImps(true), simDataFilename(""), simNetFilename("") {};

	~ImputeEstimateRecallPrecisionTask() {};

	//general task methods
	void initialiseTask();
	void outputTaskHeader();
	void outputTaskDetails();
	void doTask();

	//est recall and pre impute methods
	void setNetworkName(const string & nm) {networkName = nm;};
	void setRandomRestarts(const unsigned int & rr) {randomRestarts = rr;};
	void setJitterRestarts(const unsigned int & jr) {jitterRestarts = jr;};
	void setIterations(const unsigned int & i) {iterations = i;};
	void setDoImps(const bool & di) {doImps = di;};	
	void setSimDataFilename(const string & s) {simDataFilename = s;};
	void outputSimData(const unsigned int & it);
	void setSimNetworkFilename(const string & s) {simNetFilename = s;};
	void outputSimNetwork(const unsigned int & it);
};

//! Calculates the recall and precision for a network and another network.
class CalculateRecallPrecisionTask : public Task
{
private:
	string networkName;
	Network * network;
	string trueNetworkName;
	Network * trueNetwork;
	double recall;
	double precision;
	string filename;

public:

	CalculateRecallPrecisionTask() : network(0), recall(0),	precision(0), trueNetwork(0), trueNetworkName(""), filename("") {};

	~CalculateRecallPrecisionTask() {};

	//general task methods
	void initialiseTask();
	void outputTaskHeader();
	void outputTaskDetails();
	void doTask();

	//est recall and pre impute methods
	void setNetworkName(const string & nm) {networkName = nm;};
	void setTrueNetworkName(const string & tnm) {trueNetworkName = tnm;};
	void setFilename(const string & fnm) {filename = fnm;};
};

//! Analysis the Robustness of the network fitted where data has Measurement Error.
class MeasurementErrorRobustnessTask : public Task
{
private:	
	string networkName;
	string filename;
	string thresholdFilename;
	string igraphPrefix;
	Network * network;
	Network * bootstrapNetwork;
	string famFilename;
	AverageNetworksTask averageNetworksTask;

	double arcThreshold;
	bool arcThresholdSet;
	bool useBootstraps;
	unsigned int noIterations;
	unsigned int randomRestarts;
	unsigned int jitterRestarts;
	bool useEquivNets;

	double measurementError;
	double measurementErrorMultiple;
	map<unsigned int, double> measurementErrors;
	map<unsigned int, double> measurementErrorMultiples;
	map<string, double> measurementStrErrors;
	map<string, double> measurementStrErrorMultiples;
	string measurementErrorFile;
	string measurementErrorMultipleFile;
	bool useEffectSize;

public:

	MeasurementErrorRobustnessTask() : networkName(""), filename("robustness"), thresholdFilename(""), igraphPrefix(""), arcThreshold(0), arcThresholdSet(false), useBootstraps(false), noIterations(100), randomRestarts(0), jitterRestarts(0), useEquivNets(false), famFilename(""),
		bootstrapNetwork(0), network(0), measurementError(0), measurementErrorMultiple(0), measurementStrErrors(), measurementErrorFile(""), measurementErrorMultipleFile(""), useEffectSize(false) {};

	~MeasurementErrorRobustnessTask() {};

	//general task methods
	void initialiseTask();
	void outputTaskHeader();
	void outputTaskDetails();
	void doTask();

	//Measurement Error Robustness Task methods	
	void setNetworkName(const string & nm) { networkName = nm; };
	void setFilename(const string & fn) { filename = fn; };
	void setThresholdFilename(const string & fn) { thresholdFilename = fn; };
	void setFamilyFile(const string & ff) { famFilename = ff; };
	void setArcThreshold(const double & th) { arcThreshold = th; arcThresholdSet = true; };
	void setUseBootstraps(const bool & ub) { useBootstraps = ub; };
	void setNoIterations(const unsigned int & nb) { noIterations = nb; };
	void setRandomRestarts(const unsigned int & rr) { randomRestarts = rr; };
	void setJitterRestarts(const unsigned int & jr) { jitterRestarts = jr; };
	void setRFilenamePrefix(const string & rfn) { igraphPrefix = rfn; };
	void setMeasurementError(const double & me) { measurementError = me; };
	void setMeasurementError(const string & nd, const double & me);
	void setMeasurementErrorFile(const string & mef) { measurementErrorFile = mef; };
	void setMeasurementErrorMultiple(const double & mem) { measurementErrorMultiple = mem; };
	void setEffectSizeError(const bool & e) { useEffectSize = e; };
	void setMeasurementErrorMultiple(const string & nd, const double & me);
	void setMeasurementErrorMultipleFile(const string & mef) { measurementErrorMultipleFile = mef; };	
	void setUseEquivNets(const bool & dqc) { useEquivNets = dqc; };
};

#endif
