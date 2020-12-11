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


/*! \file Nodes.h
    \brief This file organises the nodes of a network.
    
*/

#ifndef __NODES
#define __NODES

#include "Data.h"
#include "Model.h"

#include <string>
#include <list>
#include <set>

class Network;
class NetworkDeal;
struct CtsMultiDistJoint;

//! Results of linear regressions from cts nodes.
struct CtsPriLocalDist {

	double intercept; //cts node id, mean
	map<unsigned int, double> coeffs; //cts node id of parent node, linear regression coeff for parent node. 
	double variance;  //fitted variance from linear regression
	double stDev;
	double mean;

	CtsPriLocalDist() : intercept(0), variance(0), stDev(0), mean(0) {};

	~CtsPriLocalDist() {};
};

//! Master multivariate distribution used in deal.
struct CtsMultiDistMaster {
	list<double> means; //means, ordered s.t. node then parents in order of node ID
	list< list<double> > covariances; //(node ID then parent node IDs) covarianace matrix 
	double jointDiscretePriorParaSum;

	CtsMultiDistMaster() : jointDiscretePriorParaSum(0) {};

	~CtsMultiDistMaster() {};
};

//! Local parameters used in deal.
struct CtsLocalPara {

	list<double> mu;
	double phi;
	double rho;
	list< list<double> > tau;
	list< list<double> > tauInverse;

	CtsLocalPara(): phi(0), rho(0) {};

	~CtsLocalPara()
	{
		
	};
};

//! A node with its parents.
class Node
{
private:

	string name; //column name of data 
	string displayName; //to allow name to appear the same for simulating data
	
	bool visited; //for detecting loops
	map<unsigned int, pair<double, unsigned int> > edgeSignifs; //parent node ID, chi sq, df
	
protected:

	unsigned int nodeID;
	map<unsigned int, Node *> parents;
	map<unsigned int, string> localParaNames; //discrete parent node group, name
	map<unsigned int, string> nodeParentNames; //discrete node+parent node group, name	
	bool isDiscreteNode;
	bool warningDone;

	NetworkMissingData * networkMissingData; //keeps track of if data should be considered as missing in the network

public:

	Node(const string & nm) : name(nm), displayName(""), visited(false), nodeID(0), isDiscreteNode(false), warningDone(false), networkMissingData(0)
	{ 
		
	};

	//! Delete 
	virtual ~Node()
	{
		
	};

	//virtual methods
	virtual unsigned int getFileNo() {return 0;}; 
	virtual void setDefaultInitPrior() {};
	virtual unsigned int getNoLevels() const {return 1;};
	virtual unsigned int getLevel() {return 1;};
	virtual unsigned int getRandomLevel() {return 1;};

	virtual void setDataValue() {};
	virtual void setInitialDataValue() {};
	virtual void nextDataValueIter() {};
	virtual bool nodeDataIsMissing() {return false;};
	virtual void setInitialDataValue2() {}; //second versions return missing data on a nodes rather than network wide
	virtual void nextDataValueIter2() {};
	virtual bool nodeDataIsMissing2() {return false;};
	virtual double getNodeDataCtsValue2() {return 0;};
	virtual unsigned int getNodeDataDisValue2() {return 1;};
	virtual void setInitialDataValue3() {}; //third version for another loop thro the data
	virtual void nextDataValueIter3() {};
	virtual bool nodeDataIsMissing3() {return false;};
	virtual double getNodeDataCtsValue3() {return 0;};
	virtual unsigned int getNodeDataDisValue3() {return 1;};
	virtual void setImputedDataCts(const double & val, const bool & updateImpImmed) {};
	virtual void setImputedDataDis(const unsigned int & val, const bool & updateImpImmed) {};
	virtual void setStDev() {};
	virtual double getStDev() {return 0;};
	virtual double getMean() {return 0;};
	virtual double getRandomValue() {return 0;};
	virtual double getNodeDataCtsValueAdj(const bool & impIndivOrNN, const set<unsigned int> & nodesOfInterest, const list<unsigned int> & bootNodeIDsAdjVars, const list<unsigned int> & adjMethods, Network * network, Network * bootNetwork, const bool & doRandomDisGroup, const unsigned int & grp) {return 0;};	
	virtual double getNodeDataStDevAdj(const set<unsigned int> & nodesOfInterest, const list<unsigned int> & bootNodeIDsAdjVars, const list<unsigned int> & adjMethods, Network * network, Network * bootstrapNetwork) {return 0;};
	virtual void randomlyFillMissingValues() {};	
	virtual void setupWasImputed() {};
	virtual void updateImputedDataAsNonMissing() {};

	virtual CtsPriLocalDist * getCtsPriLocalDist(const unsigned int & groupNodeID) {return 0;};
	virtual void setLevel(const unsigned int & lv) {};
	virtual void setLevel(const string & lv) {};
	virtual void setFixed(const bool & fix) {};
	virtual string getLevelName(const unsigned int & lv = 0) {return "";};
	virtual unsigned int getLevelNo(const string & strLevel) {return 0;};
	virtual void setIsFactorNode(const unsigned int & b) {};
	virtual void setIsChildFactorNode(const unsigned int & b) {};

	virtual double calcScoreBit(Network * network) { return 0; };
	virtual double calcUpdatedScoreBit(Network * network) {return 0;};
	virtual string getNodeType() {return "NA";};
	
	virtual void updateNotInGroupOrMissing(list<bool> & notInGroupOrMissingData) {};
	virtual void updateNetworkMissingData(NetworkMissingData * networkMissingData) {};
	virtual unsigned int getAmountData() const {return 0;};	
	virtual CtsData * getCtsData() {return 0;};
	virtual DiscreteData * getDiscreteData() {return 0;};
	virtual Data * getData() {return 0;};
	virtual bool getIsFactorNode() const {return false;};
	virtual bool getIsFactorChildNode() const {return false;};
	virtual string getFactorNodeName() const {return "";};
	virtual string getFactorName() const {return "";};
	virtual void clearCache() {};
	
	virtual void updateCtsJointProbDist(CtsMultiDistJoint *) {};

	virtual void clearPriorsAndPosteriors() {};
	virtual bool getIsSNPNode() {return false;};
	virtual void setIsSNPNode(const bool & isn) {};

	virtual void simData(bool & snpDataIntegers) {};
	virtual void simDataCopyParas(const string & simTaskName, Node * copyingNode, Network * origNetwork, Network * copyingNetwork, const bool & createDifferentNodeData) {};
	virtual void simDataCopyDisLevels(Node * copyingNode) {};
	virtual unsigned int getSimNodeType() {return 0;};
	virtual unsigned int getSimNoLevels() const {return 1;};	
	virtual void setSimNoLevels(const int & nl) {};
	virtual void setDefaultParameters() {};
	virtual map<string, unsigned int> getLevels() {return map<string, unsigned int>();};
	virtual void setCtsParas(map<unsigned int, CtsPriLocalDist *> & ctsParas, const string & simTaskName, Network * origNetwork, Network * copyingNetwork, const bool & createDifferentNodeData) {};
	virtual void setCtsParas(unsigned int & nodeParentGroupNo, CtsPriLocalDist * ctsParas) {};
	virtual void setDiscreteParas(map<unsigned int, unsigned int> & ndDataCnts, map<unsigned int, unsigned int> & parentDataCnts, unsigned int & totData) {};
	virtual void setSimDataCumLevelProb(unsigned int & parentGroupLevel, map<unsigned int, double> & probs) {};
	virtual void setValue(const double & v) {};
	virtual double getValue() const {return 0;};
	//virtual void outputPriorLocalProbDistsTEST() {};

	//general node methods
	string getName() {return name;}; 
	string getDisplayName() {if(displayName!="") return displayName; else if(getIsFactorChildNode() || getIsFactorNode()) return getFactorNodeName(); else return name;};
	bool getIsDiscreteNode() const {return isDiscreteNode;};
	void setDisplayName(const string & nm) {displayName = nm;};
	unsigned int getNodeID() {return nodeID;};
	void setNodeID(const unsigned int & no) {nodeID = no;};
	void setNetworkMissingData(NetworkMissingData * nmd) {networkMissingData = nmd;};
	NetworkMissingData * getNetworkMissingData() {return networkMissingData;};
	void updateMissingDataAtNodeLevel();
	void addParent(Node * parent) {parents[parent->getNodeID()] = parent;};
	void copyParents(Node * copyingNode, Network * copyingNetwork);
	void simDataCopyParents(const string & simTaskName, Node * copyingNode, Network * copyingNetwork, const bool & createDifferentNodeData);
	void copyMissingness(Node * otherNode);

	void setInitialLevels(); //set node and parent levels to 1
	bool advanceAllLevels(); //advance node and parent levels
	bool advanceParentLevels(); //advance parent levels

	bool parentExists(const unsigned int & no) {return parents.find(no) != parents.end();};
	void removeParent(const unsigned int & nodeNo);
	void removeAllParents() { if(freeClearMemory) map<unsigned int, Node*>().swap(parents); else parents.clear(); };
	unsigned int getNoParents() {return parents.size();};
	map<unsigned int, Node *> getParents() {return parents;};
	unsigned int getParentsGroupNo();
	unsigned int getNodeAndParentsGroupNo();
	unsigned int getLevelFromNodeAndParentsGroupNo(const unsigned int & groupNo);
	unsigned int getNodeAndParentsGroupNo2(Network * network);
	unsigned int getNodeAndParentsGroupNo3(Network * network);
	unsigned int getRandomNodeAndParentsGroup2(Network * network);
	unsigned int getRandomNodeAndParentsGroup3(Network * network);
	string calcDiscreteParentsName();
	void calcLocalParaNames();
	void calcDiscreteNodeParentNames();
	string getDiscreteParentName(const bool & includeNode = false);
	void outDiscreteParents(ofstream & outFile, const bool & includeNode = false);
	pair<unsigned int, unsigned int> getNoMissingNotMissing();

	//double getUpdatedScoreBit(Network * network);
	void setVisited(const bool & v) {visited =v;};
	bool getVisited() {return visited;};
	bool allParentsVisited();
	string getCondName();
	void outputEdges(ofstream & outfile, const unsigned int & nodeID);
	void outputEdgeNames(ofstream & outfile);
	void calcEdgeSignifs(Network * network);
	double getEdgeProb(const unsigned int & parentNodeID);
	pair<double, unsigned int> getEdgeSign(const unsigned int & parentNodeID);	
	virtual unsigned int getNumberOfParameters() {return 0;};

	//used for testing
	//virtual void outputPriorsTEST() {};
};

//! Node with no data, used when simulating data.
class NoDataNode : public Node
{
private:
	
	unsigned int simNodeType; //0 = cts, 1 = discrete
	unsigned int simNoLevels; //no of levels if node will become discrete 

protected:

	NoData * data;

public:

	NoDataNode(const string & nm, NoData * cd) : Node(nm), data(cd), simNodeType(0), simNoLevels(1)
	{ 
		isDiscreteNode = false;
	};

	//! Delete 
	virtual ~NoDataNode()
	{		
				
	};

	//no data node methods
	unsigned int getFileNo() {return data->fileNo;}; 
	void setSimNodeType(const unsigned int & nt) {simNodeType = nt;};	
	void setSimNoLevels(const int & nl) {simNoLevels = nl;};
	unsigned int getSimNoLevels() const {return simNoLevels;};

	//common node methods	
	string getNodeType() {return "n";};
	unsigned int getSimNodeType() {return simNodeType;};
	void setIsSNPNode(const bool & isn) {data->isSNPData = isn;};
	bool getIsSNPNode() {return data->isSNPData;};
	Data * getData() {return data;};
};

//! Cts Node.
class CtsNode : public Node
{
private:
	

protected:

	//local priors for every poss. parents setting, parent group no. given like mileometer
	//Dist given by Normal dist
	map<unsigned int, CtsPriLocalDist *> ctsPriLocalDists; //node(always group 1 as cts) + parents group number;  prob 
	map<string, map<unsigned int, CtsPriLocalDist *> > ctsPriLocalDistsCache;
	list<CtsPriLocalDist *> ctsPriLocalDistsToDel;

	list<unsigned int> nodeParentIDsForData; //nodeParent ID corresponding to the cts data, used in post. calc., set in local dist prob calc

	double value;  //set for calc'ing stuff
	CtsData * data;

	list<bool>::const_iterator missIter; //network wide

	list<double>::iterator dataIter2; //for iteratoring thro' data
	list<bool>::iterator missIter2; //just for this node
	list<bool>::iterator wasImpIter2; //just for this node
	list<double>::const_iterator dataIter3; //for iteratoring thro' data
	list<bool>::const_iterator missIter3; //just for this node
	double stDev; //for normalizing NN distances
	bool meanSet;
	double mean; //for adjusting NN distances when missing data probs

public:

	CtsNode(const string & nm, CtsData * cd) : Node(nm), data(cd), value(0), stDev(0), mean(0), meanSet(false)
	{ 
		isDiscreteNode = false;
	};

	//! Delete 
	virtual ~CtsNode()
	{			
		for(list<CtsPriLocalDist *>::iterator i = ctsPriLocalDistsToDel.begin(); i != ctsPriLocalDistsToDel.end(); ++i)
		{
			delete *i;
		};
	};

	//cts node methods
	unsigned int getFileNo() {return data->fileNo;}; 
	void clearLocalPriors() { if(freeClearMemory) map<unsigned int, CtsPriLocalDist*>().swap(ctsPriLocalDists); else ctsPriLocalDists.clear(); };
	void setValue(const double & v) {value = v;};
	double getValue() const {return value;};

	//common node methods	
	string getNodeType() {return "c";};
	unsigned int getNumberOfParameters();
	bool getIsSNPNode() {return data->isSNPData;};
	void setIsSNPNode(const bool & isn) {data->isSNPData = isn;};
	void setDefaultInitPrior();
	void setParentData(LinearRegModel & linearRegModel);
	void updateParentData(LinearRegModel & linearRegModel);
	void updateNetworkMissingData(NetworkMissingData * networkMissingData);
	unsigned int getAmountData() const {return data->getAmountOfData();};	
	bool getIsFactorNode() const {return data->isFactorNode;};
	bool getIsFactorChildNode() const {return data->isFactorChildNode;};
	void setIsFactorNode(const bool & b) {data->isFactorNode = b;};
	void setIsChildFactorNode(const bool & b) {data->isFactorChildNode = b;};
	Data * getData() {return data;};

	string getFactorNodeName() const {return data->nodeName;};
	string getFactorName() const {return data->factorName;};
	string getFactorLevelName(const unsigned int & lv = 0) {return data->getFactorLevelName(lv);};
	void clearCache() { if(freeClearMemory) map<string, map<unsigned int, CtsPriLocalDist*> >().swap(ctsPriLocalDistsCache); else ctsPriLocalDistsCache.clear(); };

	CtsData * getCtsData() {return data;};
	list<unsigned int> getCtsParents();
	void outputPriorLocalProbDists(ofstream & priorsFile);
	//void outputPriorLocalProbDistsTEST();

	void setInitialDataValue2(); //second versions return missing data on a nodes rather than network wide
	void nextDataValueIter2();
	bool nodeDataIsMissing2();
	double getNodeDataCtsValue2() {return *dataIter2;};
	void setInitialDataValue3();
	void nextDataValueIter3();
	bool nodeDataIsMissing3();
	double getNodeDataCtsValue3() {return *dataIter3;};
	void setStDev();
	double getStDev() {return stDev;};
	double getMean();
	double getRandomValue();
	void setImputedDataCts(const double & val, const bool & updateImpImmed);
	double getNodeDataCtsValueAdj(const bool & impIndivOrNN, const set<unsigned int> & nodesOfInterest, const list<unsigned int> & bootNodeIDsAdjVars, const list<unsigned int> & adjMethods, Network * network, Network * bootNetwork, const bool & doRandomDisGroup, const unsigned int & grp);
	void randomlyFillMissingValues();
	void setupWasImputed() {data->setupWasImputed();};
	void updateImputedDataAsNonMissing() {data->updateImputedDataAsNonMissing();};
	//double getNodeDataCtsValueAdj2(const set<unsigned int> & nodesOfInterest, const set<unsigned int> & bootNodeIDsAdjVars, Network * network, const unsigned int & adjMethod = 0);
	//double getNodeDataCtsValueAdj3(const set<unsigned int> & nodesOfInterest, const set<unsigned int> & bootNodeIDsAdjVars, Network * network, const unsigned int & adjMethod = 0);
	double getNodeDataStDevAdj(const set<unsigned int> & nodesOfInterest, const list<unsigned int> & bootNodeIDsAdjVars, const list<unsigned int> & adjMethods, Network * network, Network * bootstrapNetwork);
	double getNodeDataMeanOrStDevAdj(const set<unsigned int> & nodesOfInterest, const list<unsigned int> & bootNodeIDsAdjVars, const list<unsigned int> & adjMethods, Network * network, Network * bootstrapNetwork, const bool & doStDev, const double & mean);
	CtsPriLocalDist * getCtsPriLocalDist(const unsigned int & groupNodeID);
};

//! Continuous node used in deal.
class CtsDealNode : public CtsNode
{
private:
	map<unsigned int, CtsMultiDistMaster *> ctsMasterPriDists; //node(always group 1 as cts) + parents group number;  prob 

	map<unsigned int, CtsLocalPara * > ctsLocalParaPris;    //node(always group 1 as cts) + parents group number
	map<unsigned int, CtsLocalPara * >  ctsLocalParaPosts;  //node(always group 1 as cts) + parents group number

	
	map<unsigned int, unsigned int> totalInData;
	map<unsigned int, list< list<double> > > parentDataMat;
	map<unsigned int, list< list<double> > > parentDataMatTrans;
	map<unsigned int, list<double> > nodeData;

public:

	CtsDealNode(const string & nm, CtsData * cd) : CtsNode(nm, cd)
	{

	};

	//! Delete 
	virtual ~CtsDealNode()
	{
		for (map<unsigned int, CtsMultiDistMaster *>::iterator cmd = ctsMasterPriDists.begin(); cmd != ctsMasterPriDists.end(); ++cmd) delete cmd->second;
		for (map<unsigned int, CtsLocalPara *>::iterator ppri = ctsLocalParaPris.begin(); ppri != ctsLocalParaPris.end(); ++ppri) delete ppri->second;
		for (map<unsigned int, CtsLocalPara *>::iterator ppost = ctsLocalParaPosts.begin(); ppost != ctsLocalParaPosts.end(); ++ppost) delete ppost->second;
		ctsMasterPriDists.clear();
		ctsLocalParaPris.clear();
		ctsLocalParaPosts.clear();
	};

	//common node methods
	void clearPriorsAndPosteriors();
	double calcUpdatedScoreBit(Network * network);
	double calcScoreBit(Network * network);	

	void updateCtsJointProbDist(CtsMultiDistJoint * ctsJointDist); 

	//cts deal node methods
	void calcCtsMasterPrior(NetworkDeal * network);
	void calcCtsLocalParameterPrior(NetworkDeal * network);
	void calcCtsLocalParameterPost(NetworkDeal * network);
	void outputMasterPriors(ofstream & priorsFile);
	void outputLocalParameterPriors(ofstream & priorsFile);
	void outputLocalParameterPosteriors(ofstream & posteriorsFile);
	CtsLocalPara * getCtsLocalParaPrior(unsigned int & nodeParentsGroup);
	void setupNodeAndParentData();
	
	
	//void outputPriorsTEST();
};

//! Discrete node class.
class DiscreteNode : public Node
{
private:

protected:
	DiscreteData * data;

	list<unsigned int>::iterator dataIter;
	list<bool>::const_iterator missIter;
	list<unsigned int>::iterator dataIter2; //for iteratoring thro' data
	list<bool>::iterator missIter2;
	list<bool>::iterator wasImpIter2; 
	list<unsigned int>::const_iterator dataIter3; //for iteratoring thro' data
	list<bool>::const_iterator missIter3;

	unsigned int level;  //set for calc'ing stuff
	bool fixed; //for calc'ing priors

	//local priors for every poss. parents setting, parent group no. given like mileometer
	double uniformProb;
	
public:

	DiscreteNode(const string & nm, DiscreteData * dd) : Node(nm), data(dd), level(0), fixed(false)
	{ 
		
		isDiscreteNode = true;
		uniformProb = -1; //denotes not set
	};

	
	//! Delete 
	virtual ~DiscreteNode()
	{
		
	};

	unsigned int getFileNo() {return data->fileNo;}; 
	
	//common methods	
	string getNodeType() { return "d";};
	unsigned int getNumberOfParameters();
	bool getIsSNPNode() {return data->isSNPData;};
	void setIsSNPNode(const bool & isn) {data->isSNPData = isn;};
	void setDefaultInitPrior();
	unsigned int getNoLevels() const {return data->levels.size();};	
	unsigned int getLevel() {return level;};
	unsigned int getRandomLevel() {return data->getRandomLevel();};
	Data * getData() {return data;};

	//discrete node methods
	void setLevel(const unsigned int & lv) {level = lv;};
	void setLevel(const string & lv);
	void setFixed(const bool & fix) {fixed = fix;};
	bool getFixed() {return fixed;};
	string getLevelName(const unsigned int & lv = 0) {if(lv==0) return data->getLevelName(level); else return data->getLevelName(lv);};
	unsigned int getLevelNo(const string & strLevel);	
	DiscreteData * getDiscreteData() {return data;};
	void setDataValue();
	void setInitialDataValue();
	void nextDataValueIter();
	bool nodeDataIsMissing();
	void setInitialDataValue2(); //second versions return missing data on a nodes rather than network wide
	void nextDataValueIter2();
	bool nodeDataIsMissing2();
	unsigned int getNodeDataDisValue2(){return *dataIter2;};
	void setInitialDataValue3(); //second versions return missing data on a nodes rather than network wide
	void nextDataValueIter3();
	bool nodeDataIsMissing3();
	unsigned int getNodeDataDisValue3() {return *dataIter3;};
	void setImputedDataDis(const unsigned int & val, const bool & updateImpImmed);
	void randomlyFillMissingValues();
	void setupWasImputed() {data->setupWasImputed();};
	void updateImputedDataAsNonMissing() {data->updateImputedDataAsNonMissing();};
	unsigned int getAmountData() const {return data->getAmountOfData();};
	void updateNotInGroupOrMissing(list<bool> & notInGroupOrMissingData);
	void updateNetworkMissingData(NetworkMissingData * networkMissingData);

	virtual double getPriLocalProb();

	void outputPriorLocalProbDists(ofstream & priorsFile);
};

//! Discrete node for a deal network.
class DiscreteDealNode : public DiscreteNode
{
private:

	map<unsigned int, double> discreteMasterParaPri; //node+parents group number; Dirichlet dist parameter

	//parents group number; (node level; Dirichlet dist parameter)<-- one Dirichlet dist for each parents config
	//-- if no parents, grp = 0
	map<unsigned int, map<unsigned int, double> > discreteLocalParaPris;
	map<unsigned int, map<unsigned int, double> > discreteLocalParaPosts;

public:

	DiscreteDealNode(const string & nm, DiscreteData * dd) : DiscreteNode(nm, dd)
	{

	};

	//! Delete 
	virtual ~DiscreteDealNode()
	{

	};

	//common node methods
	void clearPriorsAndPosteriors();
	double calcScoreBit(Network * network);
	double calcUpdatedScoreBit(Network * network);
	
	double getMasterPriPara();

	//discrete deal node methods
	void outputMasterPriors(ofstream & priorsFile);
	void outputLocalParameterPriors(ofstream & priorsFile);
	void outputLocalParameterPosteriors(ofstream & posteriorsFile);
	void calcDiscreteMasterPrior(NetworkDeal * network);
	void calcDiscreteLocalParameterPrior();
	void calcDiscreteLocalParameterPost();

};

//! Discrete node for bnlearn network.
class DiscreteBNLearnNode : public DiscreteNode
{
private:

	//-- if no parents, grp = 0
	//map<unsigned int, double> discretePost; // node+parents group number, log probability - prior given by uniform constant (for each poss)
	map<unsigned int, unsigned int> nodeParentDataCounts; // node+parents group number, count of data with this config
	map<unsigned int, unsigned int> parentDataCounts; // parents group number, count of data with this config
	unsigned int totalData;

	map<unsigned int, map<unsigned int, double> > simDataCumLevelProbs; //parents group number, level, prob; cached for sim'ing data
	unsigned int simNoLevels; //no of levels if node will become discrete 

public:

	DiscreteBNLearnNode(const string & nm, DiscreteData * dd) : DiscreteNode(nm, dd), simNoLevels(0), totalData(0)
	{
		
	};

	//! Delete 
	virtual ~DiscreteBNLearnNode()
	{

	};

	//common node methods
	void clearPriorsAndPosteriors();
	double calcScoreBit(Network * network);
	double calcUpdatedScoreBit(Network * network);
	void simDataCopyParas(const string & simTaskName, Node * copyingNode, Network * origNetwork, Network * copyingNetwork, const bool & createDifferentNodeData);
	void setSimDataCumLevelProb(unsigned int & parentGroupLevel, map<unsigned int, double> & probs) {simDataCumLevelProbs[parentGroupLevel] = probs;};
	void setDiscreteParas(map<unsigned int, unsigned int> & ndDataCnts, map<unsigned int, unsigned int> & parentDataCnts, unsigned int & totData);
	void simDataCopyDisLevels(Node * copyingNode);
	void setDefaultParameters();
	map<string, unsigned int> getLevels();
	void simData(bool & snpDataIntegers);
	void setSimNoLevels(const int & nl) {simNoLevels = nl;};
	unsigned int getSimNoLevels() const {return simNoLevels;};

	//discrete bnlearn node methods	
	void updateDataCounts();	
	void outputPosteriors(ofstream & posteriorsFile);	
	unsigned int getNodeParentDataCount(unsigned int & nodeParentsGroup);
	unsigned int getParentDataCount(unsigned int & parentsGroup);
};

//! Continuous node for a bnlearn network.
class CtsBNLearnNode : public CtsNode
{
private:
	

public:

	CtsBNLearnNode(const string & nm, CtsData * cd) : CtsNode(nm, cd)
	{
		
	};

	//! Delete 
	virtual ~CtsBNLearnNode()
	{
		
	};

	//common node methods
	void clearPriorsAndPosteriors();
	double calcUpdatedScoreBit(Network * network);
	double calcScoreBit(Network * network);
	void simDataCopyParas(const string & simTaskName, Node * copyingNode, Network * origNetwork, Network * copyingNetwork, const bool & createDifferentNodeData);
	void setCtsParas(map<unsigned int, CtsPriLocalDist *> & ctsParas, const string & simTaskName, Network * origNetwork, Network * copyingNetwork, const bool & createDifferentNodeData);
	void setCtsParas(unsigned int & nodeParentGroupNo, CtsPriLocalDist * ctsParas);
	void setDefaultParameters();
	void simData(bool & snpDataIntegers);

	//used for testing
	//void outputPriorsTEST();
};


#endif

