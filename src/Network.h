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


/*! \file Network.h
    \brief This file organises the network.
    
*/

#ifndef __NETWORK
#define __NETWORK

#include "main.h"
#include "Nodes.h"

#include <string>
#include <set>
#include <stack>

using namespace std; // initiates the "std" or "standard" namespace

//! Multivariate distribution.
struct CtsMultiDistJoint {
	map<unsigned int, double> means; //cts node id, mean
	map<unsigned int, double> variances; //cts node id, variance for this node. No covariances between nodes defined, ie assumed to be 0
	map<unsigned int, map<unsigned int, double> > covariances; //cts node id1, cts node id2, covariance

	string name;
	
	CtsMultiDistJoint() : means(), variances(), covariances(), name("") {};

	~CtsMultiDistJoint() {};

	double getCovariance(const unsigned int &nd1, const unsigned int &nd2);
};

//! Class for storing potential nearest neighbours.
struct ClosebyIndividual {
	double distance;
	map<unsigned int, double> ctsVariableValues; // cts node ID, value
	map<unsigned int, unsigned int> disVariableValues; // discrete node ID, value

	ClosebyIndividual() : distance(0), ctsVariableValues(), disVariableValues() {};

	ClosebyIndividual(const double & dist, map<unsigned int, double> ctsVV, map<unsigned int, unsigned int> & disVV) : distance(dist), ctsVariableValues(ctsVV), disVariableValues(disVV) {};

	~ClosebyIndividual() {};
};

//! Organises the network.
class Network
{
private:

	
	//what must be in network and must not be in network
	set<pair<unsigned int, unsigned int> > allowedEdgeTypes; // fileNo1, fileNo2
	set<unsigned int> whiteNodes;
	set<unsigned int> blackNodes;
	set<pair<unsigned int, unsigned int> > whiteEdges;
	set<pair<unsigned int, unsigned int> > blackEdges;
	set<unsigned int> noParentsNodes;  //these nodes are not allowed parents
	set<unsigned int> noChildrenNodes; //these nodes are not allowed children

	unsigned int scoreType; //0 = loglik (or default), 1 = BIC, 2 = AIC, 3 = bayes with iss - not all available depending on network type, 4 = cost on edges, 5 = BIC with prob edges
	unsigned int scoreFix;  //0 = no fix, 1 ave log lik, 2 skip bad pt 
	map<pair<unsigned int, unsigned int>, double> costEdges; //when score type = cost
	map<pair<unsigned int, unsigned int>, double> costEdgeTypes;

	list<pair<string, string> > nodesSwappableInit; //for equiv networks

	//for calculating all network scores
	list<pair<unsigned int, unsigned int> > nodePairs; // node ID 1, node ID 2
	list<unsigned int> nodePairConnections;  // 1 = no edge; 2 = node1 --> node2; 3 = node2 --> node1

	//for mapping IDs
	map<unsigned int, unsigned int> nodeIDmap;

	bool isImputedData; //set to true if data is imputed for plotting graph later
	bool isBootstrapData; //set to true if using bootstrap data to avoid 0 variance warnings when unsual data is chosen

protected:
	
	multimap<double, pair<unsigned int, unsigned int> > evaluations;
	
	map<unsigned int, Node *> allNetworkNodes;  //all network nodes - edges stored as node parents
	map<unsigned int, DiscreteNode *> discreteNodes; //only discrete nodes
	map<unsigned int, CtsNode *> ctsNodes;
	map<string, Node *> noDataNetworkNodes; //no data, used for sim'ing data

	AllNodeData * allNodeData; //pointer to all node data

	double imaginarySampleSize;

public:

	Network(AllNodeData * ad, const unsigned int & st, const unsigned int & sf) : allNodeData(ad), scoreType(st), scoreFix(sf), isImputedData(false), isBootstrapData(false), imaginarySampleSize(10)
	{ 
		
	};

	//! Delete network things
	virtual ~Network()
	{		
		for(map<unsigned int, Node *>::iterator an = allNetworkNodes.begin(); an != allNetworkNodes.end(); ++an) delete an->second;		
	};

	virtual double calcScore();	
	virtual void calculatePriorsAndPosteriors() {};
	virtual void calculatePriorsAndPosteriorsDeal() {};
	virtual void calculatePriorsAndPosteriorsBNL() {};
	virtual void setupDefaultInitPriors();
	virtual void setPriorPostNames() {};
	virtual void outputPriors(string & fileName) {};
	virtual void outputPosteriors(string & fileName) {};
	virtual string getNetworkType() { return "?"; };
	virtual string getScoreTypeName() { if(scoreType == 0) return "log likelihood";	else if(scoreType == 1) return "BIC";
			else if(scoreType == 2) return "AIC"; else if(scoreType == 3) return "bayes"; else if(scoreType == 4) return "cost"; else if(scoreType == 5) return "BICprob";	return "?"; };
	virtual string getScoreFixName() { if(scoreFix == 0) return "none"; else if(scoreFix == 1) return "average";
			else if(scoreFix == 2) return "skip"; return "?"; };

	AllNodeData * getAllNodeData() const {return allNodeData;};
	map<unsigned int, DiscreteNode *> getDiscreteNodes() {return discreteNodes;};
	map<unsigned int, CtsNode *> getCtsNodes() {return ctsNodes;};	
	
	//used for setting up the network
	void setUpAllowedEdgeTypes(set<pair<unsigned int, unsigned int> > & allowedEdgeTypesFileNos, set<pair<unsigned int, unsigned int> > & notAllowedEdgeTypesFileNos);
	void addNetworkEdgeType(const unsigned int & fileNo1, const unsigned int & fileNo2);
	void setCostEdgeType(const unsigned int & fileNo1, const unsigned int & fileNo2, const double & cost);
	void setCostEdge(const string & nodeName1, const string & nodeName2, const double & cost);
	double getEdgeCost(const unsigned int & parentNodeNo, const unsigned int & childNodeNo);
	map<pair<unsigned int, unsigned int>, double> getCostEdges() const {return costEdges;};
	map<pair<unsigned int, unsigned int>, double> getCostEdgeTypes() const {return costEdgeTypes;};
	void setCostEdges(map<pair<unsigned int, unsigned int>, double> & costEdges, map<pair<unsigned int, unsigned int>, double> & costEdgeTypes);
	Node * getNetworkNode(unsigned int nodeNo);
	map<unsigned int, Node *> getNetworkNodes() {return allNetworkNodes;};
	Node * getNoDataNetworkNode(const string & nodeName) const;
	bool hasNodesWithNoData() const {return !noDataNetworkNodes.empty();};
	unsigned int getNoNodes() {return allNetworkNodes.size();};
	void addNodeInit(const string & nodeName);
	void addNoDataNode(const string & nodeName, NoData * noDataData);
	virtual void addNode(const unsigned int & nodeNumber) {};	
	void addEdgeInit(const string & nodeName1, const string & nodeName2);
	void removeEdgeInit(const string & nodeName1, const string & nodeName2);
	bool nodeExistsInit(const string & nodeName) const;
	bool edgeExistsInit(const string & nodeName1, const string & nodeName2) const;
	void addWhiteNode(const string & nodeName);
	void addWhiteEdge(const string & nodeName1, const string & nodeName2);
	void addBlackNode(const string & nodeName);
	void addBlackEdge(const string & nodeName1, const string & nodeName2);
	void addNoParentsNode(const string & nodeName);
	void addNoChildrenNode(const string & nodeName);
	bool isNodeCts(const unsigned int & nodeNo) const { return ctsNodes.find(nodeNo) != ctsNodes.end(); };
	bool isNodeDiscrete(const unsigned int & nodeNo) const { return discreteNodes.find(nodeNo) != discreteNodes.end(); };
	void setSNPNodesAsNoParentsNodes();
	bool edgeAllowed(const unsigned int & nodeNo1, const unsigned int & nodeNo2);
	void initialise(); //take away blacklist and add white list
	void setImaginarySampleSize(const double & i) {imaginarySampleSize = i;};
	double getImaginarySampleSize() const {return imaginarySampleSize;};
	NetworkMissingData * updateNetworkMissingData(); //also add ref to nodes
	NetworkMissingData * updateNetworkMissingDataSimNet(const unsigned int & noIndivs);
	pair<unsigned int, unsigned int> getNoMissingNotMissing();
	void setDefaultParameters();
	bool checkNetworkIsValid();
	void clearCache();

	unsigned int getNodesGroupNo();
	void setAllNodesUnfixed();
	void calcEdgeSignifs();
	double getEdgeProb(const string & aParent, const string & nodeStr);
	unsigned int getNumberOfParents(const string & nodeStr);
	unsigned int getNumberOfChildren(const string & nodeStr);
	pair<double, unsigned int> getEdgeSign(const string & aParent, const string & nodeStr);
	unsigned int getNumberOfParameters();
	void setScoreType(const unsigned int & dst) {scoreType = dst;};
	unsigned int getScoreType() {return scoreType; };
	unsigned int getScoreFix() {return scoreFix; };

	//black and white methods
	void updateBlackWhite(Network * network);
	void updateBlackWhiteDifferentData(Network * network, const string & preNodeName);
	set<pair<unsigned int, unsigned int> >  getAllowedEdgeTypes() const {return allowedEdgeTypes;};
	set<unsigned int>  getNoParentsNodes() const {return noParentsNodes;};
	set<unsigned int>  getNoChildrenNodes() const {return noChildrenNodes;};
	set<unsigned int>  getWhiteNodes() const {return whiteNodes;};
	set<unsigned int>  getBlackNodes() const {return blackNodes;};
	set<pair<unsigned int, unsigned int> > getWhiteEdges() const {return whiteEdges;};
	set<pair<unsigned int, unsigned int> > getBlackEdges() const {return blackEdges;};

	//used for modifying network
	void setNetwork(const string & networkStr, const string & preNodeName = "");
	void addEdgesFactorChildren(const unsigned int & nodeNo1, Node * node2);
	void addEdge(Node * node1, Node * node2);
	void addEdge(const unsigned int & nodeNo1, const unsigned int & nodeNo2);
	void addRandomEdge();
	void changeRandomEdge();
	bool isWhiteEdge(const unsigned int & nodeNo1, const unsigned int & nodeNo2);	
	bool isBlackNode(const unsigned int & nodeNo);
	bool nodeExists(const unsigned int & nodeNo) const {map<unsigned int, Node *>::const_iterator nd = allNetworkNodes.find(nodeNo); return nd != allNetworkNodes.end();};
	bool edgeExists(const unsigned int & nodeNo1, const unsigned int & nodeNo2) const;
	void removeNode(const unsigned int & nodeNo);
	void removeEdge(const unsigned int & nodeNo1, const unsigned int & nodeNo2);
	void removeEdge(const unsigned int & nodeNo1, Node * node2);
	void reverseEdge(const unsigned int & nodeNo1, const unsigned int & nodeNo2);
	void reverseEdge(Node * node1, Node * node2);
	void removeAllEdges();
	void deleteAllNetworkNodeData(); //needed for tidy up of nodes after using bootstrap data

	string getNetworkString(const unsigned int & charLimit = 300);
	string getPDAGNetworkString();
	void getPDAGData(const string & PDAGNetwork);
	void calcRecallPrecision(Network * otherNetwork, double & recallOrPrecision, const string & preNodeNameOther, const bool & isRecall);
	void copyMissingness(Network * otherNetwork, const string & preNodeNameOther, const bool & reverse = false);
	void getNumberNodesAndEdges(unsigned int & noDisNodes, unsigned int & noFactorNodes, unsigned int & noCtsNodes, unsigned int & noNoNodes, unsigned int & noEdges);
	NetworkMissingData * getNetworkMissingData();
	void setNetworkMissingData(Network * network);	
	unsigned int getDataSize() const;

	//Markov blanket methods
	set<Node *> getMarkovBlanket(const string & nodeName);
	Network * getMarkovBlanketSubNetwork(const string & nodeName);

	//for calculating all scores
	unsigned int calculateAllScores(const string & allScoresFilename, double & bestNetworkScore, string & bestNetwork); //loop through every poss, return no of networks evaluated
	void calculateAllScoresSetup();	
	void updateNextNetwork(bool & updated);

	//for searches
	virtual void evaluateAllEdgeChanges();
	virtual void evaluateAllEdgeChangesDiffs() {};
#ifdef USING_OPEN_MPI
	virtual bool updateNetworkParallel(double& networkScore) {return false;};
	virtual void evaluateAllEdgeChangesDiffsParallel() {};
#endif
	virtual double cacheAllNodes() {return calcScore();};
	virtual void setScoreBitInCache(const unsigned int & nodeID, const double & value) {};
	virtual double getScoreFromCache() {return calcScore();};
	virtual double getNodeScoreFromCache(const unsigned int & nodeID) {return calcScore();};
	virtual double calcScore2DifferentNodes(const double & origScore, const unsigned int & nodeID1, const unsigned int & nodeID2) {return calcScore();};
	virtual void setScoreBitInCache(const unsigned int & nodeID) {};

	virtual bool updateNetwork(double & networkScore);
	void updateNetworkEdge(unsigned int & parentNo, unsigned int & childNo);
	bool hasLoop();
	void clearVisitedNodes();
	void displayLoop();

	//for imputation
	void setInitialDataIterators2();
	void advanceDataIterators2();
	void setInitialDataIterators3();
	void advanceDataIterators3();
	list<unsigned int> getNodesWithMissingData();
	unsigned int imputeDataNN(const map<unsigned int, set<unsigned int> > & groupNodesWithMissingData, const map<unsigned int, set<unsigned int> > & groupNodesWithCompleteData, Network * bootstrapNetwork, const bool & doNotDoAdjust, const bool & updateImpImmed);
	unsigned int imputeDataMean(const map<unsigned int, set<unsigned int> > & groupNodesWithMissingData, Network * bootstrapNetwork, const bool & updateImpImmed);
	void setStandardDevs();
	list<unsigned int> getConnectedNodes(const unsigned int & nodeID);
	void getNoEdgesNotMissing(const map<unsigned int, set<unsigned int> > & groupNodesWithMissingData, Network * network, unsigned int & noEdgesNotMissing, unsigned int & totalEdges, unsigned int & noNonMissingSingltonNodes, unsigned int & totalSingltonNodes);
	void setNodeMap(const string & preNodeName, Network * bootstrapNetwork);
	unsigned int convertID(const unsigned int & nodeID);
	void setMapID(const unsigned int & thisNetID, const unsigned int & otherNetID);
	void setIsImputedData() {isImputedData = true;};
	bool getIsImputedData() {return isImputedData;};
	void setIsBootstrapData() {isBootstrapData = true;};
	bool getIsBootstrapData() {return isBootstrapData;};
	void randomlyFillMissingValues();
	void setupWasImputed();
	void updateImputedDataAsNonMissing();
	list<pair<string, string> > getNodesSwappable() {return nodesSwappableInit;};

	//for output
	void outputNetwork(string & filename);
	void outputNetworkIGraphEdges(string & filename);
	void outputNetworkIGraphNodes(string & filename);
	void outputNodeData(const string & nodeDataFilePrefix, const bool & outputBedFile, bool & wasDisData, bool & wasCtsData, bool & wasSNPDisData, bool & wasSNPCtsData, bool & wasSNPData, const unsigned int & startIndivNo, const unsigned int & endIndivNo, const unsigned int & jobNo);
	void outputDiscreteNodeData(const string & nodeDataFilePrefix, const bool & outputBedFile, const bool & outSNPData, bool & wasDisData, bool & wasSNPDisData, const unsigned int & startIndivNo, const unsigned int & endIndivNo, const unsigned int & jobNo);
	void outputCtsNodeData(const string & nodeDataFilePrefix, const bool & outputBedFile, const bool & outSNPData, bool & wasCtsData, bool & wasSNPCtsData, const unsigned int & startIndivNo, const unsigned int & endIndivNo, const unsigned int & jobNo);
	void outputSNPNodeData(const string & nodeDataFilePrefix, const bool & outputBedFile, bool & wasSNPData, const unsigned int & startIndivNo, const unsigned int & endIndivNo);
	unsigned int outputEquivalentNetworks(const bool & returnList, list<string> & eqNetworks, const string & nodeNamePre);

};

//! Class for Deal type calculations.
class NetworkDeal : public Network
{
private:

	map<unsigned int, double> discreteJointProbDist; //all discrete nodes group number, prob
	map<unsigned int, string> discreteJointNames;

	//all values will be fixed unless different parents give different marginal prob. dist, we'll assume this

	map<unsigned int, CtsMultiDistJoint *> ctsJointProbDist; //all discrete nodes group number, normal dists, one for each node group
	map<unsigned int, string> meanNames;

	//cts dist paras; used to calc joint prior
	set<DiscreteDealNode *> discreteDealNodes; //only discrete nodes
	set<CtsDealNode *> ctsDealNodes;

public:

	NetworkDeal(AllNodeData * ad, const unsigned int & sc) : Network(ad, sc, 0)
	{

	};

	//! Delete network things.
	~NetworkDeal()
	{
		for (map<unsigned int, CtsMultiDistJoint *>::iterator cjp = ctsJointProbDist.begin(); cjp != ctsJointProbDist.end(); ++cjp) delete cjp->second;
		ctsJointProbDist.clear();		
	};

	//Common network methods
	void calculatePriorsAndPosteriors();
	void addNode(const unsigned int & nodeNumber);
	void outputPriors(string & fileName);
	void outputPosteriors(string & fileName);
	string getNetworkType() { return "deal"; };
	string getScoreTypeName() { return "deal standard"; };
	double calcScore();

	//used for priors, posts and scores
	void calculatePriorsAndPosteriorsDeal() {calculatePriorsAndPosteriors();};
	void clearPriorsAndPosteriors();	
	void setupLocalPriProbs(const string & localProbFilename);
	void calcDiscreteJointPriorDist();
	double calcDiscreteJointProbDist();
	string calcDiscreteJointName();
	double getDiscreteJointPriProb(unsigned int & nodeGroupNo);
	
	double getDiscreteMasterPriorPara();
	void calcDiscreteMasterPriors();
	void calcDiscreteLocalParaPris();
	void calcDiscreteLocalParaPosts();
	CtsMultiDistJoint * getCtsJointProbDist(unsigned int & allNodesGrpNo);
	CtsMultiDistMaster * getCtsMasterPriorDist(CtsNode * ctsNode);
	void calcCtsJointPriorDist();
	void calcCtsMasterPriors();
	void calcCtsLocalParaPris();
	void calcCtsLocalParaPosts();
	
	
	//set prior/post names
	void setPriorPostNames();
	
	void calcDiscreteJointPriorNames();
	void calcCtsJointPriorNames();

	//TEST output
	//void outputPriorsTEST();
};

#ifdef USING_OPEN_MPI
struct BestEval
{
  double eval;
  unsigned int parentID;
  unsigned int childID;
};
#endif

//! Class for bnlearn type calculations and extensions.
class NetworkBNLearn : public Network
{
private:

	set<DiscreteBNLearnNode *> discreteBNLearnNodes; //only discrete nodes
	set<CtsBNLearnNode *> ctsBNLearnNodes;
	list<Node *> orderedNodes; // nodes ordered as in order of dependency
	map<unsigned int, double> nodeScoreCache; //node ID, score bit; for searches

#ifdef USING_OPEN_MPI
	unsigned int currentEvalNumber;	
#endif
	
public:

	NetworkBNLearn(AllNodeData * ad, const unsigned int & sc, const unsigned int & sf) : Network(ad, sc, sf)
	{
		
	};

	//Delete network things
	~NetworkBNLearn()
	{
		
	};

	//Common network methods	
	void calculatePriorsAndPosteriors();
	void calculatePriorsAndPosteriorsBNL() {calculatePriorsAndPosteriors();};
	void addNode(const unsigned int & nodeNumber);
	void outputPriors(string & fileName);
	void outputPosteriors(string & fileName);
	string getNetworkType() { return "bnlearn"; };

	//for searches
	void evaluateAllEdgeChanges();
	void evaluateAllEdgeChangesDiffs();
	double evaluateEdgeChangeDiff(const double & networkScore, const unsigned int & nodeID1, const unsigned int & nodeID2);
	void addEdgeChangeDiffOneNode(const double & networkScore, const unsigned int & nodeID);
	void addEdgeChangeDiffOneEdge(const double & networkScore, const unsigned int & nodeFromID, const unsigned int & nodeToID);
	void evaluateAllEdgeChangesDiffsUpdate(const unsigned int & parentID, const unsigned int & childID, unsigned int & updateType);
#ifdef USING_OPEN_MPI
	bool updateNetworkParallel(double & networkScore);
	void evaluateAllEdgeChangesDiffsParallel();
	void evaluateAllEdgeChangesDiffsUpdateParallel(const unsigned int & parentID, const unsigned int & childID, unsigned int & updateType);
	void addEdgeChangeDiffOneNodeParallel(const double & networkScore, const unsigned int & nodeID);
	void addEdgeChangeDiffOneEdgeParallel(const double & networkScore, const unsigned int & nodeFromID, const unsigned int & nodeToID);
	bool checkLoopForUpdateNetworkEdgeParallel(const unsigned int & parentNo, const unsigned int & childNo, unsigned int & updateType);
	void updateBestEvalParallel(double & bestEval, unsigned int & bestParentID, unsigned int & bestChildID, unsigned int & updateType);
#endif
	double cacheAllNodes();
	void setScoreBitInCache(const unsigned int & nodeID, const double & value);
	void setScoreBitInCache(const unsigned int & nodeID);
	double getNodeScoreFromCache(const unsigned int & nodeID);
	double getScoreFromCache();
	double calcScore2DifferentNodes(const double & origScore, const unsigned int & nodeID1, const unsigned int & nodeID2);
	bool updateNetwork(double & networkScore);
	bool updateNetworkEdgeLoopCheck(const unsigned int & parentNo, const unsigned int & childNo, unsigned int & updateType);

	//set prior/post names
	void setPriorPostNames();

	//bnlearn methods
	void calculateDiscretePosteriors();
	void calculateCtsPosteriors();

	//for simulating data
	void simulateData(unsigned int & noSims, bool & snpDataIntegers);
	void simulateDataOne(bool & snpDataIntegers);
	void setSomeMissingData(const map<string, unsigned int> & missingFirst, map<string, unsigned int> & missingLast, map<string, double> & missingProb);
};

#endif
