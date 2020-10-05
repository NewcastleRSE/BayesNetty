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


/*! \file Nodes.cpp
    \brief This file contains the methods the nodes.
    
*/
#include <iostream>
#include <sstream>
#include <limits>

using namespace std; // initiates the "std" or "standard" namespace

#include "Nodes.h"
#include "main.h"
#include "Task.h"
#include "cdflib.h"
#include "Model.h"
#include "Utils.h"
#include <math.h>

//! Removes parents.
void Node::removeParent(const unsigned int & nodeNo)
{
	map<unsigned int, Node *>::iterator pa = parents.find(nodeNo);
	if(pa != parents.end()) parents.erase(pa);
};

//! Copies parents.
void Node::copyParents(Node * copyingNode, Network * copyingNetwork)
{
	map<unsigned int, Node *>::const_iterator pa = parents.begin();
	Node * parentNode;
	unsigned int nodeID;

	while(pa != parents.end())
	{
		nodeID = pa->second->getNodeID();
		parentNode = copyingNetwork->getNetworkNode(nodeID); //node ids are the same between networks, as id set for the data no

		//copyingNode->addParent(parentNode);		
		copyingNetwork->addEdge(parentNode, copyingNode);

		++pa;
	};

};

//! Copies parents for simulating data.
void Node::simDataCopyParents(const string & simTaskName, Node * copyingNode, Network * copyingNetwork, const bool & createDifferentNodeData)
{
	map<unsigned int, Node *>::const_iterator pa = parents.begin();
	Node * parentNode;
	string parentNodeName;
	AllNodeData * allNodeData = copyingNetwork->getAllNodeData();
	unsigned int nodeID;

	while(pa != parents.end())
	{
		if(createDifferentNodeData) parentNodeName = simTaskName + "-" + pa->second->getName(); //get coressponding node in copying network
		else parentNodeName = pa->second->getName();

		nodeID = allNodeData->getNodeDataNumber(parentNodeName);
		parentNode = copyingNetwork->getNetworkNode(nodeID); //node ids are different between networks (for sim data, bcause data is different), use name modified for sim network

		copyingNetwork->addEdge(parentNode, copyingNode);

		++pa;
	};

};


//! Copies missingness pattern from other node to this one.
void Node::copyMissingness(Node * otherNode)
{
	Data * data = getData();
	Data * otherData = otherNode->getData();

	list<bool>::iterator omv = otherData->missingValues.begin();
	for(list<bool>::const_iterator mv = data->missingValues.begin(); mv != data->missingValues.end(); ++mv, ++omv)
	{
		*omv = *mv;		
	};
};

//! Updates missing data at node level for this node - need to redo for each node each time.
void Node::updateMissingDataAtNodeLevel()
{
	
	unsigned int sizeData = getAmountData(); //use this node for the amount of data
	
	if(freeClearMemory) list<bool>().swap(networkMissingData->missing); else networkMissingData->missing.clear();

	//set up initial missing data as no missing data, then update as missing if any found
	for(unsigned int i = 1; i <= sizeData; ++i) networkMissingData->missing.push_back(false);

	//do for this node first
	updateNetworkMissingData(networkMissingData);

	//update missing for each parent
	for(map<unsigned int, Node *>::iterator pa = parents.begin(); pa != parents.end(); ++pa)
	{		
		pa->second->updateNetworkMissingData(networkMissingData);
	};

	
};

//! Sets node and parent levels to 1.
void Node::setInitialLevels() 
{
	setLevel(1);

	for(map<unsigned int, Node *>::iterator pa = parents.begin(); pa != parents.end(); ++pa)
	{		
		pa->second->setLevel(1);
	};

};

//! Advance node and parent levels to loop thro' all combos.
bool Node::advanceAllLevels()
{
	//advance to next node+parents config
	bool updated = false;
	//try to update this node first
	if(getLevel() < getNoLevels())
	{
		if(isDiscreteNode) setLevel(getLevel() + 1);
		updated = true;
	}
	else
	{
		setLevel(1); //reset this node to 1, beginning level
		//try to update one of the parents now
		for(map<unsigned int, Node *>::iterator pa = parents.begin(); pa != parents.end(); ++pa)
		{
			if(pa->second->getLevel() < pa->second->getNoLevels())
			{
				pa->second->setLevel(pa->second->getLevel() + 1);
				//set prev nodes to 1
				updated = true;
				break;
			}
			else
			{
				pa->second->setLevel(1);
			};
		};
	};

	return updated;
};

//! Advance parent levels to loop thro' all combos.
bool Node::advanceParentLevels()
{
	bool updated = false;
		
	//try to update one of the parents
	for(map<unsigned int, Node *>::iterator pa = parents.begin(); pa != parents.end(); ++pa)
	{				
		if(pa->second->getLevel() < pa->second->getNoLevels())
		{
			pa->second->setLevel(pa->second->getLevel()+1);
			//set prev nodes to 1
			updated = true;
			break;
		}
		else
		{
			pa->second->setLevel(1);
		};				
	};
		
	return updated;
};

//! Returns nodes group number for a certain node level combo.
unsigned int Node::getParentsGroupNo()
{
	if(parents.empty()) return 0;

	unsigned int groupNo = 1;
	unsigned int prevLevels = 1;

	map<unsigned int, Node *>::const_iterator pa = parents.begin();

	while(pa != parents.end())
	{
		groupNo += (pa->second->getLevel() - 1)*prevLevels;
		
		prevLevels *= (pa->second->getNoLevels());

		++pa;
	};

	return groupNo;
};

//! Returns nodes group number for a certain node level combo.
unsigned int Node::getNodeAndParentsGroupNo()
{
	unsigned int groupNo = getLevel();
	unsigned int prevLevels = getNoLevels();

	map<unsigned int, Node *>::const_iterator pa = parents.begin();

	while(pa != parents.end())
	{
		groupNo += (pa->second->getLevel() - 1)*prevLevels;
		
		prevLevels *= (pa->second->getNoLevels());

		++pa;
	};

	return groupNo;
};

//! Returns level from nodes group number.
unsigned int Node::getLevelFromNodeAndParentsGroupNo(const unsigned int & groupNo)
{
	unsigned int noLevels = getNoLevels();
	unsigned int theLevel = groupNo % noLevels;

	return theLevel;
};

//! Returns nodes group number for a certain node level combo for iter2.
unsigned int Node::getNodeAndParentsGroupNo2(Network * network)
{
	unsigned int groupNo = network->getNetworkNode(network->convertID(nodeID))->getNodeDataDisValue2(); 
	unsigned int prevLevels = getNoLevels();

	map<unsigned int, Node *>::const_iterator pa = parents.begin();

	while(pa != parents.end())
	{
		groupNo += (network->getNetworkNode(network->convertID(pa->first))->getNodeDataDisValue2() - 1)*prevLevels;
		
		prevLevels *= (pa->second->getNoLevels());

		++pa;
	};

	return groupNo;
};

//! Returns nodes group number for a certain node level combo for iter3.
unsigned int Node::getNodeAndParentsGroupNo3(Network * network)
{
	unsigned int groupNo = network->getNetworkNode(network->convertID(nodeID))->getNodeDataDisValue3(); 
	unsigned int prevLevels = getNoLevels();

	map<unsigned int, Node *>::const_iterator pa = parents.begin();

	while(pa != parents.end())
	{
		groupNo += (network->getNetworkNode(network->convertID(pa->first))->getNodeDataDisValue3() - 1)*prevLevels;
		
		prevLevels *= (pa->second->getNoLevels());

		++pa;
	};

	return groupNo;
};

//! Returns a random node and parents group, used if unknown when imputing data.
unsigned int Node::getRandomNodeAndParentsGroup2(Network * network)
{
	
	Node * aNode = network->getNetworkNode(network->convertID(nodeID));

	unsigned int groupNo;
	if(!aNode->nodeDataIsMissing2()) groupNo = aNode->getNodeDataDisValue2(); else groupNo = aNode->getRandomLevel(); //should never be missing really

	unsigned int prevLevels = aNode->getNoLevels();
	unsigned int aLevel;
	Node * parentNode;

	map<unsigned int, Node *>::const_iterator pa = parents.begin();

	while(pa != parents.end())
	{

		parentNode = network->getNetworkNode(network->convertID(pa->first));
		
		if(!parentNode->nodeDataIsMissing2()) aLevel = parentNode->getNodeDataDisValue2();
		else
		{
			aLevel = parentNode->getRandomLevel();			
		};

		groupNo += (aLevel - 1)*prevLevels;
		
		prevLevels *= (pa->second->getNoLevels());

		++pa;
	};

	return groupNo;
};

//! Returns a random node and parents group, used if unknown when imputing data.
unsigned int Node::getRandomNodeAndParentsGroup3(Network * network)
{
	Node * aNode = network->getNetworkNode(network->convertID(nodeID));
	unsigned int groupNo;
	if(!aNode->nodeDataIsMissing3()) groupNo = aNode->getNodeDataDisValue3();
	else
	{
		//return 0;
		groupNo = aNode->getRandomLevel();
	};

	unsigned int prevLevels = aNode->getNoLevels();
	unsigned int aLevel;
	Node * parentNode;

	map<unsigned int, Node *>::const_iterator pa = parents.begin();

	while(pa != parents.end())
	{

		parentNode = network->getNetworkNode(network->convertID(pa->first));
		
		if(!parentNode->nodeDataIsMissing3()) aLevel = parentNode->getNodeDataDisValue3();
		else
		{
			//return 0; 
			aLevel = parentNode->getRandomLevel();
		};

		groupNo += (aLevel - 1)*prevLevels;
		
		prevLevels *= (pa->second->getNoLevels());

		++pa;
	};

	return groupNo;
};

//! Output priors for bnlearn cts node for testing.
//void CtsBNLearnNode::outputPriorsTEST()
//{
//	cout << "PRIORS NODE: " << getName() << "\n";
//	for(map<unsigned int, CtsPriLocalDist *>::const_iterator cp = ctsPriLocalDists.begin(); cp != ctsPriLocalDists.end(); ++cp)
//	{
//		cout << " Group no: " << cp->first << "\n";
//		cout << "  Intercept: " << cp->second->intercept << "\n";
//		cout << "  Coefficients: ";
//		for(map<unsigned int, double>::const_iterator cf = cp->second->coeffs.begin(); cf != cp->second->coeffs.end(); ++cf)
//		{
//			cout << cf->first << ": " << cf->second << " ";
//		};
//		cout << "\n";
//		cout << "  Variance: " << cp->second->variance << "\n";
//
//	};
//
//	cout << "\n";
//};

//! Output priors for deal cts node for testing.
//void CtsDealNode::outputPriorsTEST()
//{
//	cout << "PRIORS NODE: " << getName() << "\n";
//	for(map<unsigned int, CtsPriLocalDist *>::const_iterator cp = ctsPriLocalDists.begin(); cp != ctsPriLocalDists.end(); ++cp)
//	{
//		cout << " Group no: " << cp->first << "\n";
//		cout << "  Intercept: " << cp->second->intercept << "\n";
//		cout << "  Coefficients: ";
//		for(map<unsigned int, double>::const_iterator cf = cp->second->coeffs.begin(); cf != cp->second->coeffs.end(); ++cf)
//		{
//			cout << cf->first << ": " << cf->second << " ";
//		};
//		cout << "\n";
//		cout << "  Variance: " << cp->second->variance << "\n";
//	};
//
//	cout << "\n";
//};

//! Return whether the parents have been visited.
bool Node::allParentsVisited()
{
	for(map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end(); ++pa)
	{
		if(!pa->second->getVisited()) return false;
	};

	return true;
};

//! Get name of parents conditional on parents, e.g. [apples|banana,pears].
string Node::getCondName()
{
	string nodeCondName = "[" + getDisplayName();
	bool first = true;
	set<string> parentStrings;

	for(map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end(); ++pa)
	{
		if(!pa->second->getIsFactorChildNode())	parentStrings.insert(pa->second->getDisplayName());
	};

	if(parents.empty()) nodeCondName += "]";
	else
	{
		nodeCondName += "|";

		set<string>::const_iterator pas = parentStrings.begin();

		do{
			if(!first) nodeCondName += ":";
			nodeCondName += *pas;
			first = false;
			
			++pas;
		}while(pas != parentStrings.end());

		nodeCondName += "]";
	};

	return nodeCondName;
};

//! Calculates discrete node+parent name for current discrete node levels.
void Node::calcDiscreteNodeParentNames()
{
	if(freeClearMemory) map<unsigned int, string>().swap(nodeParentNames); else nodeParentNames.clear();

	unsigned int nodeParentsGroup = 1;
	
	setInitialLevels(); 

	//loop thro' node possibilities for node and parent
	do{
		nodeParentsGroup = getNodeAndParentsGroupNo();

		//set the master prior for this node level and parents config
		nodeParentNames[nodeParentsGroup] = getDiscreteParentName(true);

	}while(advanceAllLevels());

};


//! Get name of parents conditional on parents, e.g. [apples|banana,pears].
string Node::getDiscreteParentName(const bool & includeNode)
{
	string name = "";
	if(includeNode && getIsDiscreteNode()) name += getLevelName();
	
	bool doColon = includeNode && getNodeType() == "d";

	for(map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end();)
	{
		if(pa->second->getIsDiscreteNode())
		{
			if(doColon) {name += ":"; doColon = false;};
			name += pa->second->getLevelName();
			++pa;
			if(pa != parents.end() && pa->second->getIsDiscreteNode()) name += ":";
		}
		else
		{
			++pa;
		};
			
	};

	return name;
};

//! Return no missing data in network and total data.
pair<unsigned int, unsigned int> Node::getNoMissingNotMissing()
{
	unsigned int totalData = networkMissingData->missing.size();
	unsigned int missingData = networkMissingData->getNoMissing();

	return make_pair(missingData, totalData-missingData); 
};

//! Get name of parents conditional on parents, e.g. [apples|banana,pears].
void Node::outDiscreteParents(ofstream & outFile, const bool & includeNode)
{
	if(includeNode)
	{
		outFile << " DISCRETE PARENTS: ";
		if(getNodeType() == "d") outFile << getDisplayName();
	}
	else outFile << " DISCRETE PARENTS: ";

	bool doColon = includeNode && getNodeType() == "d";

	for(map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end();)
	{
		if(pa->second->getIsDiscreteNode())
		{
			if(doColon) {outFile << ":"; doColon = false;};
			outFile << pa->second->getDisplayName();
			++pa;
			if(pa != parents.end() && pa->second->getIsDiscreteNode()) outFile << ":";
		}
		else
		{
			++pa;
		};
			
	};
	outFile << "\n";
};

//! Calculates the signifiance of edges given by the parents of the node.
void Node::calcEdgeSignifs(Network * network)
{
	bool isDeal = true;
	NetworkDeal* networkDeal = dynamic_cast<NetworkDeal*>(network);
	if(networkDeal == 0) isDeal = false;

	double initScore;
	if(isDeal) initScore = network->calcScore(); else initScore = calcUpdatedScoreBit(network); //network->getNodeScoreFromCache(nodeID);

	unsigned int initNoParmeters = network->getNumberOfParameters();
	double score2;

	map<unsigned int, Node *> parentsCopy = parents;

	//clear any old calculations from possible previous use
	if(freeClearMemory) map<unsigned int, pair<double, unsigned int> >().swap(edgeSignifs); else edgeSignifs.clear();

	//loop thro' parents and calc signif of each
	for(map<unsigned int, Node *>::const_iterator pa = parentsCopy.begin(); pa != parentsCopy.end(); ++pa)
	{
		//remove edge
		network->removeEdge(pa->first, this);

		if(isDeal) score2 = network->calcScore(); else score2 = calcUpdatedScoreBit(network);

		//record how much better the model is with the edge
		edgeSignifs[pa->first] = make_pair(2*(initScore - score2), initNoParmeters - network->getNumberOfParameters());

		//put back parent
		network->addEdge(pa->second, this);
	};

};

//! Calculates discrete nodes current group level name.
string Node::calcDiscreteParentsName()
{
	string name = "";
	bool added = false;

	//loop thro nodes and set probs if all parents are given, if not skip
	for(map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end(); )
	{
		added = false;
		if(pa->second->getIsDiscreteNode())
		{
			name += pa->second->getLevelName();
			++pa;
			if(pa != parents.end() && pa->second->getIsDiscreteNode()) name += ":";
		}
		else
		{
			++pa;
		};
			
	};

	return name;
};


//! Calculates the local parameter names.
void Node::calcLocalParaNames()
{
	
	if(freeClearMemory) map<unsigned int, string>().swap(localParaNames); else localParaNames.clear();

	unsigned int parentsGroup = 1;
	map<unsigned int, double> localParas;

	//set initial levels
	setInitialLevels(); 
	
	//loop thro' node possibilities for parents
	do{
		parentsGroup = getParentsGroupNo();
	
		localParaNames[parentsGroup] = calcDiscreteParentsName();	
	}while(advanceParentLevels()); //advance to next parents config

};

//! Updates missing data list if missing data in discrete parent data or not in correct group.
void DiscreteNode::updateNotInGroupOrMissing(list<bool> & notInGroupOrMissingData)
{
	if(notInGroupOrMissingData.size() != data->values.size())
	{
		exitErr("Amount of data in a cts node is different to that in a discrete node parent!");
	};

	//set as missing if discrete data is missing or discrete data is not in the level set in the discrete node
	list<bool>::iterator nigomd = notInGroupOrMissingData.begin();
	
	for(list<unsigned int>::const_iterator lv = data->values.begin(); lv != data->values.end(); ++lv, ++nigomd)
	{
		*nigomd = *nigomd || level != *lv;		
	};
};

//! Updates missing data list if missing data in this node. 
void DiscreteNode::updateNetworkMissingData(NetworkMissingData * networkMissingData)
{
	if(networkMissingData->missing.size() != data->missingValues.size())
	{
		exitErr("Amount of data in a discrete node is different to that previously recorded!");
	};

	//set as missing if discrete data is missing or discrete data is not in the level set in the discrete node
	list<bool>::iterator nmd = networkMissingData->missing.begin();
	
	for(list<bool>::const_iterator mv = data->missingValues.begin(); mv != data->missingValues.end(); ++mv, ++nmd)
	{
		*nmd = *nmd || *mv;		
	};
};

//! Updates missing data list if missing data in this node. 
void CtsNode::updateNetworkMissingData(NetworkMissingData * networkMissingData)
{
	if(networkMissingData->missing.size() != data->missingValues.size())
	{		
		exitErr("Amount of data in a continuous node is different to that previously recorded!");
	};

	//set as missing if discrete data is missing or discrete data is not in the level set in the discrete node
	list<bool>::iterator nmd = networkMissingData->missing.begin();

	for(list<bool>::const_iterator mv = data->missingValues.begin(); mv != data->missingValues.end(); ++mv, ++nmd)
	{
		*nmd = *nmd || *mv;		
	};
};

//! Returns the IDs of the cts node parents.
list<unsigned int> CtsNode::getCtsParents()
{
	list<unsigned int> ctsParents;

	for(map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end(); ++pa)
	{
		//add cts data if a cts node
		if(pa->second->getNodeType() == "c") ctsParents.push_back(pa->first);
	};

	return ctsParents;
};

//! Sets node and parental data for a cts node. 
void CtsNode::setParentData(LinearRegModel & linearRegModel)
{
	CovariateData * covariateData = linearRegModel.getCovariateData();

	covariateData->clear(); //clear data for the next time

	covariateData->nodeValues = data; //link node data from cts node to data for regression

	//set initial missing data
	covariateData->notInGroupOrMissingData = networkMissingData->missing; //should be any data missing in network

	//add parent data and update missing or in parent group
	for(map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end(); ++pa)
	{
		pa->second->updateNotInGroupOrMissing(covariateData->notInGroupOrMissingData);

		//add cts data if a cts node
		if(pa->second->getNodeType() == "c") covariateData->covariateDataAllSubjects.push_back(pa->second->getCtsData());
	};
	
};

//! Updates node and parental data for a cts node. 
void CtsNode::updateParentData(LinearRegModel & linearRegModel)
{
	CovariateData * covariateData = linearRegModel.getCovariateData();

	//set initial missing data
	covariateData->notInGroupOrMissingData = networkMissingData->missing;

	//add parent data and update missing or in parent group
	for(map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end(); ++pa)
	{
		pa->second->updateNotInGroupOrMissing(covariateData->notInGroupOrMissingData);
	};

};


//! Sets the discrete data to the initial value.
void CtsNode::setInitialDataValue2()
{
	dataIter2 = data->values.begin(); 
	missIter2 = data->missingValues.begin();
	wasImpIter2 = data->wasImputed.begin();
};

//! Advances data iterators.
void CtsNode::nextDataValueIter2()
{
	dataIter2++;
	missIter2++;
	wasImpIter2++;
};

//! Returns true if node data is missing.
bool CtsNode::nodeDataIsMissing2()
{
	if(*missIter2) return true;

	return false;
};

//! Sets the discrete data to the initial value.
void CtsNode::setInitialDataValue3()
{
	dataIter3 = data->values.begin(); 
	missIter3 = data->missingValues.begin();
};

//! Advances data iterators.
void CtsNode::nextDataValueIter3()
{
	dataIter3++;
	missIter3++;
};

//! Returns true if node data is missing.
bool CtsNode::nodeDataIsMissing3()
{
	if(*missIter3) return true;

	return false;
};

//! Sets estimate of the standard dev
void CtsNode::setStDev()
{
	stDev = data->getStDev();
};

double CtsNode::getMean()
{
	if(!meanSet)
	{
		mean = data->getMean();
		meanSet = true;
	};

	return mean;
};

double CtsNode::getRandomValue()
{
	unsigned int noData = data->getAmountOfData();
	unsigned int randomNo;
	unsigned int i;
	list<bool>::const_iterator m;
	unsigned int tries = 0;

	do{

		randomNo = rand()%200 + 1001; //1801;
		i = 1;
		m = data->missingValues.begin();
		for(list<double>::const_iterator v = data->values.begin(); v != data->values.end(); ++v, ++m)
		{
			if(i == randomNo)
			{
				if(!(*m)) 
				{
					return *v;
				}
				else break;
			};
			++i;
		};

		++tries;
	}while(tries < 1e4);

	out("Failed to get random no for "); out(this->getDisplayName()); out("\n");

	return 0;
};

//! Sets imputed value.
void CtsNode::setImputedDataCts(const double & val, const bool & updateImpImmed)
{
	*dataIter2 = val;
	if(updateImpImmed) *missIter2 = false; //set to non missing later
	*wasImpIter2 = true;
};

//! Returns adjusted data value of node for impuation for imputed indiv or the NN, call with bootstrap Network node.
double CtsNode::getNodeDataCtsValueAdj(const bool & impIndivOrNN, const set<unsigned int> & nodesOfInterest, const list<unsigned int> & bootNodeIDsAdjVars, const list<unsigned int> & adjMethods, Network * network, Network * bootNetwork, const bool & doRandomDisGroup, const unsigned int & grp)
{ 
	//need to update CtsNode::getNodeDataMeanOrStDevAdj to match if this is updated

	unsigned int nodeParentsGroup;
	//impIndivOrNN, true if indiv that is being imputed, false if a possible NN it is being compared to

	if(doRandomDisGroup)
	{
		
		if(impIndivOrNN) nodeParentsGroup = grp; //getRandomNodeAndParentsGroup2(network);
		else
		{
			nodeParentsGroup = getRandomNodeAndParentsGroup3(network);
			if(nodeParentsGroup == 0)
			{				
				return std::numeric_limits<double>::infinity();
			};
		};
	}
	else
	{
		if(impIndivOrNN) nodeParentsGroup = getNodeAndParentsGroupNo2(network); else nodeParentsGroup = getNodeAndParentsGroupNo3(network);  //structure from bootstrap network, data from orig
	};

	CtsPriLocalDist * aCtsPriLocalDist;
	unsigned int aNodeID = network->convertID(getNodeID());
	Node * aNode = network->getNetworkNode(aNodeID);
	map<unsigned int, Node *> theParents = getParents();
	double valueInit, value;
	// get original data for *dataIter2;
	if(impIndivOrNN) valueInit = aNode->getNodeDataCtsValue2(); else valueInit = aNode->getNodeDataCtsValue3();
	set<unsigned int>::const_iterator noi;
	
	map<unsigned int, CtsPriLocalDist *>::const_iterator cpld = ctsPriLocalDists.find(nodeParentsGroup);
	unsigned int bNodeID;
	Node * bNode;	
	
	unsigned int numMissParents = 0;
	bool hasDiscreteParents = false;
	unsigned int missBootID, missBootParentID;

	if(cpld != ctsPriLocalDists.end())
	{
		aCtsPriLocalDist = cpld->second;

		value = valueInit - aCtsPriLocalDist->intercept;
	
		list<unsigned int>::const_iterator am = adjMethods.begin();
		for(list<unsigned int>::const_iterator bv = bootNodeIDsAdjVars.begin(); bv != bootNodeIDsAdjVars.end(); ++bv, ++am)
		{
			bNodeID = network->convertID(*bv); //ensure original data is used		
			bNode = network->getNetworkNode(bNodeID);			
			noi = nodesOfInterest.find(bNodeID);		
			if(noi == nodesOfInterest.end())
			{
				if(!bNode->getIsDiscreteNode())
				{
					if(*am==0)
					{
						if(impIndivOrNN) value -= aCtsPriLocalDist->coeffs[*bv] * bNode->getNodeDataCtsValue2();	//use other variable value as it is not missing
						else value -= aCtsPriLocalDist->coeffs[*bv] * bNode->getNodeDataCtsValue3();
					}
					//else if(*am==1) value -= 0;	//no adjust
					else if(*am==2) value -= aCtsPriLocalDist->coeffs[*bv] * bNode->getMean();	//use other variable mean as it is missing
					else if(*am==3) value -= aCtsPriLocalDist->coeffs[*bv] * bNode->getRandomValue();	//use other variable mean as it is missing
				}
				else
				{
					hasDiscreteParents = true;
				};
			};
			
		};

		
		if(hasDiscreteParents)
		{
			numMissParents = 0;
			for(set<unsigned int>::const_iterator ni = nodesOfInterest.begin(); ni != nodesOfInterest.end(); ++ni)
			{
				//check if a parent of this node in the bootstrap network	
				missBootID = parentExists(bootNetwork->convertID(*ni));
				if(parentExists(missBootID))
				{
					missBootParentID = missBootID;
					numMissParents++;	
				};
			};

			//estimate C as C = (B - beta_0^i)/beta_C^i, accounting for discrete parents by dividing by coeff for C which may be different. [A]->(B)<-(C), C is missing, A discrete, B and C cts
			if(numMissParents == 1) value /= aCtsPriLocalDist->coeffs[missBootParentID]; 
			
		};

		
	}
	else
	{	
		if(impIndivOrNN) exitErr("Problem adjusting continuous variables when imputing data (value 2)!");
		else exitErr("Problem adjusting continuous variables when imputing data (value 3)!");	
	};

	return value;
};

//! Replaces missing values with randomly sampled values that are not missing.
void CtsNode::randomlyFillMissingValues()
{
	list<bool>::iterator m = data->missingValues.begin();
	for(list<double>::iterator v = data->values.begin(); v != data->values.end(); ++v, ++m)
	{
		if(*m)
		{
			*v = data->getRandomCompleteValue();
			*m = false;
		};
	};
};

CtsPriLocalDist * CtsNode::getCtsPriLocalDist(const unsigned int & groupNodeID)
{
	map<unsigned int, CtsPriLocalDist *>::const_iterator cp = ctsPriLocalDists.find(groupNodeID);

	if(cp == ctsPriLocalDists.end())
	{
		string msg = "Problem retrieving coefficients for " + getDisplayName() + "!\n";
		exitErr( msg);
	};

	return cp->second;
};

//! Returns mean for adjusted cts data, called from original data network.
double CtsNode::getNodeDataMeanOrStDevAdj(const set<unsigned int> & nodesOfInterest, const list<unsigned int> & bootNodeIDsAdjVars, const list<unsigned int> & adjMethods, Network * network, Network * bootstrapNetwork, const bool & doStDev, const double & mean)
{
	unsigned int nodeParentsGroup;
	CtsPriLocalDist * aCtsPriLocalDist;
	double value = 0;
	
	map<unsigned int, Node *> theParents;
	list<unsigned int>::const_iterator am; //adjust methods

	//cts parent iterators
	list< list<double>::const_iterator > parCtsIts;
	list< list<double>::const_iterator >::const_iterator parCtsItsIt;
	list<double>::const_iterator aParCtsIt;

	//discrete parent iterators
	list< list<unsigned int>::const_iterator > parDisIts;
	list< list<unsigned int>::const_iterator >::const_iterator parDisItsIt;
	list<unsigned int>::const_iterator aParDisIt;
	set<unsigned int>::const_iterator noi;
	
	unsigned int aNodeID;

	//convert to data for this network
	for(list<unsigned int>::const_iterator bv = bootNodeIDsAdjVars.begin(); bv != bootNodeIDsAdjVars.end(); ++bv)
	{
		aNodeID = network->convertID(*bv);
		theParents[aNodeID] = network->getNetworkNode(aNodeID);
	};

	for(map<unsigned int, Node *>::const_iterator pa = theParents.begin(); pa != theParents.end(); ++pa)
	{
		if(pa->second->getNodeType() == "c") parCtsIts.push_back(pa->second->getCtsData()->values.begin());
		else if(pa->second->getNodeType() == "d") parDisIts.push_back(pa->second->getDiscreteData()->values.begin()); 
	};

	double total = 0;
	double count = 0;
	double diff;
	
	unsigned int aBootNodeID = bootstrapNetwork->convertID(getNodeID());
	Node * aBootNode = bootstrapNetwork->getNetworkNode(aBootNodeID); //corresponding
	
	unsigned int paBootNodeID;
	Node * paBootNode;
	Node * bNode;
	Node * bBootNode;

	map<unsigned int, list<unsigned int>::const_iterator > noiLevIts; // bootstrap node ID, data values iterators - for discrete nodes 

	//update values in bootstrap network for discrete nodes so that thr group can be calculated
	for(set<unsigned int>::const_iterator noi = nodesOfInterest.begin(); noi != nodesOfInterest.end(); ++noi)
	{
		bNode = network->getNetworkNode(*noi);
		if(bNode->getIsDiscreteNode())
		{			
			noiLevIts[bootstrapNetwork->convertID(*noi)] = bNode->getDiscreteData()->values.begin();
		};

	};

	list<bool>::const_iterator mv = network->getNetworkMissingData()->missing.begin();

	unsigned int numMissParents = 0;
	bool hasDiscreteParents = false;
	unsigned int missBootID, missBootParentID;
	double valueInit;

	for(list<double>::const_iterator v = data->values.begin(); v != data->values.end(); ++v, ++mv)
	{
		//set parent data... loop thro' cts and dis parents and set values
		parCtsItsIt = parCtsIts.begin();
		parDisItsIt = parDisIts.begin();

		for(map<unsigned int, Node *>::const_iterator pa = theParents.begin(); pa != theParents.end(); ++pa)
		{
			if(pa->second->getNodeType() == "c")
			{
				aParCtsIt = *parCtsItsIt;
				pa->second->setValue(*aParCtsIt);
				parCtsItsIt++;
			}
			else if(pa->second->getNodeType() == "d")
			{
				aParDisIt = *parDisItsIt;				
				pa->second->setLevel(*aParDisIt);
				parDisItsIt++;

				//set discrete value in bootstrap network node, may have different groups if discrete parents
				paBootNodeID = bootstrapNetwork->convertID(pa->first);
				paBootNode = bootstrapNetwork->getNetworkNode(paBootNodeID);
				paBootNode->setLevel(*aParDisIt);
			};
		};

		if(!*mv)
		{	
			
			for(map<unsigned int, list<unsigned int>::const_iterator >::iterator noii = noiLevIts.begin(); noii != noiLevIts.end(); ++noii)
			{
				bBootNode = bootstrapNetwork->getNetworkNode(noii->first);
				bBootNode->setLevel(*(noii->second));

			};
			

			nodeParentsGroup = aBootNode->getNodeAndParentsGroupNo();

			aCtsPriLocalDist = aBootNode->getCtsPriLocalDist(nodeParentsGroup);

			valueInit = *v;

			if(aCtsPriLocalDist != 0)
			{
				
				value = valueInit - aCtsPriLocalDist->intercept;
				hasDiscreteParents = false;

		
				am = adjMethods.begin();
			
				for(map<unsigned int, Node *>::const_iterator pa = theParents.begin(); pa != theParents.end(); ++pa)
				{							
					noi = nodesOfInterest.find(pa->first);		
					if(noi == nodesOfInterest.end())
					{
						if(!pa->second->getIsDiscreteNode())
						{
							if(*am==0) value -= aCtsPriLocalDist->coeffs[bootstrapNetwork->convertID(pa->first)] * pa->second->getValue();	//use other variable value as it is not missing
							//else if(*am==1) value -= 0; //no adjust
							else if(*am==2) value -= aCtsPriLocalDist->coeffs[bootstrapNetwork->convertID(pa->first)] * pa->second->getMean();	//use other variable mean as it is missing
							else if(*am==3) value -= aCtsPriLocalDist->coeffs[bootstrapNetwork->convertID(pa->first)] * pa->second->getRandomValue();	//use other variable random value as it is missing
						}
						else
						{
							hasDiscreteParents = true;
						};
					};
			
					
				};

				if(hasDiscreteParents)
				{

					numMissParents = 0;
					for(set<unsigned int>::const_iterator ni = nodesOfInterest.begin(); ni != nodesOfInterest.end(); ++ni)
					{
						//check if a parent of this node in the bootstrap network	
						missBootID = parentExists(bootstrapNetwork->convertID(*ni));
						if(parentExists(missBootID))
						{
							missBootParentID = missBootID;
							numMissParents++;	
						};
					};

					//estimate C as C = (B - beta_0^i)/beta_C^i, accounting for discrete parents by dividing by coeff for C which may be different. [A]->(B)<-(C), C is missing, A discrete, B and C cts
					if(numMissParents == 1) value /= aCtsPriLocalDist->coeffs[missBootParentID]; 
					
				};

				
			}
			else
			{
				exitErr("Problem adjusting continuous variables when imputing data (mean)!");
			};

			if(doStDev)
			{
				diff = value - mean;
				total += diff*diff;
			}
			else total += value;

			count++;
		};


		//advance iterators for cts parent data
		for(list< list<double>::const_iterator >::iterator paic = parCtsIts.begin(); paic != parCtsIts.end(); ++paic)
		{
			++(*paic);
		};

		//advance iterators for dis parent data
		for(list< list<unsigned int>::const_iterator >::iterator paid = parDisIts.begin(); paid != parDisIts.end(); ++paid)
		{
			++(*paid);
		};

		//advance iterators for dis nodes of interest
		for(map<unsigned int, list<unsigned int>::const_iterator >::iterator noii = noiLevIts.begin(); noii != noiLevIts.end(); ++noii)
		{
			++(noii->second);
		};
	};

	if(!doStDev) return total/count;

	return sqrt(total/(count-1));
};

//! Returns st dev for adjusted cts data, called from original data network.
double CtsNode::getNodeDataStDevAdj(const set<unsigned int> & nodesOfInterest, const list<unsigned int> & bootNodeIDsAdjVars, const list<unsigned int> & adjMethods, Network * network, Network * bootstrapNetwork)
{
	double mean = getNodeDataMeanOrStDevAdj(nodesOfInterest, bootNodeIDsAdjVars, adjMethods, network, bootstrapNetwork, false, 0);

	return getNodeDataMeanOrStDevAdj(nodesOfInterest, bootNodeIDsAdjVars, adjMethods, network, bootstrapNetwork, true, mean);
};

//! Delete and clear prior and posterior data.
void CtsDealNode::clearPriorsAndPosteriors()
{
	clearLocalPriors();
	for(map<unsigned int, CtsMultiDistMaster *>::iterator cmd = ctsMasterPriDists.begin(); cmd != ctsMasterPriDists.end(); ++cmd) delete cmd->second;
	for(map<unsigned int, CtsLocalPara *>::iterator ppri = ctsLocalParaPris.begin(); ppri != ctsLocalParaPris.end(); ++ppri) delete ppri->second;
	for(map<unsigned int, CtsLocalPara *>::iterator ppost = ctsLocalParaPosts.begin(); ppost != ctsLocalParaPosts.end(); ++ppost) delete ppost->second;

	if(freeClearMemory)
	{
		map<unsigned int, CtsMultiDistMaster*>().swap(ctsMasterPriDists); 
		map<unsigned int, CtsLocalPara* >().swap(ctsLocalParaPris);  
		map<unsigned int, CtsLocalPara* >().swap(ctsLocalParaPosts);  
		list<unsigned int>().swap(nodeParentIDsForData);
		
		map<unsigned int, unsigned int>().swap(totalInData);
		map<unsigned int, list< list<double> > >().swap(parentDataMat);
		map<unsigned int, list< list<double> > >().swap(parentDataMatTrans);
		map<unsigned int, list<double> >().swap(nodeData);
	}
	else
	{
		ctsMasterPriDists.clear();
		ctsLocalParaPris.clear();
		ctsLocalParaPosts.clear();
		nodeParentIDsForData.clear();

		//clear data used to calc score
		totalInData.clear();
		parentDataMat.clear();
		parentDataMatTrans.clear();
		nodeData.clear();
	};

};

//! Delete and clear prior and posterior data.
void DiscreteDealNode::clearPriorsAndPosteriors()
{
	if(freeClearMemory)
	{
		map<unsigned int, double>().swap(discreteMasterParaPri); 
		map<unsigned int, map<unsigned int, double> >().swap(discreteLocalParaPris);
		map<unsigned int, map<unsigned int, double> >().swap(discreteLocalParaPosts);
	}
	else
	{
		discreteMasterParaPri.clear();
		discreteLocalParaPris.clear();
		discreteLocalParaPosts.clear();
	};
};


//! Sets up local prior dist as normal dist using linear regression on parents.
void CtsNode::setDefaultInitPrior()
{
	//look in cache first
	bool linearRegDone = false;

	//Some problem with keeping a cache? Too big...?
	/*string parentsName = "";
	for(map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end(); )
	{
		parentsName += toString(pa->first);
		++pa;
		if(pa != parents.end()) parentsName += ":";
	};

	map<string, map<unsigned int, CtsPriLocalDist *> >::const_iterator cpld = ctsPriLocalDistsCache.find(parentsName);
	
	if(cpld != ctsPriLocalDistsCache.end())
	{
		ctsPriLocalDists = cpld->second;
		linearRegDone = true;
	};
	*/

	//loop thro' each config of the discrete parents
	unsigned int nodeParentsGroup = 1;
	CtsPriLocalDist * aCtsPriLocalDist;
	LinearRegModel linearRegModel;
	CovariateData * covariateData = linearRegModel.getCovariateData();

	setInitialLevels(); 

	bool updated = false;
	double rss;

	nodeParentsGroup = getNodeAndParentsGroupNo();

	//set up nodeParentIDsForData, list of node+parent IDs corresponding to the group 
	for(list<bool>::const_iterator mv = networkMissingData->missing.begin(); mv != networkMissingData->missing.end(); ++mv)
	{
		if(*mv) nodeParentIDsForData.push_back(0); //missing data
		else nodeParentIDsForData.push_back(nodeParentsGroup);
	};

	//setup model data with subset of data for the set of discrete parents	
	setParentData(linearRegModel);

	list<unsigned int>::iterator npi;

	//loop thro' node possibilities for node and parent
	do{
		
		if(!linearRegDone) 
		{
			aCtsPriLocalDist = new CtsPriLocalDist();
			ctsPriLocalDistsToDel.push_back(aCtsPriLocalDist);

			//do linear regression
			linearRegModel.fitModelCovar(rss);
			aCtsPriLocalDist->intercept = linearRegModel.getParameter(1);
			unsigned int paraNo = 2;
			for(map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end(); ++pa)
			{			
				//add cts node para if a cts node
				if(pa->second->getNodeType() == "c")
				{
					aCtsPriLocalDist->coeffs[pa->first] = linearRegModel.getParameter(paraNo);
					paraNo++;
				};
			};

			aCtsPriLocalDist->variance = linearRegModel.getVariance();
			aCtsPriLocalDist->mean = linearRegModel.getMean();
			aCtsPriLocalDist->stDev = sqrt(aCtsPriLocalDist->variance);

			//save results of linear regression
			ctsPriLocalDists[nodeParentsGroup] = aCtsPriLocalDist;

		}; //end of linear reg

	    if(!advanceParentLevels()) break; //advance to next parents config break;

		nodeParentsGroup = getNodeAndParentsGroupNo();

		//setup model data with subset of data for the set of discrete parents
		updateParentData(linearRegModel);

		//update list of node+parent IDs for posterior calculation later
		npi = nodeParentIDsForData.begin();
		for(list<bool>::const_iterator cd = covariateData->notInGroupOrMissingData.begin(); cd != covariateData->notInGroupOrMissingData.end(); ++cd, ++npi)
		{
			if(!*cd) *npi = nodeParentsGroup;			
		};

	}while(true); //end of loop for node poss'es

};

//! Sets up local prior dist as uniform dists for discrete nodes.
void DiscreteNode::setDefaultInitPrior()
{
	uniformProb = 1.0/(double)(this->getNoLevels()); 
};

//! Sets priJointBitProb if parents (if any) have defined probs for calc'ing joint prior.
double DiscreteNode::getPriLocalProb()
{
	if(uniformProb >= 0) return uniformProb;

	return uniformProb; //no other options set up
};

//! Sets priJointBitProb if parents (if any) have defined probs for calc'ing joint prior.
double DiscreteDealNode::getMasterPriPara()
{
	//get pri local prob given parents
	unsigned int nodeAndParentGrpNo = getNodeAndParentsGroupNo();
	map<unsigned int, double>::const_iterator dmpp = discreteMasterParaPri.find(nodeAndParentGrpNo);

	if(dmpp == discreteMasterParaPri.end())
	{
		string mess = "Master parameter prior not set up for node " + this->getDisplayName() + "!";
		exitErr(mess);
	};

	return dmpp->second;
};

//! Updates the cts joint prob dist for this node using prev calculations.
void CtsDealNode::updateCtsJointProbDist(CtsMultiDistJoint * ctsJointDist)
{
	double mean, variance, covariance;
	map<unsigned int, double> covariances;

	map<unsigned int, double>::const_iterator paMe; //find parent mean
	map<unsigned int, double>::const_iterator paVa; //find parent variance

	unsigned int nodeParentsGrpNo = getNodeAndParentsGroupNo();
	map<unsigned int, CtsPriLocalDist *>::const_iterator cpl = ctsPriLocalDists.find(nodeParentsGrpNo);
	if(cpl == ctsPriLocalDists.end())
	{
		out("Continuous prior local distribution not found for node "); out(getDisplayName()); out("!\n");
		exit(1);
	};

	//calculate mean, variance and covariances for the cts joint prob dist for this cts node
	mean = cpl->second->intercept;
	variance = cpl->second->variance;
	
	for(map<unsigned int, double>::const_iterator co = cpl->second->coeffs.begin(); co != cpl->second->coeffs.end(); ++co)
	{

		paMe = ctsJointDist->means.find(co->first);
		if(paMe == ctsJointDist->means.end())
		{
			out("Mean not found in continuous joint prior distribution for a parent node of node: "); out(getDisplayName()); out("!\n");
			exit(1);
		};
		
		mean += (co->second)*(paMe->second);

		covariance = 0;
		for(map<unsigned int, double>::const_iterator co2 = cpl->second->coeffs.begin(); co2 != cpl->second->coeffs.end(); ++co2)
		{			
			covariance += (co2->second)*ctsJointDist->getCovariance(co->first, co2->first);
		};

		variance += (co->second)*covariance;
		covariances[co->first] = covariance;
	};

	ctsJointDist->means[nodeID] = mean;
	ctsJointDist->variances[nodeID] = variance;
	ctsJointDist->covariances[nodeID] = covariances;
};


//! Calculates the master prior for a node.
void DiscreteDealNode::calcDiscreteMasterPrior(NetworkDeal * network)
{
	unsigned int nodeParentsGroup = 1;

	//fix nodes for this node and parents, and initial level, others not fixed.
	network->setAllNodesUnfixed();
	setFixed(true); //fix this node
	level = 1;
	for(map<unsigned int, Node *>::iterator pa = parents.begin(); pa != parents.end(); ++pa)
	{
		pa->second->setFixed(true);
		pa->second->setLevel(1);
	};

	//loop thro' node possibilities for node and parent
	do{
		nodeParentsGroup = getNodeAndParentsGroupNo();
	
		//set the master prior for this node level and parents config
		discreteMasterParaPri[nodeParentsGroup] = network->getDiscreteMasterPriorPara();

	}while(advanceAllLevels()); //advance to next node+parents config

};


//! Calculates master prior for a cts node.
void CtsDealNode::calcCtsMasterPrior(NetworkDeal * network)
{
	//loop thro' the discrete node configs

	unsigned int nodeParentsGroup = 1;

	//fix nodes for this node and parents, and initial level, others not fixed.
	network->setAllNodesUnfixed();
	
	for(map<unsigned int, Node *>::iterator pa = parents.begin(); pa != parents.end(); ++pa)
	{
		pa->second->setFixed(true);
		pa->second->setLevel(1);
	};

	//loop thro' node possibilities for node and parent
	do{
		nodeParentsGroup = getNodeAndParentsGroupNo();

		//set the master prior for this node level and parents config
		ctsMasterPriDists[nodeParentsGroup] = network->getCtsMasterPriorDist(this);

	}while(advanceParentLevels());
};

//! Calculates the local parameter prior for a cts node.
void CtsDealNode::calcCtsLocalParameterPrior(NetworkDeal * network)
{
	unsigned int nodeParentsGroup = 1;

	//fix nodes for this node and parents, and initial level, others not fixed.
	network->setAllNodesUnfixed();
	
	for(map<unsigned int, Node *>::iterator pa = parents.begin(); pa != parents.end(); ++pa)
	{
		pa->second->setFixed(true);
		pa->second->setLevel(1);
	};

	//loop thro' node possibilities for node and parent
	do{
		nodeParentsGroup = getNodeAndParentsGroupNo();

		//set the local parameter prior for this node level and parents config	
		ctsLocalParaPris[nodeParentsGroup] = getCtsLocalParaPrior(nodeParentsGroup);
	
	}while(advanceParentLevels()); //advance to next parent config

};

//! Calculates the local parameter prior for a cts node and a given node+parents group.
CtsLocalPara * CtsDealNode::getCtsLocalParaPrior(unsigned int & nodeParentsGroup)
{
	CtsLocalPara * ctsLocalPara = new CtsLocalPara();

	//decompose covariance matrix
	double phi1by1;
	list<double> phiVec;
	list<list<double> > phiMatpbyp;
	list<list<double> > phiMatpbypInverse;
	CtsMultiDistMaster * ctsMultiDistMaster;
	double mu1by1;
	list<double> muVecpby1;

	map<unsigned int, CtsMultiDistMaster *>::const_iterator cm = ctsMasterPriDists.find(nodeParentsGroup);

	if(cm != ctsMasterPriDists.end()) ctsMultiDistMaster = cm->second;
	else {
		exitErr("Cannot find master prior of a continuous node!");
		exit(0); //only to avoid build error... as exit in above function
	};

	list<double> firstRow = *ctsMultiDistMaster->covariances.begin();

	if(firstRow.empty()) exitErr("Problem with empty master prior covariance matrix!");

	phi1by1 = *firstRow.begin();
	firstRow.pop_front();
	phiVec = firstRow;

	phiMatpbyp = ctsMultiDistMaster->covariances;
	phiMatpbyp.pop_front(); //knock off first row

	for(list<list<double> >::iterator r = phiMatpbyp.begin(); r != phiMatpbyp.end(); ++r)
	{		
		r->pop_front(); //knock off first column
	};

	muVecpby1 = ctsMultiDistMaster->means;
	mu1by1 = *muVecpby1.begin();
	muVecpby1.pop_front();

	//calcualte inverse
	unsigned int matSize = phiMatpbyp.size();

	//set inverse matrix as identity first
	list<double> aRow;
	for(unsigned int ri = 1; ri <= matSize; ++ri)
	{
		aRow.clear();
		for(unsigned int ci = 1; ci <= matSize; ++ci)
		{
			if(ri != ci) aRow.push_back(0);
			else aRow.push_back(1);
		};
		phiMatpbypInverse.push_back(aRow);
	};

	getInverseMatrix(phiMatpbyp, phiMatpbypInverse);

	//calcualte mu
	double newmu1;
	list<double> phiInvbyMupa;
	double newmu1bit2;

	list<double> phibyPhiInv;
	getVecMatrixMulti(phiVec, phiMatpbypInverse, phibyPhiInv);
	getVecVecMulti(phibyPhiInv, muVecpby1, newmu1bit2);

	ctsLocalPara->mu = phibyPhiInv;

	newmu1 = mu1by1 - newmu1bit2; 
	ctsLocalPara->mu.push_front(newmu1);

	//calculate phi
	double phiBit2;
	getVecVecMulti(phibyPhiInv, phiVec, phiBit2);

	double newPhi = phi1by1 - phiBit2;
	ctsLocalPara->phi = newPhi;

	//calculate rho
	ctsLocalPara->rho = ctsMultiDistMaster->jointDiscretePriorParaSum + phiVec.size(); 

	//calculate tau
	list<double> tauFirstRow;

	if(muVecpby1.empty()) //simple for 1 dimension tau
	{
		tauFirstRow.push_back(ctsMultiDistMaster->jointDiscretePriorParaSum);
		ctsLocalPara->tau.push_back(tauFirstRow);

		//set inverse
		tauFirstRow.push_back(1.0/(*tauFirstRow.begin()));
		tauFirstRow.pop_front();
		ctsLocalPara->tauInverse.push_back(tauFirstRow);
	}
	else
	{
		list<list<double> > tauMatInverse, tauMat; //setup inverse tau first and inverse tau inverse to get tau
		double tau1by1, tauBit2;
		list<double> tauMatRow;//set up tau as identity for use with inverse method
		list<double> phiInvmu;
		
		getMatrixVecMulti(phiMatpbypInverse, muVecpby1, phiInvmu);
		getVecVecMulti(muVecpby1, phiInvmu, tauBit2);

		tau1by1 = 1.0/ctsMultiDistMaster->jointDiscretePriorParaSum + tauBit2;
	
		tauFirstRow.push_back(tau1by1);
	
		for(list<double>::const_iterator mp = phiInvmu.begin(); mp != phiInvmu.end(); ++mp)
		{
			tauFirstRow.push_back(-(*mp));
			
		};

		tauMatInverse.push_back(tauFirstRow);
		
		list<double>::const_iterator i = phiInvmu.begin();		
		for(list<list<double> >::const_iterator piv = phiMatpbypInverse.begin(); piv != phiMatpbypInverse.end(); ++piv, ++i)
		{
			aRow = *piv;
			aRow.push_front(-(*i));

			tauMatInverse.push_back(aRow);
		};

		//setting up tauMat as identity
		for(unsigned int rowNo = 1; rowNo <= tauMatInverse.size(); ++rowNo)
		{
			tauMatRow.clear();
			for(unsigned int colNo = 1; colNo <= tauMatInverse.size(); ++colNo)
			{
				if(rowNo == colNo) tauMatRow.push_back(1); else tauMatRow.push_back(0);
			};
			tauMat.push_back(tauMatRow); 
		};

		//save for later, before it is turned to the identity by following function	
		copyMatrix(tauMatInverse, ctsLocalPara->tauInverse);
		getInverseMatrix(tauMatInverse, tauMat); //inverse "inverse tau" to get tau
		copyMatrix(tauMat, ctsLocalPara->tau);
	};

	return ctsLocalPara;
};


//! Calculates the local parameter prior for a cts node.
void CtsDealNode::setupNodeAndParentData()
{
	if(freeClearMemory)
	{
		map<unsigned int, unsigned int>().swap(totalInData);
		map<unsigned int, list< list<double> > >().swap(parentDataMat);
		map<unsigned int, list< list<double> > >().swap(parentDataMatTrans);
		map<unsigned int, list<double> >().swap(nodeData);
	}
	else
	{
		totalInData.clear();
		parentDataMat.clear();
		parentDataMatTrans.clear();
		nodeData.clear();
	};

	map<unsigned int, unsigned int>::iterator tot;
	map<unsigned int, list< list<double> > >::iterator par;
	map<unsigned int, list<double> >::iterator nod;

	//parent iterators
	list< list<double>::const_iterator > parIts;

	for(map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end(); ++pa)
	{
		if(pa->second->getNodeType() == "c") parIts.push_back(pa->second->getCtsData()->values.begin());
	};

	list<double> aList;
	list< list<double> > aMat;
	list< list<double> >::iterator rowParMat;

	//initialize to empty/zero counts
	for(map<unsigned int, CtsPriLocalDist *>::const_iterator cpld = ctsPriLocalDists.begin(); cpld != ctsPriLocalDists.end(); ++cpld)
	{
		totalInData[cpld->first] = 0;

		aList.clear();	
		nodeData[cpld->first] = aList;

		aMat.clear();
		parentDataMat[cpld->first] = aMat;

		parentDataMatTrans[cpld->first] = aMat;
	};

	//loop thro' cts data, and list saying which node+parent group or missing
	list<unsigned int>::const_iterator npi = nodeParentIDsForData.begin();

	for(list<double>::const_iterator cd = data->values.begin(); cd != data->values.end(); ++cd, ++npi)
	{

		if(*npi != 0)
		{
			//do total, n
			tot = totalInData.find(*npi);
			if(tot == totalInData.end())
			{
				totalInData[*npi] = 1;
			}
			else
				totalInData[*npi] += 1;

			//do node data, y
			nod = nodeData.find(*npi);
			if(nod == nodeData.end())
			{
				aList.clear();
				aList.push_back(*cd);
				nodeData[*npi] = aList;
			}
			else
			{
				nod->second.push_back(*cd);
			};

			//do parent data, z
			aList.clear();

			//add row to z matrix
			aList.push_back(1);
			for(list< list<double>::const_iterator >::iterator pai = parIts.begin(); pai != parIts.end(); ++pai)
			{
				aList.push_back(*(*pai));
			};

			par = parentDataMat.find(*npi);
			if(par == parentDataMat.end())
			{				
				aMat.clear();				
				aMat.push_back(aList);
				parentDataMat[*npi] = aMat;
			}
			else
			{				
				par->second.push_back(aList);
			};

			//do parent data, z
			par = parentDataMatTrans.find(*npi);
			if(par == parentDataMatTrans.end() || par->second.empty())
			{
				aMat.clear(); aList.clear(); aList.push_back(1);				
				aMat.push_back(aList);
				for(list< list<double>::const_iterator >::iterator pai = parIts.begin(); pai != parIts.end(); ++pai)
				{
					aList.clear(); aList.push_back(*(*pai));
					aMat.push_back(aList);
				};

				parentDataMatTrans[*npi] = aMat;
			} else {
				rowParMat = par->second.begin();
				rowParMat->push_back(1);
				for(list< list<double>::const_iterator >::iterator pai = parIts.begin(); pai != parIts.end(); ++pai)
				{
					++rowParMat;
					rowParMat->push_back(*(*pai));
				};
			};

		};//end, if missing data

		//advance iterators for parent data
		for(list< list<double>::const_iterator >::iterator pai = parIts.begin(); pai != parIts.end(); ++pai)
		{
			++(*pai);
		};

	};//end of cts data loop

};

//! Calculates the local parameter prior for a cts node.
void CtsDealNode::calcCtsLocalParameterPost(NetworkDeal * network)
{
	if(totalInData.empty()) setupNodeAndParentData();

	CtsLocalPara * posterior;

	map<unsigned int, unsigned int>::iterator tot = totalInData.begin();
	map<unsigned int, list< list<double> > >::iterator par = parentDataMatTrans.begin();
	map<unsigned int, list<double> >::iterator nod = nodeData.begin();

	list< list<double> > zTz;
	list< list<double> > tauCopy, inverseTauPost, identity;
	list<double> aRow;
	unsigned int tauSize = ctsLocalParaPris.begin()->second->tau.size();

	for(unsigned int rowNo = 1; rowNo <= tauSize; ++rowNo)
	{
		aRow.clear();
		for(unsigned int colNo = 1; colNo <= tauSize; ++colNo)
		{
			if(rowNo == colNo) aRow.push_back(1); else aRow.push_back(0);
		};
		identity.push_back(aRow); 
	};
	
	list<double> vecSum, vecZMu, vecDiff1, vecDiffMu, vecTauMu, vecZY;
	double phiBit1, phiBit2;
	
	if(totalInData.size() != ctsLocalParaPris.size() || parentDataMatTrans.size() != ctsLocalParaPris.size() || nodeData.size() != ctsLocalParaPris.size())
	{
		exitErr("Error: Problem calculating posteriors!");
	};

	//create new posteriors
	for(map<unsigned int, CtsLocalPara * >::const_iterator clpp = ctsLocalParaPris.begin(); clpp != ctsLocalParaPris.end(); ++clpp, ++par, ++nod, ++tot)
	{
		//create post..
		posterior = new CtsLocalPara();
		ctsLocalParaPosts[clpp->first] = posterior;

		//posterior tau
		if(!par->second.empty()) //set posterior to prior if no info
		{
			getMatrixTransMatrixMulti(par->second, zTz); 
			getMatrixAddMatrix(clpp->second->tau, zTz, posterior->tau);
		}
		else
		{
			posterior->tau = clpp->second->tau;
		};

		//posterior mu
		if(!par->second.empty())  //set posterior to prior if no info
		{
			inverseTauPost = identity;
			copyMatrix(posterior->tau, tauCopy); //preserve the tau values as the inverse method would change it to the identity!
			getInverseMatrix(tauCopy, inverseTauPost);
			getMatrixVecMulti(clpp->second->tau, clpp->second->mu, vecTauMu);
			getMatrixVecMulti(par->second, nod->second, vecZY);
			getVecVecAdd(vecTauMu, vecZY, vecSum);
			getMatrixVecMulti(inverseTauPost, vecSum, posterior->mu);
		}
		else
		{
			posterior->mu = clpp->second->mu;
		};

		//posterior rho
		if(!par->second.empty())  //set posterior to prior if no info
		{
			posterior->rho = clpp->second->rho + tot->second;
			//posterior phi
			getMatrixTransVecMulti(par->second, posterior->mu, vecZMu);
			getVecVecSub(nod->second, vecZMu, vecDiff1);
			getVecVecMulti(vecDiff1, nod->second, phiBit1);
			getVecVecSub(clpp->second->mu, posterior->mu, vecDiffMu);
			getVecVecMulti(vecDiffMu, vecTauMu, phiBit2);

			posterior->phi = clpp->second->phi + phiBit1 + phiBit2;
		}
		else
		{
			posterior->rho = clpp->second->rho;
			posterior->phi = clpp->second->phi;
		};

	};//end of prior (and posterior) loop

};

//! Calculates the local parameter prior.
void DiscreteDealNode::calcDiscreteLocalParameterPrior()
{
	unsigned int parentsGroup = 1;
	//unsigned int nodeParentsGroup = 1;
	map<unsigned int, double> localParas;

	//set initial levels
	setInitialLevels(); 
	
	//loop thro' node possibilities for parents
	do{
		parentsGroup = getParentsGroupNo();
	
		//get local parameters, loop thro' levels
		for(unsigned int lv = 1; lv <= getNoLevels(); ++lv)
		{
			level = lv;
			//nodeParentsGroup = getNodeAndParentsGroupNo();
			localParas[lv] = getMasterPriPara(); //set para value for this level given the current parents config
		};

		discreteLocalParaPris[parentsGroup] = localParas;
		
	}while(advanceParentLevels()); //advance to next parents config

};

//! Sets the level using the level name.
void DiscreteNode::setLevel(const string & lv)
{	
	//find the level	
	map<string, unsigned int>::const_iterator lev = data->levels.find(lv);

	if(lev == data->levels.end())
	{
		outErr("Level "); outErr(lv); outErr(" not found for node "); outErr(getName()); outErr(" when parsing parameter file.");
		exitErr("");
	};

	level = lev->second;
};

//! Returns level number of discrete node.
unsigned int DiscreteNode::getLevelNo(const string & strLevel)
{
	map<string, unsigned int>::const_iterator dl = data->levels.find(strLevel);

	if(dl == data->levels.end())
	{
		return data->addLevel(strLevel); //adds level if not found		
	};
		
	return dl->second;
};

//! Sets the discrete data to the current value.
void DiscreteNode::setDataValue()
{
	level = *dataIter;
};

//! Sets the discrete data to the initial value.
void DiscreteNode::setInitialDataValue()
{
	dataIter = data->values.begin(); 
	missIter = networkMissingData->missing.begin();

	level = *dataIter;
};

//! Advances data iterators.
void DiscreteNode::nextDataValueIter()
{
	dataIter++;
	missIter++;
};

//! Returns true if node data or other node data is missing.
bool DiscreteNode::nodeDataIsMissing()
{
	if(*missIter) return true;

	return false;
};

//! Sets the discrete data to the initial value.
void DiscreteNode::setInitialDataValue2()
{
	dataIter2 = data->values.begin(); 
	missIter2 = data->missingValues.begin();
	wasImpIter2 = data->wasImputed.begin();
};

//! Advances data iterators.
void DiscreteNode::nextDataValueIter2()
{
	dataIter2++;
	missIter2++;
	wasImpIter2++;
};

//! Returns true if node data is missing.
bool DiscreteNode::nodeDataIsMissing2()
{
	if(*missIter2) return true;

	return false;
};

//! Sets the discrete data to the initial value.
void DiscreteNode::setInitialDataValue3()
{
	dataIter3 = data->values.begin(); 
	missIter3 = data->missingValues.begin();
};

//! Advances data iterators.
void DiscreteNode::nextDataValueIter3()
{
	dataIter3++;
	missIter3++;
};

//! Returns true if node data is missing.
bool DiscreteNode::nodeDataIsMissing3()
{
	if(*missIter3) return true;

	return false;
};

//! Sets imputed value for discrete data.
void DiscreteNode::setImputedDataDis(const unsigned int & val, const bool & updateImpImmed)
{
	*dataIter2 = val;
	if(updateImpImmed) *missIter2 = false; //update later poss
	*wasImpIter2 = true;
};

//! Replaces missing values with randomly sampled values that are not missing.
void DiscreteNode::randomlyFillMissingValues()
{
	list<bool>::iterator m = data->missingValues.begin();
	for(list<unsigned int>::iterator v = data->values.begin(); v != data->values.end(); ++v, ++m)
	{
		if(*m)
		{
			*v = data->getRandomLevel();
			*m = false;
		};
	};
};

//! Calculates the local parameter posterior.
void DiscreteDealNode::calcDiscreteLocalParameterPost()
{
	unsigned int parentsGrpNo;
	map<unsigned int, map<unsigned int, double> >::iterator dlpp;
	map<unsigned int, double>::iterator pp; 
	
	//set up initial values..
	setInitialDataValue();

	for(map<unsigned int, Node *>::iterator pai = parents.begin(); pai != parents.end(); ++pai)
	{
		pai->second->setInitialDataValue();
	};

	//firstly copy over the prior
	discreteLocalParaPosts = discreteLocalParaPris;

	//next values...
	do
	{
		//add to the parameter posterior 
		parentsGrpNo = getParentsGroupNo();

		if(!nodeDataIsMissing())
		{

			dlpp = discreteLocalParaPosts.find(parentsGrpNo);
			
			if(dlpp != discreteLocalParaPosts.end())
			{
				pp = dlpp->second.find(level);
				if(pp != dlpp->second.end())
				{
					pp->second++; //add one to this parameter
				}
				else
				{
					dlpp->second[level] = 1;
				};
			}
			else
			{
				out("Cannot find posterior values for node "); out(getDisplayName()); out("!\n");
				exitErr("");
			};
			
		};

		//advance data to next reading
		nextDataValueIter();
		
		if(dataIter != data->values.end())
		{
			setDataValue();

			//set parent data iterators and data
			for(map<unsigned int, Node *>::iterator pa = parents.begin(); pa != parents.end(); ++pa)
			{
				pa->second->nextDataValueIter();
				pa->second->setDataValue();
			};
		};

	}while(dataIter != data->values.end());

};

//! Calculates the score contribution given by this discrete node.
double DiscreteDealNode::calcScoreBit(Network * network)
{
	double scoreBit = 0;
	map<unsigned int, map<unsigned int, double> >::iterator pri = discreteLocalParaPris.begin(); 
	map<unsigned int, map<unsigned int, double> >::iterator post = discreteLocalParaPosts.begin();
	map<unsigned int, double>::iterator prilv;
	map<unsigned int, double>::iterator postlv;
	double priorParaSum;
	double postParaSum;

	do{

		priorParaSum = 0;
		postParaSum = 0;
		prilv = pri->second.begin();
		postlv = post->second.begin();
		do{
			priorParaSum += prilv->second;
			postParaSum += postlv->second;
			scoreBit += gamln(&postlv->second) - gamln(&prilv->second);
			++prilv;
			++postlv;
		}while(prilv != pri->second.end());

		scoreBit += gamln(&priorParaSum) - gamln(&postParaSum);

		++pri;
		++post;
	}while(pri != discreteLocalParaPris.end());

	return scoreBit;
};

//! Calculates node contribution to score when a nodes parents are updated.
double DiscreteDealNode::calcUpdatedScoreBit(Network * network)
{
	NetworkDeal* networkDeal = dynamic_cast<NetworkDeal*>(network);
	if(networkDeal == 0) exitErr("Wrong kind of network detected, should be \"deal\"!");

	clearPriorsAndPosteriors();

	//update prior local distribution
	setDefaultInitPrior();

	//update the master prior
	calcDiscreteMasterPrior(networkDeal); //assumes prior dist for each node is the same for each marginal dist, ie each different set of parent levels (and different parents even)

	//update discreteLocalParaPris, the local parameter prior
	calcDiscreteLocalParameterPrior();

	//update discreteLocalParaPosts, the local parameter posterior
	calcDiscreteLocalParameterPost();

	return calcScoreBit(network);
};

//! Calculates node contribution to score when a nodes parents are updated.
double CtsDealNode::calcUpdatedScoreBit(Network * network)
{
	NetworkDeal* networkDeal = dynamic_cast<NetworkDeal*>(network);
	if(networkDeal == 0) exitErr("Wrong kind of network detected, should be \"deal\"!");

	clearPriorsAndPosteriors();

	//update prior local distribution
	setDefaultInitPrior();

	//update the master prior
	calcCtsMasterPrior(networkDeal); //assumes prior dist for each node is the same for each marginal dist, ie each different set of parent levels (and different parents even)

	//update ctsLocalParaPris, the local parameter prior
	calcCtsLocalParameterPrior(networkDeal);

	//update ctsLocalParaPosts, the local parameter posterior - no need to do this for cts nodes
	//calcCtsLocalParameterPost(network);

	return calcScoreBit(network);
};

//! Calculates the score contribution given by this cts node.
double CtsDealNode::calcScoreBit(Network * network)
{
	if(totalInData.empty()) setupNodeAndParentData();

	double scoreBit = 0;
	map<unsigned int, CtsLocalPara * >::const_iterator pri = ctsLocalParaPris.begin();
	//map<unsigned int, CtsLocalPara * >::const_iterator post = ctsLocalParaPosts.begin();
	double priRhoHalf, priRhoPlus1Half;//, postRhoHalf;

	map<unsigned int, unsigned int>::iterator tot = totalInData.begin();
	map<unsigned int, list< list<double> > >::iterator par = parentDataMat.begin();
	map<unsigned int, list<double> >::iterator nod = nodeData.begin();

	list<double>::iterator yit;
	list< list<double> >::iterator zit;
	list<double> vecZTau;

	double logpi = log(pi);
	double priorBit;
	double zTauz;
	double zMu, Q;

	do{
		priRhoHalf = 0.5*pri->second->rho;
		priRhoPlus1Half = 0.5*(pri->second->rho+1);
		priorBit = gamln(&priRhoPlus1Half) - gamln(&priRhoHalf);
	
		//loop thro data
		zit = par->second.begin();
		for(yit = nod->second.begin(); yit != nod->second.end(); ++yit, ++zit)
		{
			getVecMatrixMulti(*zit, pri->second->tauInverse, vecZTau);	
			getVecVecMulti(vecZTau, *zit, zTauz);

			scoreBit += priorBit - 0.5*(log(pri->second->phi) + logpi + log(1+zTauz));

			getVecVecMulti(*zit, pri->second->mu, zMu);

			Q = *yit - zMu;
			Q = (Q*Q)/(pri->second->phi*(1 + zTauz));

			scoreBit -= priRhoPlus1Half*log(1+Q);
		};

		//loop thro' data
		++par;
		++nod;
		++tot;
	
		++pri;		
	}while(pri != ctsLocalParaPris.end()); //end of loop thro different disrete parent groups

	return scoreBit;
};

//! Outputs the edges of a network.
void Node::outputEdges(ofstream & outfile, const unsigned int & nodeID)
{
	if(edgeSignifs.size() != parents.size())
	{
		exitErr("Problem calculating edge signifiances!");
	};

	map<unsigned int, pair<double, unsigned int> >::const_iterator es = edgeSignifs.begin();
	for(map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end(); ++pa, ++es)
	{
		if(!pa->second->getIsFactorChildNode())
		{
			outfile << pa->first << " " << nodeID << " " << es->second.first << "\n";
		};
	};

};

//! Returns the edge probability.
double Node::getEdgeProb(const unsigned int & parentNodeID)
{
	map<unsigned int, pair<double, unsigned int> >::const_iterator es = edgeSignifs.find(parentNodeID);

	if(es == edgeSignifs.end())
	{
		exitErr("Problem calculating edge probability!");
	};

	return getQvalueChiSq(es->second.first, (double)es->second.second);
};

//! Returns the edge significance.
pair<double, unsigned int> Node::getEdgeSign(const unsigned int & parentNodeID)
{
	map<unsigned int, pair<double, unsigned int> >::const_iterator es = edgeSignifs.find(parentNodeID);

	if(es == edgeSignifs.end())
	{
		exitErr("Problem calculating edge probability!");
	};

	return es->second;
};

//! Returns the number of parameters for a cts node.
unsigned int CtsNode::getNumberOfParameters()
{
	unsigned int noParameters = 2; //constant and variance
	unsigned int noLevelCombos = 1;

	for(map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end(); ++pa)
	{
		if(!pa->second->getIsDiscreteNode()) noParameters++;
		noLevelCombos *= pa->second->getNoLevels();
	};

	return noLevelCombos*noParameters;
};

//! Returns the number of parameters for a discrete node.
unsigned int DiscreteNode::getNumberOfParameters()
{
	unsigned int noLevelCombos = getNoLevels() - 1;

	for(map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end(); ++pa)
	{
		noLevelCombos *= pa->second->getNoLevels();
	};

	return noLevelCombos;
};

//! Outputs the edge names.
void Node::outputEdgeNames(ofstream & outfile)
{
	for(map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end(); ++pa)
	{
		if(!pa->second->getIsFactorChildNode())
		{
			outfile << pa->second->getDisplayName() << " " << getDisplayName() << "\n";
		};
	};
};

//! Outputs the prior local probability distributions for a discrete node.
void DiscreteNode::outputPriorLocalProbDists(ofstream & priorsFile)
{
	unsigned int noLevels = getNoLevels(); 

	if(getIsSNPNode()) priorsFile << "DISCRETE SNP NODE: ";
	else priorsFile << "DISCRETE NODE: ";
	priorsFile << getDisplayName() << "\n";
	
	bool hasDiscreteParents = localParaNames.size() > 1;

	if(hasDiscreteParents) outDiscreteParents(priorsFile, false);

	for(map<unsigned int, string>::const_iterator ln = localParaNames.begin(); ln != localParaNames.end(); ++ln)
	{
		if(hasDiscreteParents) priorsFile << " " << ln->second << ":\n";

		for(unsigned int lv = 1; lv <= noLevels; ++lv) priorsFile << "  " << getLevelName(lv) <<": " << uniformProb << "\n";
	};
	priorsFile << "\n";
};

//void CtsNode::outputPriorLocalProbDistsTEST()
//{
//	if(getIsSNPNode()) cout << "CONTINUOUS SNP NODE: ";
//	else  cout << "CONTINUOUS NODE: ";
//	 cout << getDisplayName();
//	if(getIsFactorNode() || getIsFactorChildNode())  cout << "(" << getFactorName() << ")";
//	 cout << "\n";
//
//	bool discreteGroups = !nodeParentNames.empty() && nodeParentNames.begin()->second != "";
//
//	if(nodeParentNames.size() != ctsPriLocalDists.size()) exitErr("Problem with discrete node and parents name (for a continuous node)!");
//
//	//if(discreteGroups) outDiscreteParents(priorsFile, true);
//
//	map<unsigned int, string>::const_iterator npn = nodeParentNames.begin();
//
//	map<unsigned int, Node *>::const_iterator pa;
//
//	
//	for(map<unsigned int, CtsPriLocalDist *>::const_iterator cp = ctsPriLocalDists.begin(); cp != ctsPriLocalDists.end(); ++cp, ++npn)
//	{
//		if(discreteGroups)  cout << " " << npn->second << ":\n";
//		if(getIsFactorNode() || getIsFactorChildNode())  cout << "  FACTOR: " << getFactorName() << "\n";
//		 cout << "  Intercept: " << cp->second->intercept << "\n";
//		 cout << "  Coefficients: ";
//		pa = parents.begin();
//		for(map<unsigned int, double>::const_iterator cf = cp->second->coeffs.begin(); cf != cp->second->coeffs.end(); ++cf)
//		{
//			while(pa->first != cf->first) {++pa; if(pa == parents.end()) exitErr("Problem with coefficient names of continuous node for prior local probability function!");};
//			 cout << pa->second->getDisplayName();
//			if(pa->second->getIsFactorNode() || pa->second->getIsFactorChildNode()) cout << "(" << pa->second->getFactorName() << ")";
//			cout << ": " << cf->second << " ";
//			
//		};
//		 cout << "\n";
//		 cout << "  Mean: " <<  cp->second->mean << "\n";
//		 cout << "  Variance: " << cp->second->variance << "\n";
//	};
//
//	 cout << "\n";
//
//};

//! Outputs the prior local probability distributions for a cts node.
void CtsNode::outputPriorLocalProbDists(ofstream & priorsFile)
{
	
	if(getIsSNPNode()) priorsFile << "CONTINUOUS SNP NODE: ";
	else priorsFile << "CONTINUOUS NODE: ";
	priorsFile << getDisplayName();
	if(getIsFactorNode() || getIsFactorChildNode()) priorsFile << "(" << getFactorName() << ")";
	priorsFile << "\n";

	bool discreteGroups = !nodeParentNames.empty() && nodeParentNames.begin()->second != "";

	if(nodeParentNames.size() != ctsPriLocalDists.size()) exitErr("Problem with discrete node and parents name (for a continuous node)!");

	if(discreteGroups) outDiscreteParents(priorsFile, true);

	map<unsigned int, string>::const_iterator npn = nodeParentNames.begin();

	map<unsigned int, Node *>::const_iterator pa;

	
	for(map<unsigned int, CtsPriLocalDist *>::const_iterator cp = ctsPriLocalDists.begin(); cp != ctsPriLocalDists.end(); ++cp, ++npn)
	{
		if(discreteGroups) priorsFile << " " << npn->second << ":\n";
		if(getIsFactorNode() || getIsFactorChildNode()) priorsFile << "  FACTOR: " << getFactorName() << "\n";
		priorsFile << "  Intercept: " << cp->second->intercept << "\n";
		priorsFile << "  Coefficients: ";
		pa = parents.begin();
		for(map<unsigned int, double>::const_iterator cf = cp->second->coeffs.begin(); cf != cp->second->coeffs.end(); ++cf)
		{
			while(pa->first != cf->first) {++pa; if(pa == parents.end()) exitErr("Problem with coefficient names of continuous node for prior local probability function!");};
			priorsFile << pa->second->getDisplayName();
			if(pa->second->getIsFactorNode() || pa->second->getIsFactorChildNode()) priorsFile << "(" << pa->second->getFactorName() << ")";
			priorsFile << ": " << cf->second << " ";
			
		};
		priorsFile << "\n";
		priorsFile << "  Mean: " <<  cp->second->mean << "\n";
		priorsFile << "  Variance: " << cp->second->variance << "\n";
	};

	priorsFile << "\n";
};

//! Outputs the master priors for a discrete node from a deal network.
void DiscreteDealNode::outputMasterPriors(ofstream & priorsFile)
{
	if(getIsSNPNode()) priorsFile << "DISCRETE SNP NODE: ";
	else priorsFile << "DISCRETE NODE: ";
	priorsFile << getDisplayName() << "\n";
	
	bool hasDiscreteParents = nodeParentNames.size() > 1;

	if(hasDiscreteParents) outDiscreteParents(priorsFile, true);

	if(discreteMasterParaPri.size() != nodeParentNames.size()) exitErr("Problem with discrete master priors!");
	map<unsigned int, string>::const_iterator npn = nodeParentNames.begin(); 
	
	for(map<unsigned int, double>::const_iterator dm = discreteMasterParaPri.begin(); dm != discreteMasterParaPri.end(); ++dm, ++npn)
	{
		priorsFile << " " << npn->second << ": " << dm->second << "\n";
	};
	priorsFile << "\n";
};

//! Outputs the local parameter priors for a discrete node from a deal network.
void DiscreteDealNode::outputLocalParameterPriors(ofstream & priorsFile)
{
	bool hasDisParents = false;
	if(!discreteLocalParaPris.empty() && discreteLocalParaPris.begin()->first != 0) hasDisParents = true;

	if(localParaNames.size() != discreteLocalParaPris.size()) exitErr("Problem with discrete local parameter priors!");

	map<unsigned int, string>::const_iterator dlpn = localParaNames.begin();

	if(getIsSNPNode()) priorsFile << "DISCRETE SNP NODE: ";
	else priorsFile << "DISCRETE NODE: ";
	priorsFile << getDisplayName() << "\n";

	if(hasDisParents) outDiscreteParents(priorsFile);
		
	for(map<unsigned int, map<unsigned int, double> >::const_iterator pri = discreteLocalParaPris.begin(); pri != discreteLocalParaPris.end(); ++pri, ++dlpn)
	{
		if(hasDisParents) priorsFile << " " << dlpn->second << ":\n";
		for(map<unsigned int, double>::const_iterator para = pri->second.begin(); para != pri->second.end(); ++para)
		{
			priorsFile << "  " << getLevelName(para->first) << ": "<< para->second <<"\n";
		};
		
	};
	priorsFile << "\n";
};

//! Outputs the local parameter posteriors for a discrete node from a deal network.
void DiscreteDealNode::outputLocalParameterPosteriors(ofstream & posteriorsFile)
{
	bool hasDisParents = false;
	if(!discreteLocalParaPris.empty() && discreteLocalParaPris.begin()->first != 0) hasDisParents = true;

	if(localParaNames.size() != discreteLocalParaPris.size()) exitErr("Problem with discrete local parameter posteriors!");

	map<unsigned int, string>::const_iterator dlpn = localParaNames.begin();

	if(getIsSNPNode()) posteriorsFile << "DISCRETE SNP NODE: ";
	else posteriorsFile << "DISCRETE NODE: ";
	posteriorsFile << getDisplayName() << "\n";

	if(hasDisParents) outDiscreteParents(posteriorsFile);

	for(map<unsigned int, map<unsigned int, double> >::const_iterator po = discreteLocalParaPosts.begin(); po != discreteLocalParaPosts.end(); ++po, ++dlpn)
	{
		if(hasDisParents) posteriorsFile << " " << dlpn->second << ":\n";
		for(map<unsigned int, double>::const_iterator para = po->second.begin(); para != po->second.end(); ++para)
		{
			posteriorsFile << "  " << getLevelName(para->first) << ": "<< para->second <<"\n";
		};
		
	};
	posteriorsFile << "\n";
};

//! Outputs the master priors for a cts node from a deal network.
void CtsDealNode::outputMasterPriors(ofstream & priorsFile)
{
	if(getIsSNPNode()) priorsFile << "CONTINUOUS SNP NODE: ";
	else priorsFile << "CONTINUOUS NODE: ";
	priorsFile << getDisplayName() << "\n";
	
	bool hasDiscreteParents = !nodeParentNames.empty() && nodeParentNames.begin()->second != "";

	if(hasDiscreteParents) outDiscreteParents(priorsFile, true);

	if(ctsMasterPriDists.size() != nodeParentNames.size()) exitErr("Problem with continous master priors!");
	map<unsigned int, string>::const_iterator npn = nodeParentNames.begin(); 
	
	map<unsigned int, Node *>::const_iterator pa;		
	list<double>::const_iterator m;

	for(map<unsigned int, CtsMultiDistMaster *>::const_iterator cmd = ctsMasterPriDists.begin(); cmd != ctsMasterPriDists.end(); ++cmd, ++npn)
	{
		if(hasDiscreteParents) priorsFile << npn->second << ":\n";

		priorsFile << "  Means:\n   ";
		pa = parents.begin();
		m = cmd->second->means.begin();
		priorsFile << getDisplayName() << ": " << *m << " ";
		++m;
		for(; m != cmd->second->means.end(); ++m)
		{
			while(pa->second->getNodeType() != "c") {++pa; if(pa == parents.end()) exitErr("Problem with coefficient names of continuous node for master prior!");};
			priorsFile << pa->second->getDisplayName() << ": " << *m << " ";
			++pa;
		};

		priorsFile << "\n";
		priorsFile << "  Covariance matrix:\n";
		for(list< list<double> >::const_iterator co1 = cmd->second->covariances.begin(); co1 != cmd->second->covariances.end(); ++co1)
		{
			priorsFile << "   ";
			for(list<double>::const_iterator co2 = co1->begin(); co2 != co1->end(); ++co2)
			{
				priorsFile << *co2 << " ";
			};
			priorsFile << "\n";
		};

	};
	priorsFile << "\n";
};

//! Write to priors file prior paramters of a cts node. 
void CtsDealNode::outputLocalParameterPriors(ofstream & priorsFile)
{
	
	bool hasDisParents = false;
	if(!ctsLocalParaPris.empty() && localParaNames.begin()->second != "") hasDisParents = true;

	if(localParaNames.size() != ctsLocalParaPris.size()) exitErr("Problem with continuous local parameter priors!");

	map<unsigned int, string>::const_iterator dlpn = localParaNames.begin();

	if(getIsSNPNode()) priorsFile << "CONTINUOUS SNP NODE: ";
	else priorsFile << "CONTINUOUS NODE: ";
	priorsFile << getDisplayName() << "\n";

	if(hasDisParents) outDiscreteParents(priorsFile);

	for(map<unsigned int, CtsLocalPara * >::const_iterator pri = ctsLocalParaPris.begin(); pri != ctsLocalParaPris.end(); ++pri, ++dlpn)
	{
		if(hasDisParents) priorsFile << " " << dlpn->second << ":\n";		
		priorsFile << "  mu: ";
		for(list<double>::const_iterator m = pri->second->mu.begin(); m != pri->second->mu.end(); ++m)
		{
			priorsFile << *m << " ";
		};
		priorsFile << "\n";
		priorsFile << "  phi: " << pri->second->phi << "\n";
		priorsFile << "  rho: " << pri->second->rho << "\n";
		priorsFile << "  tau:\n";
		for(list<list<double> >::const_iterator t = pri->second->tau.begin(); t != pri->second->tau.end(); ++t)
		{
			priorsFile << "   ";
			for(list<double>::const_iterator tt = t->begin(); tt != t->end(); ++tt)
			{
				priorsFile << *tt << " ";
			};		
			priorsFile << "\n";
		};

	};
	priorsFile << "\n";
};

//! Outputs the local parameter posteriors for a cts node from a deal network.
void CtsDealNode::outputLocalParameterPosteriors(ofstream & posteriorsFile)
{
	bool hasDisParents = false;
	if(!ctsLocalParaPris.empty() && localParaNames.begin()->second != "") hasDisParents = true;

	if(localParaNames.size() != ctsLocalParaPris.size()) exitErr("Problem with continuous local parameter priors!");

	map<unsigned int, string>::const_iterator dlpn = localParaNames.begin();

	if(getIsSNPNode()) posteriorsFile << "CONTINUOUS SNP NODE: ";
	else posteriorsFile << "CONTINUOUS NODE: ";
	posteriorsFile << getDisplayName() << "\n";

	if(hasDisParents) outDiscreteParents(posteriorsFile);

	for(map<unsigned int, CtsLocalPara * >::const_iterator post = ctsLocalParaPosts.begin(); post != ctsLocalParaPosts.end(); ++post, ++dlpn)
	{
		if(hasDisParents) posteriorsFile << " " << dlpn->second << ":\n";		
		posteriorsFile << "  mu: ";
		for(list<double>::const_iterator m = post->second->mu.begin(); m != post->second->mu.end(); ++m)
		{
			posteriorsFile << *m << " ";
		};
		posteriorsFile << "\n";
		posteriorsFile << "  phi: " << post->second->phi << "\n";
		posteriorsFile << "  rho: " << post->second->rho << "\n";
		posteriorsFile << "  tau:\n";
		for(list<list<double> >::const_iterator t = post->second->tau.begin(); t != post->second->tau.end(); ++t)
		{
			posteriorsFile << "   ";
			for(list<double>::const_iterator tt = t->begin(); tt != t->end(); ++tt)
			{
				posteriorsFile << *tt << " ";
			};		
			posteriorsFile << "\n";
		};

	};

	posteriorsFile << "\n";
};

//! Clears priors and posteriors for BNlearn discrete node.
void DiscreteBNLearnNode::clearPriorsAndPosteriors()
{
	if(freeClearMemory)
	{
		map<unsigned int, unsigned int>().swap(parentDataCounts);
		map<unsigned int, unsigned int>().swap(nodeParentDataCounts); 	
	}
	else
	{
		parentDataCounts.clear();
		nodeParentDataCounts.clear();
	};
};

//! Updates calculation of score bit for BNlearn discrete node.
double DiscreteBNLearnNode::calcUpdatedScoreBit(Network * network)
{
	NetworkBNLearn* networkDeal = dynamic_cast<NetworkBNLearn*>(network);
	if (networkDeal == 0) exitErr("Wrong kind of network detected, should be \"bnlearn\"!");

	clearPriorsAndPosteriors();

	//if(network->getUseMissingDataAtNodeLevel()) updateMissingDataAtNodeLevel(); //obselete

	//update prior local distribution
	setDefaultInitPrior();

	//update data counts
	updateDataCounts();
	
	return calcScoreBit(network);
};

//! Updates data count for node+parent config for discete bnlearn node.
void DiscreteBNLearnNode::updateDataCounts()
{
	unsigned int parentsGrpNo;
	unsigned int nodeParentsGrpNo;
	map<unsigned int, unsigned int>::iterator dc;
	totalData = 0;

	//set up initial values..
	setInitialDataValue();

	for(map<unsigned int, Node *>::iterator pai = parents.begin(); pai != parents.end(); ++pai)
	{
		pai->second->setInitialDataValue();
	};

	//next values...
	while(dataIter != data->values.end())
	{	
		//add to the parameter posterior 
		parentsGrpNo = getParentsGroupNo();
		nodeParentsGrpNo = getNodeAndParentsGroupNo();
	
		//add to data counts
		if(!nodeDataIsMissing()) //decide if missing data
		{
			dc = parentDataCounts.find(parentsGrpNo);

			if (dc != parentDataCounts.end())
			{
				parentDataCounts[parentsGrpNo]++; 			
			}
			else 
			{
				parentDataCounts[parentsGrpNo] = 1;			
			};

			dc = nodeParentDataCounts.find(nodeParentsGrpNo);

			if(dc != nodeParentDataCounts.end())
			{
				nodeParentDataCounts[nodeParentsGrpNo]++; 	
			}
			else
			{
				nodeParentDataCounts[nodeParentsGrpNo] = 1;
			};
			
			totalData++;
		};

		//advance data to next reading
		nextDataValueIter();

		if(dataIter != data->values.end())
		{
			setDataValue();

			//set parent data iterators and data
			for(map<unsigned int, Node *>::iterator pa = parents.begin(); pa != parents.end(); ++pa)
			{
				pa->second->nextDataValueIter();
				pa->second->setDataValue();
			};
		};

	};

};

//! Calculates score bit for BNlearn cts node.
double DiscreteBNLearnNode::calcScoreBit(Network * network)
{
	
	double ans = 0;
	double n = (double)totalData;
	double nodePaISS = 0;// These set firstly for AIC 
	double paISS = 0;	
	double iss = 0;       
	double k = 1;    //scoreType == 2 //AIC     

	unsigned int scoreType = network->getScoreType();

	if(scoreType == 0) //loglik
	{
		k = 0;
	}
	else if(scoreType == 1) //BIC
	{
		k = 0.5*log(n);
	}
	//else if(scoreType == 2) //AIC
	//{
	//	k = 1;
	//}
	else if(scoreType == 5) //BIC prob
	{
		k = 0.5*log(n);
	}
	else if(scoreType == 3) //bayes
	{
		iss = network->getImaginarySampleSize();
		paISS = 1; 
		for(map<unsigned int, Node *>::iterator pa = parents.begin(); pa != parents.end(); ++pa) paISS *= (double)pa->second->getNoLevels();
		paISS = iss/paISS;
		nodePaISS = paISS/(double)getNoLevels();	
		k = 0;
	};

	unsigned int nodeParentsGroup = 1;
	unsigned int parentsGroup = 1;
	unsigned int noParas = getNoLevels() - 1;

	setInitialLevels(); 
	
	for(map<unsigned int, Node *>::iterator pa = parents.begin(); pa != parents.end(); ++pa)
	{		
		noParas *= pa->second->getNoLevels();
	};

	unsigned int nodeParentCount;
	unsigned int parentCount;
	
	//loop thro' node possibilities for node and parent
	do{
		nodeParentsGroup = getNodeAndParentsGroupNo();
		parentsGroup = getParentsGroupNo();

		
		nodeParentCount = getNodeParentDataCount(nodeParentsGroup);
		parentCount = getParentDataCount(parentsGroup);

		if(nodeParentCount > 0 && parentCount > 0) ans += (double)nodeParentCount * log((((double)nodeParentCount + nodePaISS) / ((double)parentCount + paISS)));

	}while(advanceAllLevels()); //advance to next node+parents config


	ans -= (noParas)*k; //penalty for parameters (1/2)*(no paras)*log(sample size)

	if(scoreType == 5) //BIC + log(prob of edges), BICprob
	{
		for(map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end(); ++pa)
		{
			ans += log(network->getEdgeCost(pa->first, nodeID)); //it is a prob here
		};
	};

	return ans;
};

//! Gets data count for a certain config of node and parents.
unsigned int DiscreteBNLearnNode::getNodeParentDataCount(unsigned int & nodeParentsGroup)
{
	map<unsigned int, unsigned int>::const_iterator dc = nodeParentDataCounts.find(nodeParentsGroup);

	if(dc != nodeParentDataCounts.end()) return dc->second;

	return 0;
};


//! Gets data count for a certain config of node and parents.
unsigned int DiscreteBNLearnNode::getParentDataCount(unsigned int & parentsGroup)
{
	map<unsigned int, unsigned int>::const_iterator dc = parentDataCounts.find(parentsGroup);

	if (dc != parentDataCounts.end()) return dc->second;

	return 0;
};

//! Clears priors and posteriors for BNlearn cts node.
void CtsBNLearnNode::clearPriorsAndPosteriors()
{
	clearLocalPriors();
	
	if(freeClearMemory) list<unsigned int>().swap(nodeParentIDsForData); else nodeParentIDsForData.clear();
};

//! Updates calculation of score bit for BNlearn cts node.
double CtsBNLearnNode::calcUpdatedScoreBit(Network * network)
{
	if(/*getIsSNPNode() ||*/ getIsFactorNode() || getIsFactorChildNode())
	{
		return 0; //may be factors with zero variance if taking subsamples etc. Factor nodes and SNPs do not have parents so no need to calc proper score
	};
		 
	clearPriorsAndPosteriors();

	//update prior local distribution
	setDefaultInitPrior(); //linear reg

	return calcScoreBit(network);
};

//! Calculates score bit for BNlearn cts node.
double CtsBNLearnNode::calcScoreBit(Network * network)
{
	double ans = 0;
	unsigned int noParas = 2; //constant and variance
		
	unsigned int scoreType = network->getScoreType();
	unsigned int scoreFix = network->getScoreFix();

	setInitialLevels(); 

	//parent iterators
	list< list<double>::const_iterator > parIts;
	
	for (map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end(); ++pa)
	{
		if (pa->second->getNodeType() == "c")
		{
			parIts.push_back(pa->second->getCtsData()->values.begin());
			noParas++;
		};
	};

	//loop thro' cts data, and list saying which node+parent group or missing
	list<unsigned int>::const_iterator npi = nodeParentIDsForData.begin();
	map<unsigned int, CtsPriLocalDist *>::const_iterator cpld;  
	map<unsigned int, double>::const_iterator co; 
	double mean, diff, scoreBit;
	unsigned int noBadPts = 0;
	map<double, unsigned int> logVars;

	for(list<double>::const_iterator cd = data->values.begin(); cd != data->values.end(); ++cd, ++npi)
	{
		if(*npi != 0)
		{
			//firstly find fitted normal dist for this discrete parent config
			cpld = ctsPriLocalDists.find(*npi);
			if (cpld == ctsPriLocalDists.end()) exitErr("Problem with calculating cts node, "+ getDisplayName()+ ", score for bnlearn network!");

			if(cpld->second->coeffs.size() != parIts.size()) exitErr("Problem with calculating cts node, " + getDisplayName() + ", score (coeffs) for bnlearn network!");
			
			mean = cpld->second->intercept;
			co = cpld->second->coeffs.begin();

			//parent data
			for (list< list<double>::const_iterator >::iterator pai = parIts.begin(); pai != parIts.end(); ++pai, ++co)
			{
				mean += (*(*pai))*(co->second);
			};
			
			diff = *cd - mean;
	
			scoreBit = (diff*diff*0.5) / cpld->second->variance;

			if(cpld->second->variance == 0)
			{
				if(!warningDone && !network->getIsBootstrapData())
				{
					out("*** WARNING: variable ");
					out(getDisplayName());
					out(" has 0 variance, therefore a network score cannot be calculated! ***\n");
					warningDone = true;
				};			
			};

			if(scoreFix != 0 && scoreBit * 0 != 0)
			{
				noBadPts++;
			}
			else
			{
				ans -= scoreBit; //cache log value?

				//count how many times for this variance
				logVars[cpld->second->variance]++;
			};

		};//end, if missing data
		
		//advance iterators for parent data
		for(list< list<double>::const_iterator >::iterator pai = parIts.begin(); pai != parIts.end(); ++pai)
		{
			++(*pai);
		};

	};//end of cts data loop

	//add bit to log like for variance, the number of times it occurs in data
	for(map<double, unsigned int>::const_iterator lv = logVars.begin(); lv != logVars.end(); ++lv)
	{
		ans -= 0.5*log(6.28318530718 * lv->first)*(double)(lv->second); 
	};

	//double ans = 0;
	double n = getNoMissingNotMissing().second; //need the same number of points, n, at each node to weight 
	double nodePaISS = 0;// These set firstly for AIC 
	double paISS = 0;	
	double iss = 0;        
	double k = 1;         

	//adjust if bad pts
	if(scoreFix == 1 && noBadPts > 0)
	{
		ans *= n/(n-(double)noBadPts); //scale so that missing bad points are given the average contribution of other points
	};

	if(scoreType == 4)
	{
		k = 0.5*log(n); // as BIC

		for(map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end(); ++pa)
		{
			ans -= network->getEdgeCost(pa->first, nodeID)*k;
		};
	}
	else
	{
		//get number of parameters
		for (map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end(); ++pa)
		{
			if (pa->second->getNodeType() == "d")
			{	
				noParas *= pa->second->getNoLevels();
			};
		};
	
		if(scoreType == 0) //loglik
		{
			k = 0;
		}
		else if(scoreType == 1) //BIC
		{
			k = 0.5*log(n);
		}
		//else if(scoreType == 2) //AIC
		//{
		//	k = 1;
		//}
		else if(scoreType == 5) //BIC + log(prob of edges)
		{
			k = 0.5*log(n); // as BIC
		}
		else if(scoreType == 3) //bayes
		{
			//not actually done as too complex - use BIC instead
			iss = network->getImaginarySampleSize();
			paISS = 1; 
			for(map<unsigned int, Node *>::iterator pa = parents.begin(); pa != parents.end(); ++pa) paISS *= (double)pa->second->getNoLevels();
			paISS = iss/paISS;
			nodePaISS = paISS/(double)getNoLevels();			
			k = 0;
		};

		ans -= (noParas)*k; //penalty for parameters (1/2)*(no paras)*log(sample size)
	
		if(scoreType == 5) //BIC + log(prob of edges) //BICprob
		{
			for(map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end(); ++pa)
			{
				ans += log(network->getEdgeCost(pa->first, nodeID)); //it is a prob here
				//double prob = network->getEdgeCost(pa->first, nodeID);
				//ans += 0.5*(log(prob) - log(1.0-prob));
			};
		};

	};//end of score penalty - not cost
	
	return ans;
};

//! Copies fitted parameters sim'd network
void CtsBNLearnNode::simDataCopyParas(const string & simTaskName, Node * copyingNode, Network * origNetwork, Network * copyingNetwork, const bool & createDifferentNodeData)
{
	copyingNode->setCtsParas(ctsPriLocalDists, simTaskName, origNetwork, copyingNetwork, createDifferentNodeData);
};

//! Simulates some data for one obs for a cts node.
void CtsBNLearnNode::simData(bool & snpDataIntegers)
{
	double value = 0;
	unsigned int parentsGroup = getParentsGroupNo();
	unsigned int nodeParentsGroup = getNodeAndParentsGroupNo();

	map<unsigned int, CtsPriLocalDist *>::const_iterator cpld;  
	double mean = 0;
	double sd;
	map<double, unsigned int> logVars;
	
	//firstly find fitted normal dist for this discrete parent config
	cpld = ctsPriLocalDists.find(nodeParentsGroup);
	if(cpld == ctsPriLocalDists.end()) exitErr("Problem with simulating data for cts node, "+ getDisplayName()+ "!\n If it has discrete parents then there may be a rare combination that provides no cts distribution to simulate from.\n If you are estimating recall and precision then try removing some nodes.");
		
	mean = cpld->second->intercept;
	sd = cpld->second->stDev;

	map<unsigned int, double>::const_iterator co = cpld->second->coeffs.begin();
	map<unsigned int, Node *>::const_iterator pa = parents.begin();

	//parent data
	while(pa != parents.end())
	{
		//advance to next cts node
		while(pa != parents.end() && pa->second->getIsDiscreteNode())
		{
			++pa;
		};

		if(pa == parents.end()) break;

		if(co == cpld->second->coeffs.end()) exitErr("Problem with simulating data for cts node, " + getDisplayName() + ", number of coefficients and parents mismatch!");

		if(co->first != pa->first)
		{
			exitErr("Problem with simulating data for cts node, " + getDisplayName() + ", coefficients do not match parents!");
		};

		mean += (pa->second->getValue())*(co->second);

		++pa;
		++co;
	};
			
	bool missing = false;

	if(mean*0 != 0 || sd*0 != 0)
	{
		missing = true;
	}
	else if(sd==0)
	{
		value = mean;
	}
	else
	{
		int which = 2; //set for drawing random number
		double p = (double)rand()/(double)RAND_MAX;
		double q = 1 - p;
		int status = 0;
		double bound = 0;
	
		cdfnor(&which, &p, &q, &value, &mean, &sd, &status, &bound);

		missing = (status != 0);
	};

	if(getIsSNPNode() && snpDataIntegers)
	{
		if(value < 0) value = 0;
		else if(value > 2) value = 2;
	
		value = double((int)(value + 0.5)); //round to nearest integer but restricted to 0, 1 or 2
	};

	//now set the data, not missing
	data->addData(value, missing);
	setValue(value);
};

//! Copies cts parameters from another node, used for setting up sim network
void CtsBNLearnNode::setCtsParas(map<unsigned int, CtsPriLocalDist *> & ctsParas, const string & simTaskName, Network * origNetwork, Network * copyingNetwork, const bool & createDifferentNodeData)
{
	CtsPriLocalDist * copiedCtsParas;// = new CtsPriLocalDist();
	map<unsigned int, CtsPriLocalDist *>::iterator ccp;

	AllNodeData * allNodeData = copyingNetwork->getAllNodeData();
	string parentNodeName;
	unsigned int parentCoeffID;
	unsigned int parentNodeID;

	for(map<unsigned int, CtsPriLocalDist *>::const_iterator cp = ctsParas.begin(); cp != ctsParas.end(); ++cp)
	{
		ccp = ctsPriLocalDists.find(cp->first);

		if(ccp != ctsPriLocalDists.end())
		{
			if(freeClearMemory) map<unsigned int, double>().swap(ccp->second->coeffs); else ccp->second->coeffs.clear();

			copiedCtsParas = ccp->second;
		}
		else
		{
			copiedCtsParas = new CtsPriLocalDist();
			ctsPriLocalDistsToDel.push_back(copiedCtsParas);
		};

		copiedCtsParas->intercept = cp->second->intercept;

		for(map<unsigned int, double>::const_iterator cf = cp->second->coeffs.begin(); cf != cp->second->coeffs.end(); ++cf)
		{		
			if(createDifferentNodeData) parentNodeName = simTaskName + "-" + origNetwork->getNetworkNode(cf->first)->getName(); //get coressponding node in copying network
			else parentNodeName = origNetwork->getNetworkNode(cf->first)->getName();

			parentNodeID = allNodeData->getNodeDataNumber(parentNodeName);
			parentCoeffID = copyingNetwork->getNetworkNode(parentNodeID)->getNodeID(); //node ids are different between networks, use name modified for sim network
			copiedCtsParas->coeffs[parentCoeffID] = cf->second;		
		};

		copiedCtsParas->variance = cp->second->variance;
		copiedCtsParas->stDev = cp->second->stDev;

		
		ctsPriLocalDists[cp->first] = copiedCtsParas;
	};

	//add to cache so that it gets deleted later
	string parentsName = "";
	for(map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end(); )
	{
		parentsName += toString(pa->first);
		++pa;
		if(pa != parents.end()) parentsName += ":";
	};

};

//! Adds some cts parameters as given by parameter file for sim'ing data.
void CtsBNLearnNode::setCtsParas(unsigned int & nodeParentGroupNo, CtsPriLocalDist * ctsParas)
{
	ctsPriLocalDists[nodeParentGroupNo] = ctsParas;
	ctsPriLocalDistsToDel.push_back(ctsParas);
};

//! Sets default cts parameters.
void CtsBNLearnNode::setDefaultParameters()
{
	CtsPriLocalDist * ctsParas;

	//default discrete nodes have 3 levels unless otherwise specified
	unsigned int noDisParentGroups = 1;
	set<unsigned int> ctsParentIDs;

	for(map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end(); ++pa)
	{
		if(pa->second->getIsDiscreteNode()) noDisParentGroups *= pa->second->getNoLevels();	
		else
		{
			ctsParentIDs.insert(pa->second->getNodeID());
		};
	};
	
	for(unsigned int ndParentGroupNo = 1; ndParentGroupNo <= noDisParentGroups; ++ndParentGroupNo)
	{
		ctsParas = new CtsPriLocalDist();
		ctsPriLocalDistsToDel.push_back(ctsParas);

		ctsParas->intercept = 10;
		ctsParas->variance = 1;
		ctsParas->stDev = 1;

		for(set<unsigned int>::const_iterator cp = ctsParentIDs.begin(); cp != ctsParentIDs.end(); ++cp)
		{
			ctsParas->coeffs[*cp] = 1;
		};

		ctsPriLocalDists[ndParentGroupNo] = ctsParas;
	};
};

//! Copies fitted parameters sim'd network
void DiscreteBNLearnNode::simDataCopyParas(const string & simTaskName, Node * copyingNode, Network * origNetwork, Network * copyingNetwork, const bool & createDifferentNodeData)
{
	copyingNode->setDiscreteParas(nodeParentDataCounts, parentDataCounts, totalData);
};

//! Copies the discrete data level names
void DiscreteBNLearnNode::simDataCopyDisLevels(Node * copyingNode)
{
	data->levels = copyingNode->getLevels();
};

//! Returns the name and id of the levels
map<string, unsigned int> DiscreteBNLearnNode::getLevels()
{
	return data->levels;
}

//! Copies fitted paras for discrete nodes
void DiscreteBNLearnNode::setDiscreteParas(map<unsigned int, unsigned int> & ndDataCnts, map<unsigned int, unsigned int> & parentDataCnts, unsigned int & totData)
{
	nodeParentDataCounts = ndDataCnts;
	parentDataCounts = parentDataCnts; 
	totalData = 0;// no actual data so don't copy this number over, not needed for sims
};

//! Sets default levels and probs for sim'ing data.
void DiscreteBNLearnNode::setDefaultParameters()
{

	for(unsigned int lv = 0; lv < simNoLevels; ++lv) data->levels[toString(lv)] = lv + 1;

	map<unsigned int, double> cumProbs;
	if(simNoLevels == 3 && getIsSNPNode())
	{
		cumProbs[1] = 0.25;
		cumProbs[2] = 0.75;
	}
	else
	{
		double aProb = 1.0/(double)(simNoLevels);
		for(unsigned int lv = 1; lv < simNoLevels; ++lv) cumProbs[lv] = aProb*lv;
	};

	unsigned int noDisParentGroups = 1;

	for(map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end(); ++pa)
	{
		if(pa->second->getIsDiscreteNode()) noDisParentGroups *= pa->second->getSimNoLevels();
	};

	if(noDisParentGroups == 1)
	{
		simDataCumLevelProbs[0] = cumProbs;
	}
	else
	{
		for(unsigned int ndParentGroupNo = 1; ndParentGroupNo <= noDisParentGroups; ++ndParentGroupNo)
		{
			simDataCumLevelProbs[ndParentGroupNo] = cumProbs; 
		};
	};
};

//! Simulates some data for one obs for a discrete node.
void DiscreteBNLearnNode::simData(bool & snpDataIntegers)
{
	map<unsigned int, double> cumLevelProbs; //level, prob -- depends on parents

	unsigned int parentsGroup = getParentsGroupNo();

	//try and find in the cache firstly
	map<unsigned int, map<unsigned int, double> >::const_iterator sdclp = simDataCumLevelProbs.find(parentsGroup);

	if(sdclp != simDataCumLevelProbs.end())
	{
		cumLevelProbs = sdclp->second;
	}
	else
	{
		unsigned int nodeParentsGroup;
		unsigned int parentCount;
		unsigned int nodeParentCount; 
		unsigned int noLevels = getNoLevels();
		double cumProb;
		double prevCumProb = 0;

		map<unsigned int, unsigned int>::const_iterator pdc = parentDataCounts.find(parentsGroup);
		if(pdc != parentDataCounts.end())
		{
			parentCount = pdc->second;
			
			for(unsigned int lv = 1; lv < noLevels; ++lv)
			{
				setLevel(lv);
				nodeParentsGroup = getNodeAndParentsGroupNo();
				nodeParentCount = getNodeParentDataCount(nodeParentsGroup);
				cumProb = (((double)nodeParentCount) / ((double)parentCount));	
				cumLevelProbs[lv] = prevCumProb + cumProb;
				prevCumProb = cumProb;
			};
		}
		else
		{			
			//exitErr("Could not find parent data count when simulating discrete data!");
			//use overall totals for any parent group as no info on this parent combo
			map<unsigned int, unsigned int> levelCounts;
			map<unsigned int, unsigned int>::iterator lc;
			unsigned int theLevel;
			unsigned int levelCount;
			unsigned int totalCounts = 0;

			for(map<unsigned int, unsigned int>::const_iterator pdc2 = nodeParentDataCounts.begin(); pdc2 != nodeParentDataCounts.end(); pdc2++)
			{
					theLevel = getLevelFromNodeAndParentsGroupNo(pdc2->first);
					lc = levelCounts.find(theLevel);
					if(lc == levelCounts.end()) levelCounts[theLevel] = pdc2->second;
					else lc->second += pdc2->second;

					totalCounts += pdc2->second;
			};
				
			for(unsigned int lv = 1; lv < noLevels; ++lv)
			{
				lc = levelCounts.find(lv);
				if(lc == levelCounts.end()) levelCount = 0; else levelCount = lc->second;

				cumProb = (((double)levelCount) / ((double)totalCounts));	
				cumLevelProbs[lv] = prevCumProb + cumProb;
				prevCumProb = cumProb;
			};
		};

		//add to cache
		simDataCumLevelProbs[parentsGroup] = cumLevelProbs;
	};

	//set the level randomly
	double ranNum = (double)rand()/(double)RAND_MAX;
	unsigned int randomLevel = 1;
	unsigned int prevLevel = 0;

	for(map<unsigned int, double>::const_iterator lvp = cumLevelProbs.begin(); lvp != cumLevelProbs.end(); )
	{
		if(ranNum <= lvp->second)
		{
			randomLevel = lvp->first;			
			break;
		};

		prevLevel = lvp->first;

		++lvp;

		if(lvp == cumLevelProbs.end()) {randomLevel = prevLevel + 1; break;};
	};

	//now set the data, not missing
	data->addSimData(randomLevel, false);
	setLevel(randomLevel);
};

//! Output posteriors for BNlearn discrete node.
void DiscreteBNLearnNode::outputPosteriors(ofstream & posteriorsFile)
{
	bool hasDisParents = localParaNames.size() > 1;
	
	if(getIsSNPNode()) posteriorsFile << "DISCRETE SNP NODE: ";
	else posteriorsFile << "DISCRETE NODE: ";
	posteriorsFile << getDisplayName() << "\n";

	if(hasDisParents) outDiscreteParents(posteriorsFile);

	unsigned int nodeParentsGroup = 1;
	unsigned int parentsGroup = 1;
	
	setInitialLevels(); 	

	unsigned int nodeParentCount;
	unsigned int parentCount;

	map<unsigned int, string>::const_iterator dlpn;
	string prevDisParentsName = "";

	//loop thro' node possibilities for node and parent
	do{
		nodeParentsGroup = getNodeAndParentsGroupNo();
		parentsGroup = getParentsGroupNo();

		//check if beginnning a new discrete parents group
		if(hasDisParents)
		{
			dlpn = localParaNames.find(parentsGroup);
			if(dlpn == localParaNames.end())
			{
				outErr("Problem outputting posterior for discrete node "); outErr(getName()); outErr("!");
				exitErr("");
			};

			if(dlpn->second != prevDisParentsName)
			{
				posteriorsFile << " " << dlpn->second << ":\n";
				prevDisParentsName = dlpn->second;
			};
		};

		nodeParentCount = getNodeParentDataCount(nodeParentsGroup);
		parentCount = getParentDataCount(parentsGroup);
		posteriorsFile << "  " << getLevelName(level) << ": ";
		if(nodeParentCount == 0) posteriorsFile << "0\n";
		else posteriorsFile << (((double)nodeParentCount) / ((double)parentCount)) <<"\n";

	}while(advanceAllLevels()); //advance to next node+parents config

	posteriorsFile << "\n";
};
