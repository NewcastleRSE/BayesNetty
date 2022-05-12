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


/*! \file Network.cpp
    \brief This file contains the methods the various analyse.   
*/

#include <iostream>
#include <sstream>

using namespace std; // initiates the "std" or "standard" namespace

#include "Network.h"
#include "main.h"
#include "Task.h"
#include "cdflib.h"
#include "Model.h"
#include "Utils.h"
#include <math.h>
#include <vector>

#ifdef USING_OPEN_MPI
#include "mpi.h"
#endif 

//! Returns (variance or) covariance for a two nodes.
double CtsMultiDistJoint::getCovariance(const unsigned int &nd1, const unsigned int &nd2)
{
	if(nd1==nd2)
	{
		map<unsigned int, double>::const_iterator v = variances.find(nd1);
		if(v != variances.end()) return v->second;
		else
		{
			exitErr("Problem finding variance in joint prior!");
		};
	};

	map<unsigned int, map<unsigned int, double> >::const_iterator cv1 = covariances.find(nd1);
	if(cv1 != covariances.end())
	{
		map<unsigned int, double>::const_iterator cv2 = cv1->second.find(nd2);

		if(cv2 != cv1->second.end()) return cv2->second;
	};	
	
	cv1 = covariances.find(nd2);
	if(cv1 != covariances.end())
	{
		map<unsigned int, double>::const_iterator cv22 = cv1->second.find(nd1);

		if(cv22 != cv1->second.end()) return cv22->second;
	};	

	return 0;
};

//! Delete and clear prior and posterior data.
void NetworkDeal::clearPriorsAndPosteriors()
{	
	if(freeClearMemory) map<unsigned int, double>().swap(discreteJointProbDist);  else discreteJointProbDist.clear();

	for(map<unsigned int, CtsMultiDistJoint *>::iterator cjp = ctsJointProbDist.begin(); cjp != ctsJointProbDist.end(); ++cjp) delete cjp->second;	

	if(freeClearMemory) map<unsigned int, CtsMultiDistJoint*>().swap(ctsJointProbDist); else ctsJointProbDist.clear();

	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		nd->second->clearPriorsAndPosteriors();
	};
};

//! Gets network node given the number.
Node * Network::getNetworkNode(unsigned int nodeNo) //const
{
	map<unsigned int, Node *>::const_iterator nd = allNetworkNodes.find(nodeNo);

	if(nd == allNetworkNodes.end())
	{
		//try to find in all nodes
		string nodeName = allNodeData->getNodeDataNameError(nodeNo);
		string mess;
		if(nodeName != "") mess = "Attempt to get network node, \"" + nodeName + "\", but no such node exists in the network!\nCheck that any blacklists (and such) only include nodes in the current network.\n";
		else mess = "Attempt to get network node but no such node exists!\n";
		
		//unsigned int arghh = nd->second->getFileNo() + 1; //deliberate crash for gdb debugging 

		exitErr(mess);
	};

	return nd->second;
};

//! Gets network node given the name.
Node * Network::getNoDataNetworkNode(const string & nodeName) const
{
	map<string, Node *>::const_iterator nd = noDataNetworkNodes.find(nodeName);

	if(nd == noDataNetworkNodes.end())
	{
		string mess = "Attempt to get network node, " + nodeName + ", but no such node exists!\n";
		exitErr(mess);
	};

	return nd->second;
};

//! Sets network structure given network string in format [A][B][C|A:B].
void Network::setNetwork(const string & networkStr, const string & preNodeName)
{
	
	bool foundNode;	
	string nodeParentsStr;	
	string nodeStr;
	string parentsStr, aParent;
	string networkString = networkStr;
	unsigned int pos1, pos2, posBar, posColon;
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
			if(posBar < nodeParentsStr.length()) nodeStr = preNodeName + nodeParentsStr.substr(0, posBar);
			else nodeStr = preNodeName + nodeParentsStr;

			if(preNodeName == "")
			{
				if(nodeExistsInit(nodeStr))
				{
					outErr("Node "); outErr(nodeStr); out(" is repeated in the network file!\n");
					exit(1);
				};
				
				addNodeInit(nodeStr);
			};

			//check if parents exist
			if(posBar < nodeParentsStr.length())
			{
				parentsStr = nodeParentsStr.substr(posBar+1);

				do{					
					posColon = (unsigned int)parentsStr.find_first_of(':');

					if(posColon < parentsStr.length()) aParent = preNodeName + parentsStr.substr(0, posColon);
					else aParent = preNodeName + parentsStr;

					//add edge
					if(!edgeExistsInit(aParent, nodeStr)) addEdgeInit(aParent, nodeStr);

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


	if(getNoNodes() == 0 && preNodeName == "") exitErr("Failed to load a network with any nodes!\n Check that you have set the correct input format for the network - should like this \"[A][B][C|A:B]\" for this network input option.");
};

//! Adds node to the (deal) network given the node name.
void NetworkDeal::addNode(const unsigned int & nodeNumber)
{
	Data * theData = allNodeData->getNodeData(nodeNumber);

	if(theData->dataType == 1) //discrete
	{
		DiscreteData * disData = static_cast<DiscreteData*>(theData);
		DiscreteDealNode * aNode = new DiscreteDealNode(theData->name, disData);
		allNetworkNodes[nodeNumber] = aNode;
		discreteNodes[nodeNumber] = aNode;
		aNode->setNodeID(nodeNumber);
		discreteDealNodes.insert(aNode);
	} 
	else if(theData->dataType == 2) //cts
	{
		CtsData * ctsData = static_cast<CtsData*>(theData);
		CtsDealNode * aNode = new CtsDealNode(theData->name, ctsData);
		allNetworkNodes[nodeNumber] = aNode;
		ctsNodes[nodeNumber] = aNode;
		aNode->setNodeID(nodeNumber);
		ctsDealNodes.insert(aNode);
	}
	else if(theData->dataType == 3)	addNoDataNode(theData->name, static_cast<NoData*>(theData));
	
};

//! Adds edges for a factor variable for all the other factor data.
void Network::addEdgesFactorChildren(const unsigned int & nodeNo1, Node * node2)
{
	Node * aFactorChildNode;
	
	list<unsigned int> childFactorNodeIDs = allNodeData->getFactorDataGroup(nodeNo1);

	for(list<unsigned int>::const_iterator cfni = childFactorNodeIDs.begin(); cfni != childFactorNodeIDs.end(); ++cfni)
	{
		aFactorChildNode = getNetworkNode(*cfni);
		node2->addParent(aFactorChildNode);
	};
};

//! Adds an edge to the network given the node names.
void Network::addEdgeInit(const string & nodeName1, const string & nodeName2)
{
	unsigned int nodeNumber1 = allNodeData->getNodeDataNumber(nodeName1);
	unsigned int nodeNumber2 = allNodeData->getNodeDataNumber(nodeName2);

	Node * aNode1 = getNetworkNode(nodeNumber1);
	Node * aNode2 = getNetworkNode(nodeNumber2);

	aNode2->addParent(aNode1);

	//if node is a factor node then also add the nodes for the other factors
	if(aNode1->getIsFactorNode()) addEdgesFactorChildren(nodeNumber1, aNode2);
};

//! Removes an edge to the network given the node names.
void Network::removeEdgeInit(const string & nodeName1, const string & nodeName2)
{
	unsigned int nodeNumber1 = allNodeData->getNodeDataNumber(nodeName1);
	unsigned int nodeNumber2 = allNodeData->getNodeDataNumber(nodeName2);

	removeEdge(nodeNumber1, nodeNumber2);
};

//! Adds an edge to the network given the node IDs.
void Network::addEdge(Node * node1, Node * node2)
{
	node2->addParent(node1);

	if(node1->getIsFactorNode()) addEdgesFactorChildren(node1->getNodeID(), node2);
};

//! Adds an edge to the network given the node IDs.
void Network::addEdge(const unsigned int & nodeNo1, const unsigned int & nodeNo2)
{
	Node * aNode1 = getNetworkNode(nodeNo1);
	Node * aNode2 = getNetworkNode(nodeNo2);
	aNode2->addParent(aNode1);

	if(aNode1->getIsFactorNode()) addEdgesFactorChildren(nodeNo1, aNode2);
};

//! Adds a random edge.
void Network::addRandomEdge()
{
	//do not add an existing edge, do not reverse an edge, or add a black edge
	unsigned int maxNodeID = allNetworkNodes.rbegin()->first;
	unsigned int minNodeID = allNetworkNodes.begin()->first;
	unsigned int range = maxNodeID - minNodeID + 1;

	unsigned int nodeID1, nodeID2;
	map<unsigned int, Node *>::const_iterator nd; //allNetworkNodes;
	bool edgeFound = true;
	unsigned int maxIts = 50; //make sure does not try too long to add an edge
	unsigned int iters = 0;

	do{
		
		do{
			
			do{
				nodeID1 = rand() % range + minNodeID;
				nd = allNetworkNodes.find(nodeID1);
			}while(nd == allNetworkNodes.end());

			do{
				nodeID2 = rand() % range + minNodeID;
				nd = allNetworkNodes.find(nodeID2);
			}while(nd == allNetworkNodes.end());

			//check edge is valid			
			if(edgeExists(nodeID1, nodeID2) || edgeExists(nodeID2, nodeID1)) edgeFound = false; //no existing edges	
			else if(!edgeAllowed(nodeID1, nodeID2)) edgeFound = false; //inc no black edges

			iters++;
		}while(!edgeFound && iters <= maxIts);

		if(edgeFound)
		{
			//add the edge 1 --> 2
			addEdge(nodeID1, nodeID2);

			if(hasLoop())
			{
				removeEdge(nodeID1, nodeID2);
				edgeFound = false;
			};
		};

	}while(!edgeFound && iters <= maxIts);
};

//! Changes a random edge.
void Network::changeRandomEdge()
{
	//Add an existing edge or reverse an edge, do not add a black edge
	unsigned int maxNodeID = allNetworkNodes.rbegin()->first;
	unsigned int nodeID1, nodeID2;
	map<unsigned int, Node *>::const_iterator nd; // allNetworkNodes;
	bool edgeFound = true;
	unsigned int maxIts = 10000; //50; //make sure does not try too long to change an edge
	unsigned int iters = 0;

	do{
		
		do{
			
			edgeFound = true;

			do{
				nodeID1 = rand() % maxNodeID + 1;
				nd = allNetworkNodes.find(nodeID1);
			}while(nd == allNetworkNodes.end());

			do{
				nodeID2 = rand() % maxNodeID + 1;
				nd = allNetworkNodes.find(nodeID2);
			}while(nd == allNetworkNodes.end());

			//check edge is valid			
			if(!edgeAllowed(nodeID1, nodeID2) || isWhiteEdge(nodeID1, nodeID2) || isWhiteEdge(nodeID2, nodeID1)) edgeFound = false; //cts --> discrete is not allowed etc			
			else if(edgeExists(nodeID1, nodeID2) && !edgeAllowed(nodeID2, nodeID1)) edgeFound = false; //no parents for this node 
			
			iters++;
		}while(!edgeFound && iters <= maxIts);

		if(edgeFound)
		{
			//change direction if exists
			if(edgeExists(nodeID1, nodeID2))
			{
				reverseEdge(nodeID1, nodeID2);

				if(hasLoop())
				{
					reverseEdge(nodeID2, nodeID1);
					edgeFound = false;
				};
			}
			else if(edgeExists(nodeID2, nodeID1))
			{
				reverseEdge(nodeID2, nodeID1);

				if(hasLoop())
				{
					reverseEdge(nodeID1, nodeID2);
					edgeFound = false;
				};
			}
			else
			{
				//add the edge 1 --> 2
				addEdge(nodeID1, nodeID2);

				if(hasLoop())
				{
					removeEdge(nodeID1, nodeID2);
					edgeFound = false;
				};
			};
		
		};

	}while(!edgeFound && iters <= maxIts);

};

//! Adds white node to network.
void Network::addWhiteNode(const string & nodeName)
{
	unsigned int nodeNumber = allNodeData->getNodeDataNumber(nodeName);
	whiteNodes.insert(nodeNumber);
};

//! Adds white edge node to network.
void Network::addWhiteEdge(const string & nodeName1, const string & nodeName2)
{
	unsigned int nodeNumber1 = allNodeData->getNodeDataNumber(nodeName1);
	unsigned int nodeNumber2 = allNodeData->getNodeDataNumber(nodeName2);
	whiteEdges.insert(make_pair(nodeNumber1, nodeNumber2));
};

//! Checks if an edge is white.
bool Network::isWhiteEdge(const unsigned int & nodeNo1, const unsigned int & nodeNo2)
{
	return whiteEdges.find(make_pair(nodeNo1, nodeNo2)) != whiteEdges.end();
};

//! Sets up which edges are allowed in the network.
void Network::setUpAllowedEdgeTypes(set<pair<unsigned int, unsigned int> > & allowedEdgeTypesFileNos, set<pair<unsigned int, unsigned int> > & notAllowedEdgeTypesFileNos)
{

	if(!allowedEdgeTypesFileNos.empty() && !notAllowedEdgeTypesFileNos.empty())
	{
		exitErr("You may not specify some network edge types to be OK and others not OK!"); 
	};

	set<unsigned int> discreteFileNos;
	set<unsigned int> ctsFileNos;	

	//set up sets of file nos, i.e. edge types
	for(map<unsigned int, DiscreteNode *>::const_iterator dnd = discreteNodes.begin(); dnd != discreteNodes.end(); ++dnd)
	{
		discreteFileNos.insert(dnd->second->getFileNo());		
	};

	for(map<unsigned int, CtsNode *>::const_iterator cnd = ctsNodes.begin(); cnd != ctsNodes.end(); ++cnd)
	{
		ctsFileNos.insert(cnd->second->getFileNo());		
	};

	//treat no data nodes as cts
	for(map<string, Node *>::const_iterator nnd = noDataNetworkNodes.begin(); nnd != noDataNetworkNodes.end(); ++nnd)
	{
		ctsFileNos.insert(nnd->second->getFileNo()); //will be 0	
	};

	if(!allowedEdgeTypesFileNos.empty())
	{
		//check for not allowed cts-->dis		
		for(set<unsigned int>::const_iterator cfc = ctsFileNos.begin(); cfc != ctsFileNos.end(); ++cfc)
		{		
			//cts-->dis
			for(set<unsigned int>::const_iterator dfc = discreteFileNos.begin(); dfc != discreteFileNos.end(); ++dfc)
			{
				if(allowedEdgeTypesFileNos.find(make_pair(*cfc, *dfc)) != allowedEdgeTypesFileNos.end())
				{				
					exitErr("Edges such that continuous --> discerte are invalid!");			
				};
			};
		};

		allowedEdgeTypes = allowedEdgeTypesFileNos;		
	}
	else
	{
		//add all by default except cts-->discrete
		for(set<unsigned int>::const_iterator df = discreteFileNos.begin(); df != discreteFileNos.end(); ++df)
		{
			//dis-->dis
			for(set<unsigned int>::const_iterator df2 = discreteFileNos.begin(); df2 != discreteFileNos.end(); ++df2)
			{
				allowedEdgeTypes.insert(make_pair(*df, *df2));
			};

			//dis-->cts
			for(set<unsigned int>::const_iterator cf2 = ctsFileNos.begin(); cf2 != ctsFileNos.end(); ++cf2)
			{
				allowedEdgeTypes.insert(make_pair(*df, *cf2));
			};
		};

		for(set<unsigned int>::const_iterator cf = ctsFileNos.begin(); cf != ctsFileNos.end(); ++cf)
		{
			//cts-->cts
			for(set<unsigned int>::const_iterator cf3 = ctsFileNos.begin(); cf3 != ctsFileNos.end(); ++cf3)
			{
				allowedEdgeTypes.insert(make_pair(*cf, *cf3));
			};

		};

		//remove not allowed edge types
		if(!notAllowedEdgeTypesFileNos.empty())
		{
			set<pair<unsigned int, unsigned int> >::const_iterator aet; 
			for(set<pair<unsigned int, unsigned int> >::const_iterator nae = notAllowedEdgeTypesFileNos.begin(); nae != notAllowedEdgeTypesFileNos.end(); ++nae)
			{
				aet = allowedEdgeTypes.find(*nae);
				if(aet != allowedEdgeTypes.end()) allowedEdgeTypes.erase(*aet);
			};
		};
	};

	//do not allow SNPs to have parents
	setSNPNodesAsNoParentsNodes();
};

//! Adds an edge type that is allowed in the network.
void Network::addNetworkEdgeType(const unsigned int & fileNo1, const unsigned int & fileNo2)
{
	allowedEdgeTypes.insert(make_pair(fileNo1, fileNo2));
};

//! Checks if an edge is allowed in the network or not, node1 --> node2.
bool Network::edgeAllowed(const unsigned int & nodeNo1, const unsigned int & nodeNo2)
{	
	Node * node1 = getNetworkNode(nodeNo1);
	Node * node2 = getNetworkNode(nodeNo2);

	//look up file nos for node types
	unsigned int fileNo1 = node1->getFileNo();
	unsigned int fileNo2 = node2->getFileNo();

	return nodeNo1 != nodeNo2 && allowedEdgeTypes.find(make_pair(fileNo1, fileNo2)) != allowedEdgeTypes.end()
		&& blackEdges.find(make_pair(nodeNo1, nodeNo2)) == blackEdges.end()
		&& (noChildrenNodes.find(nodeNo1) == noChildrenNodes.end() || isWhiteEdge(nodeNo1, nodeNo2))
		&& (noParentsNodes.find(nodeNo2) == noParentsNodes.end() || isWhiteEdge(nodeNo1, nodeNo2))
		&& !node1->getIsFactorChildNode() && !node2->getIsFactorChildNode()
		&& !node2->getIsFactorNode();
};

//! Adds a cost for an edge type.
void Network::setCostEdgeType(const unsigned int & fileNo1, const unsigned int & fileNo2, const double & cost)
{
	costEdgeTypes[make_pair(fileNo1, fileNo2)] = cost;
};

//! Adds a cost for a specific edge.
void Network::setCostEdge(const string & nodeName1, const string & nodeName2, const double & cost)
{
	unsigned nodeNumber1 = allNodeData->getNodeDataNumberInit(nodeName1);
	unsigned nodeNumber2 = allNodeData->getNodeDataNumberInit(nodeName2);

	costEdges[make_pair(nodeNumber1, nodeNumber2)] = cost;
};

//! Returns the cost of an edge
double Network::getEdgeCost(const unsigned int & parentNodeNo, const unsigned int & childNodeNo)
{
	//look in edges firstly
	map<pair<unsigned int, unsigned int>, double>::const_iterator ce = costEdges.find(make_pair(parentNodeNo, childNodeNo));
	if(ce != costEdges.end()) return ce->second;

	//now look in edge types
	Node * node1 = getNetworkNode(parentNodeNo);
	Node * node2 = getNetworkNode(childNodeNo);

	//look up file nos for node types
	unsigned int fileNo1 = node1->getFileNo();
	unsigned int fileNo2 = node2->getFileNo();

	map<pair<unsigned int, unsigned int>, double>::const_iterator cet = costEdgeTypes.find(make_pair(fileNo1, fileNo2));
	if(cet != costEdgeTypes.end()) return cet->second;
	
	//if(scoreType==5) //must be for prior prob of edges only way to weight edges now
	//{
		if(node1->getIsDiscreteNode() && !node2->getIsDiscreteNode()) return 1; //discrete to cts is prob 1

		return 0.5;
	//};

	//return 1; //default
};

//! Adds node with no data.
void Network::addNoDataNode(const string & nodeName, NoData * noDataData)
{
	string newNodeName = nodeName;
	list<string> nodeAndType = splitString(nodeName, "|");
	unsigned int simNodeType = 0; //0 = cts, 1 = discrete
	int noLevels = 3; //default
	bool isSNP = false;

	if(nodeAndType.size() >= 2)
	{
		list<string>::const_iterator nat = nodeAndType.begin();
		newNodeName = *nat;
		++nat;

		if(*nat == "cts" || *nat == "continuous" || *nat == "cts.snp" || *nat == "continuous.snp") simNodeType = 0;
		else if(*nat == "dis" || *nat == "discrete" || *nat == "dis.snp" || *nat == "discrete.snp") simNodeType = 1;

		if(*nat == "cts.snp" || *nat == "continuous.snp" || *nat == "dis.snp" || *nat == "discrete.snp") isSNP = true;

		if(nodeAndType.size() >= 3)
		{
			++nat;
			noLevels = atoi((*nat).c_str());
			if(noLevels==0) noLevels = 3;
		};

	};

	NoDataNode * aNode = new NoDataNode(newNodeName, noDataData);
	
	aNode->setSimNodeType(simNodeType);
	aNode->setSimNoLevels(noLevels);
	aNode->setIsSNPNode(isSNP);
	noDataNetworkNodes[newNodeName] = aNode;
	if(!allNodeData->isNodeDataAdded(noDataData->name)) allNodeData->addNodeData(noDataData);
	unsigned int nodeNo = allNodeData->getNodeDataNumber(newNodeName);
	aNode->setNodeID(nodeNo);
	allNetworkNodes[nodeNo] = aNode;
};

//! Adds node to a network.
void Network::addNodeInit(const string & nodeName)
{
	list<string> nodeAndType = splitString(nodeName, "|");
	string theNodeName;
	if(!nodeAndType.empty()) theNodeName = *nodeAndType.begin();
	else exitErr("Problem adding node "+nodeName+"!");

	unsigned nodeNumber = allNodeData->getNodeDataNumberInit(theNodeName);

	if(nodeNumber == 0)
	{ 
		NoData * noDataData = new NoData(theNodeName);
		addNoDataNode(nodeName, noDataData);
	}
	else
	{
		addNode(nodeNumber);
	};

	//add other factors variables if this node is a head factor node
	list<unsigned int> otherFactors = allNodeData->getFactorDataGroup(nodeNumber);

	for(list<unsigned int>::const_iterator of = otherFactors.begin(); of != otherFactors.end(); ++of)
	{
		addNode(*of);
	};

};

//! Adds a black node to the network.
void Network::addBlackNode(const string & nodeName)
{
	unsigned int nodeNumber = allNodeData->getNodeDataNumber(nodeName);
	blackNodes.insert(nodeNumber);
};

//! Adds a black edge to the network.
void Network::addBlackEdge(const string & nodeName1, const string & nodeName2)
{
	unsigned int nodeNumber1 = allNodeData->getNodeDataNumber(nodeName1);
	unsigned int nodeNumber2 = allNodeData->getNodeDataNumber(nodeName2);
	blackEdges.insert(make_pair(nodeNumber1, nodeNumber2));
};

//! Sets a node to not allow any parents.
void Network::addNoParentsNode(const string & nodeName)
{
	unsigned int nodeNumber = allNodeData->getNodeDataNumber(nodeName);
	noParentsNodes.insert(nodeNumber);
};

//! Sets all SNP nodes to not allow parent nodes.
void Network::setSNPNodesAsNoParentsNodes()
{
	for(map<unsigned int, Node *>::const_iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		if(nd->second->getIsSNPNode()) addNoParentsNode(nd->second->getName());
	};
};

//! Sets node to not allow children.
void Network::addNoChildrenNode(const string & nodeName)
{
	unsigned int nodeNumber = allNodeData->getNodeDataNumber(nodeName);
	noChildrenNodes.insert(nodeNumber);
};

//! Update the cost of edges.
void Network::setCostEdges(map<pair<unsigned int, unsigned int>, double> & costEdges0, map<pair<unsigned int, unsigned int>, double> & costEdgeTypes0)
{
	if(freeClearMemory)
	{
		map<pair<unsigned int, unsigned int>, double>().swap(costEdges);
		map<pair<unsigned int, unsigned int>, double>().swap(costEdgeTypes);
	}
	else
	{
		costEdges.clear();
		costEdgeTypes.clear();
	};

	for(map<pair<unsigned int, unsigned int>, double>::const_iterator ce = costEdges0.begin(); ce != costEdges0.end(); ++ce)
	{
		costEdges[ce->first] = ce->second;
	};

	for(map<pair<unsigned int, unsigned int>, double>::const_iterator cet = costEdgeTypes0.begin(); cet != costEdgeTypes0.end(); ++cet)
	{
		costEdgeTypes[cet->first] = cet->second;
	};
};

//! Updates the black and white edges.
void Network::updateBlackWhite(Network * network)
{

	whiteNodes = network->getWhiteNodes();
	blackNodes = network->getBlackNodes();
	whiteEdges = network->getWhiteEdges();
	blackEdges = network->getBlackEdges();
	allowedEdgeTypes = network->getAllowedEdgeTypes();
	noParentsNodes = network->getNoParentsNodes();
	noChildrenNodes = network->getNoChildrenNodes();

	//add white edges to network
	for(set<pair<unsigned int, unsigned int> >::const_iterator we = whiteEdges.begin(); we != whiteEdges.end(); ++we)
	{
		addEdge(we->first, we->second);
	};
};

//! Updates the black and white edges when display name of nodes are the same but data is different.
void Network::updateBlackWhiteDifferentData(Network * network, const string & preNodeName)
{

	string nodeName1, nodeName2;
	set<unsigned int> otherWhiteNodes = network->getWhiteNodes();
	for(set<unsigned int>::const_iterator own = otherWhiteNodes.begin(); own != otherWhiteNodes.end(); ++own)
	{
		nodeName1 = preNodeName + network->getNetworkNode(*own)->getName();
		addWhiteNode(nodeName1);
	};

	set<unsigned int> otherBlackNodes = network->getBlackNodes();
	for(set<unsigned int>::const_iterator obn = otherBlackNodes.begin(); obn != otherBlackNodes.end(); ++obn)
	{
		nodeName1 = preNodeName + network->getNetworkNode(*obn)->getName();
		addBlackNode(nodeName1);
	};

	set<pair<unsigned int, unsigned int> > otherWhiteEdges = network->getWhiteEdges();
	for(set<pair<unsigned int, unsigned int> >::const_iterator owe = otherWhiteEdges.begin(); owe != otherWhiteEdges.end(); ++owe)
	{
		nodeName1 = preNodeName + network->getNetworkNode(owe->first)->getName();
		nodeName2 = preNodeName + network->getNetworkNode(owe->second)->getName();

		addWhiteEdge(nodeName1, nodeName2);
	};

	set<pair<unsigned int, unsigned int> > otherBlackEdges = network->getBlackEdges();
	for(set<pair<unsigned int, unsigned int> >::const_iterator obe = otherBlackEdges.begin(); obe != otherBlackEdges.end(); ++obe)
	{
		nodeName1 = preNodeName + network->getNetworkNode(obe->first)->getName();
		nodeName2 = preNodeName + network->getNetworkNode(obe->second)->getName();
		addBlackEdge(nodeName1, nodeName2);
	};
	
	allowedEdgeTypes = network->getAllowedEdgeTypes();

	set<unsigned int> otherNoParents = network->getNoParentsNodes();
	for(set<unsigned int>::const_iterator op = otherNoParents.begin(); op != otherNoParents.end(); ++op)
	{
		nodeName1 = preNodeName + network->getNetworkNode(*op)->getName();		
		addNoParentsNode(nodeName1);
	};
	
	set<unsigned int> otherNoChildren = network->getNoChildrenNodes();
	for(set<unsigned int>::const_iterator oc = otherNoChildren.begin(); oc != otherNoChildren.end(); ++oc)
	{
		nodeName1 = preNodeName + network->getNetworkNode(*oc)->getName();		
		addNoChildrenNode(nodeName1);
	};

	//add white edges to network
	for(set<pair<unsigned int, unsigned int> >::const_iterator we = whiteEdges.begin(); we != whiteEdges.end(); ++we)
	{
		addEdge(we->first, we->second);
	};

};

//! Return whether a node is black or not.
bool Network::isBlackNode(const unsigned int & nodeNo)
{
	return blackNodes.find(nodeNo) != blackNodes.end();
};

//! Returns if node exists given the node name.
bool Network::nodeExistsInit(const string & nodeName) const
{	
	unsigned int nodeNo = allNodeData->getNodeDataNumberInit(nodeName);
	
	if(nodeNo == 0) return noDataNetworkNodes.find(nodeName) != noDataNetworkNodes.end();

	return nodeExists(nodeNo);
};

//! Returns if edge exists in network given node IDs.
bool Network::edgeExists(const unsigned int & nodeNo1, const unsigned int & nodeNo2) const
{
	Node * node1;
	Node * node2;
	map<unsigned int, Node *>::const_iterator nd1 = allNetworkNodes.find(nodeNo1);
	if(nd1 != allNetworkNodes.end()) node1 = nd1->second; else return false;

	map<unsigned int, Node *>::const_iterator nd2 = allNetworkNodes.find(nodeNo2);
	if(nd2 != allNetworkNodes.end()) node2 = nd2->second; else return false;

	return node2->parentExists(nd1->first);
};

//! Returns if edge exists in network given node names.
bool Network::edgeExistsInit(const string & nodeName1, const string & nodeName2) const
{
	unsigned int nodeNumber1 = allNodeData->getNodeDataNumber(nodeName1);
	unsigned int nodeNumber2 = allNodeData->getNodeDataNumber(nodeName2);
	return edgeExists(nodeNumber1, nodeNumber2);
};

//! Removes a node from the network. 
void Network::removeNode(const unsigned int & nodeNo)
{
	map<unsigned int, Node *>::const_iterator nd = allNetworkNodes.find(nodeNo);

	if(nd != allNetworkNodes.end())
	{
		allNetworkNodes.erase(nodeNo);
		discreteNodes.erase(nodeNo);   //also remove from the discrete list, should be there
		ctsNodes.erase(nodeNo);        //also remove from the cts list, should be there
		delete nd->second;
	};

	//remove other factors variables if this node is a head factor node
	list<unsigned int> otherFactors = allNodeData->getFactorDataGroup(nodeNo);

	for(list<unsigned int>::const_iterator of = otherFactors.begin(); of != otherFactors.end(); ++of)
	{
		removeNode(*of);
	};
};

//! Remove edge where Node 1 is parent and Node 2 is child.
void Network::removeEdge(const unsigned int & nodeNo1, Node * node2)
{	
	node2->removeParent(nodeNo1);
	
	//if this node is a head factor node then remove edges of child factor nodes 
	list<unsigned int> otherFactors = allNodeData->getFactorDataGroup(nodeNo1);

	for(list<unsigned int>::const_iterator of = otherFactors.begin(); of != otherFactors.end(); ++of)
	{				
		node2->removeParent(*of);
	};	
};

//! Remove edge where Node 1 is parent and Node 2 is child.
void Network::removeEdge(const unsigned int & nodeNo1, const unsigned int & nodeNo2)
{
	map<unsigned int, Node *>::const_iterator nd2 = allNetworkNodes.find(nodeNo2);

	if(nd2 != allNetworkNodes.end())
	{
		map<unsigned int, Node *>::const_iterator nd1 = allNetworkNodes.find(nodeNo1);
		if(nd1 != allNetworkNodes.end())
		{
			nd2->second->removeParent(nd1->first);
	
			//if this node is a head factor node then remove edges of child factor nodes 
			list<unsigned int> otherFactors = allNodeData->getFactorDataGroup(nodeNo1);

			for(list<unsigned int>::const_iterator of = otherFactors.begin(); of != otherFactors.end(); ++of)
			{				
				nd2->second->removeParent(*of);
			};
		};
	};
};

//! Removes all edges of the network.
void Network::removeAllEdges()
{
	for(map<unsigned int, Node *>::const_iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		nd->second->removeAllParents();
	};

};

//! Deletes all network node data, use with care, for tidy up after bootstrapping.
void Network::deleteAllNetworkNodeData()
{
	for(map<unsigned int, Node *>::const_iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		allNodeData->deleteNodeData(nd->first);
	};
};

//! Reverse edge where Node 1 is parent and Node 2 is child.
void Network::reverseEdge(Node * node1, Node * node2)
{
	removeEdge(node1->getNodeID(), node2);
	addEdge(node2, node1);	
};

//! Reverse edge where Node 1 is parent and Node 2 is child.
void Network::reverseEdge(const unsigned int & nodeNo1, const unsigned int & nodeNo2)
{
	map<unsigned int, Node *>::const_iterator nd2 = allNetworkNodes.find(nodeNo2);

	if(nd2 != allNetworkNodes.end())
	{
		map<unsigned int, Node *>::const_iterator nd1 = allNetworkNodes.find(nodeNo1);
		if(nd1 != allNetworkNodes.end())
		{		
			removeEdge(nd1->first, nd2->second);
			addEdge(nd2->second, nd1->second);				
		};
	};

};

//! Returns the number of edges and nodes of the network.
void Network::getNumberNodesAndEdges(unsigned int & noDisNodes, unsigned int & noFactorNodes, unsigned int & noCtsNodes, unsigned int & noNoNodes, unsigned int & noEdges)
{
	noDisNodes = discreteNodes.size();
	noCtsNodes = ctsNodes.size();
	noNoNodes = noDataNetworkNodes.size();
	noFactorNodes = 0; 
	unsigned int noFactorChildNodes = 0; 

	noEdges = 0;
	for(map<unsigned int, Node *>::const_iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		noEdges += nd->second->getNoParents();
		noFactorNodes += nd->second->getIsFactorNode();
		noFactorChildNodes += nd->second->getIsFactorChildNode();
	};

	noCtsNodes = noCtsNodes - noFactorNodes - noFactorChildNodes;
};

//! Gets network missing data ref.
NetworkMissingData * Network::getNetworkMissingData()
{
	if(allNetworkNodes.empty()) exitErr("Network has no nodes!");

	return allNetworkNodes.begin()->second->getNetworkMissingData();
};

//! Sets ref for each node for which data is missing for this network.
void Network::setNetworkMissingData(Network * network)
{
	NetworkMissingData * networkMissingData = network->getNetworkMissingData();

	for(map<unsigned int, Node *>::const_iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		 nd->second->setNetworkMissingData(networkMissingData);
	};

};


//! Sets up network, add white nodes and edges and removes the black ones.
void Network::initialise()
{

	set<unsigned int>::const_iterator wn;
	set<unsigned int>::const_iterator bn;

	//add white nodes
	for(wn = whiteNodes.begin(); wn != whiteNodes.end(); ++wn)
	{
		bn = blackNodes.find(*wn);

		if(bn != blackNodes.end())
		{
			outErr("Node "); outErr(allNodeData->getNodeDataName(*wn)); outErr(" is in the whitelist AND the blacklist!");
			exitErr("");
		};

		if(!nodeExists(*wn))
		{
			string nodeName = allNodeData->getNodeDataName(*wn);
			addNodeInit(nodeName);
		};
	};

	//remove black nodes
	for(bn = blackNodes.begin(); bn != blackNodes.end(); ++bn)
	{
		removeNode(*bn);
	};

	set<pair<unsigned int, unsigned int> >::const_iterator we;
	set<pair<unsigned int, unsigned int> >::const_iterator be;

	//add white edges
	for(we = whiteEdges.begin(); we != whiteEdges.end(); ++we)
	{
		be = blackEdges.find(*we);
		if(be != blackEdges.end())
		{
			string mess = "Edge " + allNodeData->getNodeDataName(we->first) + " --> " + allNodeData->getNodeDataName(we->second) + " is in the whitelist AND the blacklist!";
			exitErr(mess);
		};

		//check other direction does not already exist
		if(!(edgeExists(we->second, we->first) && isWhiteEdge(we->second, we->first)))
		{
			removeEdge(we->first, we->second);
			addEdge(we->first, we->second);
		};
	};

	//remove black edges
	for(be = blackEdges.begin(); be != blackEdges.end(); ++be)
	{
		removeEdge(be->first, be->second);
	};

	//check the edges are valid
	map<unsigned int, Node *> parents;

	for(map<unsigned int, Node *>::const_iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{

		parents = nd->second->getParents();
		for(map<unsigned int, Node *>::const_iterator p = parents.begin(); p != parents.end(); ++p)
		{

			if(!edgeAllowed(p->first, nd->first) && !p->second->getIsFactorChildNode()) //excuse a child factor --> some node
			{				
				outErr("Edge "); outErr(p->second->getDisplayName()); outErr(" --> "); outErr(nd->second->getDisplayName()); outErr(" is invalid!"); outErr("\n");
				exitErr("");
			};

		};

	};

};

//! Returns the data size of the network.
unsigned int Network::getDataSize() const
{
	if(!discreteNodes.empty())
	{
		return discreteNodes.begin()->second->getAmountData();
	}
	else if(!ctsNodes.empty())
	{
		return ctsNodes.begin()->second->getAmountData();
	}
	else if(!noDataNetworkNodes.empty())
	{
		return 0;
	}
	else
	{
		exitErr("The network has no nodes!\n");
	};
	
	return 0;
};

//! Updates missing data for a network.
NetworkMissingData * Network::updateNetworkMissingData()
{
	//check for undefined nodes
	if(!noDataNetworkNodes.empty())
	{
		outErr("Node "); outErr(noDataNetworkNodes.begin()->second->getDisplayName()); outErr(" has no data!\n");
		outErr("The network requires data for every node except when simulating data, or maybe calculating the Markov blanket.");
		exitErr("");
	};

	NetworkMissingData * networkMissingData = new NetworkMissingData();

	unsigned int sizeData = getDataSize(); //uses first node but will be checked along the way
	
	//set up initial missing data as no missing data, then update as missing if any found
	for(unsigned int i = 1; i <= sizeData; ++i) networkMissingData->missing.push_back(false);

	//loop thro' discrete node data
	for(map<unsigned int, DiscreteNode *>::iterator nd = discreteNodes.begin(); nd != discreteNodes.end(); ++nd)
	{
		nd->second->updateNetworkMissingData(networkMissingData);
		nd->second->setNetworkMissingData(networkMissingData);
	};

	//loop thro cts node data
	for(map<unsigned int, CtsNode *>::iterator nd = ctsNodes.begin(); nd != ctsNodes.end(); ++nd)
	{
		nd->second->updateNetworkMissingData(networkMissingData);
		nd->second->setNetworkMissingData(networkMissingData);
	};
	
	return networkMissingData;
};

//! Updates missing data for a network when sim'ing data.
NetworkMissingData * Network::updateNetworkMissingDataSimNet(const unsigned int & noIndivs)
{
	NetworkMissingData * networkMissingData = new NetworkMissingData();

	//set up initial missing data as no missing data, then update as missing if any found
	for(unsigned int i = 1; i <= noIndivs; ++i) networkMissingData->missing.push_back(false);

	//loop thro' discrete node data
	for(map<unsigned int, DiscreteNode *>::iterator nd = discreteNodes.begin(); nd != discreteNodes.end(); ++nd)
	{		
		nd->second->setNetworkMissingData(networkMissingData);
	};

	//loop thro cts node data
	for(map<unsigned int, CtsNode *>::iterator nd = ctsNodes.begin(); nd != ctsNodes.end(); ++nd)
	{		
		nd->second->setNetworkMissingData(networkMissingData);
	};

	return networkMissingData;
};

//! Returns amount of missing and not missing data for the network.
pair<unsigned int, unsigned int> Network::getNoMissingNotMissing()
{
	if(allNetworkNodes.size() > 0)
	{
		return allNetworkNodes.begin()->second->getNoMissingNotMissing();
	};

	return make_pair(0, 0);
};

//! Sets default parameters used for sim'ing data.
void Network::setDefaultParameters()
{
	for(map<unsigned int, Node *>::const_iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		nd->second->setDefaultParameters();
	};

};

//! Sets up uniform local prior probs.
void Network::setupDefaultInitPriors()
{
	//loop thro' parameters and set all local pri probs to uniform dists
	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		nd->second->setDefaultInitPrior(); //if discrete set uniform over no of levels, if cts fit linear regression with parents
	};

};

//! Returns Markov Blanket for a given node as a network.
set<Node *> Network::getMarkovBlanket(const string & nodeName)
{
	
	set<Node *> markovBlanket;
	unsigned int nodeNo = allNodeData->getNodeDataNumber(nodeName);
	Node * theNode = getNetworkNode(nodeNo);

	//add parents
	map<unsigned int, Node *> parents = theNode->getParents();
	for(map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end(); ++pa)
	{
		markovBlanket.insert(pa->second);
	};

	list<Node *> children;

	//add children
	for(map<unsigned int, Node *>::const_iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		if(nd->second->parentExists(nodeNo))
		{
			markovBlanket.insert(nd->second);
			children.push_back(nd->second);
		};
	};

	map<unsigned int, Node *> parentsCld;

	//add spouses
	for(list<Node *>::const_iterator cld = children.begin(); cld != children.end(); ++cld)
	{
		parentsCld = (*cld)->getParents();
		for(map<unsigned int, Node *>::const_iterator pac = parents.begin(); pac != parents.end(); ++pac)
		{
			if(pac->first != nodeNo) markovBlanket.insert(pac->second);
		};

	};

	return markovBlanket;
};

//! Returns Markov blanket for a node.
Network * Network::getMarkovBlanketSubNetwork(const string & nodeName)
{
	NetworkBNLearn * markovBlanketNetwork = new NetworkBNLearn(allNodeData, getScoreType(), getScoreFix());

	unsigned int nodeNo = allNodeData->getNodeDataNumber(nodeName);

	Node * theNode = getNetworkNode(nodeNo);

	markovBlanketNetwork->addNodeInit(nodeName);

	Node * theNodeInBlanket = markovBlanketNetwork->getNetworkNode(nodeNo);
	Node * nodeInBlanket; 
	string nodeNameInit;

	//add parents
	map<unsigned int, Node *> parents = theNode->getParents();
	for(map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end(); ++pa)
	{
		nodeNameInit = pa->second->getName();
		markovBlanketNetwork->addNodeInit(nodeNameInit);
		nodeInBlanket = markovBlanketNetwork->getNetworkNode(pa->first); //parent
		addEdge(nodeInBlanket, theNodeInBlanket);
	};

	list<Node *> children;

	//add children
	for(map<unsigned int, Node *>::const_iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		if(nd->second->parentExists(nodeNo))
		{
			markovBlanketNetwork->addNodeInit(nd->second->getName());
			children.push_back(nd->second);
		};
	};

	map<unsigned int, Node *> parentsCld;
	Node * nodeInBlanket2;
	unsigned int nodeID;

	//add spouses
	for(list<Node *>::const_iterator cld = children.begin(); cld != children.end(); ++cld)
	{
		nodeID = (*cld)->getNodeID();
		nodeInBlanket = markovBlanketNetwork->getNetworkNode(nodeID); //child
		parentsCld = (*cld)->getParents();

		for(map<unsigned int, Node *>::const_iterator pac = parentsCld.begin(); pac != parentsCld.end(); ++pac)
		{
			if(!markovBlanketNetwork->nodeExists(pac->first)) markovBlanketNetwork->addNodeInit(pac->second->getName());
			nodeInBlanket2 = markovBlanketNetwork->getNetworkNode(pac->first); //spouse
			addEdge(nodeInBlanket2, nodeInBlanket);
		};

	};

	return markovBlanketNetwork;
};

//! Sets up local prior probs.
void NetworkDeal::setupLocalPriProbs(const string & localProbFilename)
{
	//set up using given user prior, or otherwise use uniform priors
	if(localProbFilename == "") {setupDefaultInitPriors(); return;};

	//setup local priors using priors given in file

};

//! Returns nodes group number for a certain discrete node level combo.
unsigned int Network::getNodesGroupNo()
{
	unsigned int groupNo = 1;
	unsigned int prevLevels = 1;

	map<unsigned int, DiscreteNode *>::const_iterator nd = discreteNodes.begin();
	
	while(nd != discreteNodes.end())
	{
		groupNo += (nd->second->getLevel() - 1)*prevLevels;
		
		prevLevels *= (nd->second->getNoLevels());

		++nd;
	};

	return groupNo;
};

//! Calculates the joint prior prob for current level states of Nodes.
double NetworkDeal::calcDiscreteJointProbDist()
{
	double prob = 1;

	//loop thro nodes and set probs if all parents are given, if not skip
	for(map<unsigned int, DiscreteNode *>::iterator nd = discreteNodes.begin(); nd != discreteNodes.end(); ++nd)
	{
		prob *= nd->second->getPriLocalProb();
	};

	return prob;
};

//! Calculates discrete nodes current group level name.
string NetworkDeal::calcDiscreteJointName()
{
	string name = "";

	//loop thro nodes and set probs if all parents are given, if not skip
	for(map<unsigned int, DiscreteNode *>::iterator nd = discreteNodes.begin(); nd != discreteNodes.end(); )
	{
		name += nd->second->getLevelName();  
		++nd;
		if(nd != discreteNodes.end()) name += ":";
	};

	return name;
};

//! Returns discrete joint prior prob, p(i); alpha_i = p(i) * ISS.
double NetworkDeal::getDiscreteJointPriProb(unsigned int & nodeGroupNo)
{
	map<unsigned int, double>::const_iterator djpp = discreteJointProbDist.find(nodeGroupNo);

	if(djpp != discreteJointProbDist.end()) return djpp->second;

	out("Problem setting up master prior, unable to find node group "); out(nodeGroupNo); out("\n");
	exitErr("");

	return 0;
};

//Sets up the discrete joint prior probs.
void NetworkDeal::calcDiscreteJointPriorDist()
{
	//loop thro' all possible combos for the Nodes, like mileometer
	unsigned int totalNodes = discreteNodes.size();
	unsigned int nodeGroupNo = 1;
	map<unsigned int, DiscreteNode *>::iterator aNode = discreteNodes.begin();
	unsigned int levelNo = 1;
	
	//set up initial levels
	for(map<unsigned int, DiscreteNode *>::iterator nd = discreteNodes.begin(); nd != discreteNodes.end(); ++nd)
	{
		nd->second->setLevel(1);
	};

	bool updated = false;

	do{

		//get prob
		discreteJointProbDist[nodeGroupNo] = calcDiscreteJointProbDist();

		updated = false;
		for(map<unsigned int, DiscreteNode *>::iterator ndc = discreteNodes.begin(); ndc != discreteNodes.end(); ++ndc)
		{
			if(ndc->second->getLevel() < ndc->second->getNoLevels())
			{
				ndc->second->setLevel(ndc->second->getLevel()+1);
				//set prev nodes to 1
				updated = true;
				break;
			}
			else
			{
				ndc->second->setLevel(1);
			};

		};

		if(!updated) break;

		nodeGroupNo = getNodesGroupNo();
	}while(true);

};

//! Sets up prior and posterior names, for output.
void NetworkDeal::setPriorPostNames()
{
	calcDiscreteJointPriorNames();
	calcCtsJointPriorNames();

	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		nd->second->calcDiscreteNodeParentNames();
		nd->second->calcLocalParaNames();
	};
};

//! Sets up prior and posterior names, for output.
void NetworkBNLearn::setPriorPostNames()
{
	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		nd->second->calcDiscreteNodeParentNames();
		nd->second->calcLocalParaNames();
	};
};

//Sets up the discrete joint prior names.
void NetworkDeal::calcDiscreteJointPriorNames()
{
	if(freeClearMemory) map<unsigned int, string>().swap(discreteJointNames); else discreteJointNames.clear();

	//loop thro' all possible combos for the Nodes, like mileometer
	unsigned int totalNodes = discreteNodes.size();
	unsigned int nodeGroupNo = 1;
	map<unsigned int, DiscreteNode *>::iterator aNode = discreteNodes.begin();
	unsigned int levelNo = 1;
	
	//set up initial levels
	for(map<unsigned int, DiscreteNode *>::iterator nd = discreteNodes.begin(); nd != discreteNodes.end(); ++nd)
	{
		nd->second->setLevel(1);
	};

	bool updated = false;

	do{

		//get prob		
		discreteJointNames[nodeGroupNo] = calcDiscreteJointName();

		updated = false;
		for(map<unsigned int, DiscreteNode *>::iterator ndc = discreteNodes.begin(); ndc != discreteNodes.end(); ++ndc)
		{
			if(ndc->second->getLevel() < ndc->second->getNoLevels())
			{
				ndc->second->setLevel(ndc->second->getLevel()+1);
				//set prev nodes to 1
				updated = true;
				break;
			}
			else
			{
				ndc->second->setLevel(1);
			};

		};

		if(!updated) break;

		nodeGroupNo = getNodesGroupNo();
	}while(true);

};

//! Returns the Joint Prob dist for cts nodes for given discrete node settings, or create a new one.
CtsMultiDistJoint * NetworkDeal::getCtsJointProbDist(unsigned int & allNodesGrpNo)
{
	map<unsigned int, CtsMultiDistJoint *>::const_iterator cjd = ctsJointProbDist.find(allNodesGrpNo);

	if(cjd != ctsJointProbDist.end()) return cjd->second;

	CtsMultiDistJoint * aCtsJointDist = new CtsMultiDistJoint();
	ctsJointProbDist[allNodesGrpNo] = aCtsJointDist;
	return aCtsJointDist;
};

//! Sets up Joint priors for the cts nodes.
void NetworkDeal::calcCtsJointPriorDist()
{
	
	//loop thro' the nodes in order of dependency
	//set all nodes to unvisited
	clearVisitedNodes();

	unsigned int noVisitedThisLoop = 0;
	unsigned int totalVisited = 0;
	unsigned int totalNodes = allNetworkNodes.size();

	bool updated = false;
	unsigned int nodeGroupNo;
	map<unsigned int, DiscreteNode *>::iterator aNode;
	unsigned int levelNo;
	CtsMultiDistJoint * ctsJointDist;

	//try to visit every node if the parents are given (if any), if no loops then can visit every loop
	//output even if a string is detected, but try to do in order first
	do{
		noVisitedThisLoop = 0;

		for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
		{
			if(!nd->second->getVisited() && nd->second->allParentsVisited())
			{
				//loop thro' every discrete combo
				//loop thro' all possible combos for the Nodes, like mileometer
				
				nodeGroupNo = 1;
				aNode = discreteNodes.begin();
				levelNo = 1;
	
				//set up initial levels
				for(map<unsigned int, DiscreteNode *>::iterator ndd = discreteNodes.begin(); ndd != discreteNodes.end(); ++ndd)
				{
					ndd->second->setLevel(1);
				};

				updated = false;

				do{
						
					//update setup cts node joint dist
					//get the cts joint prob dist for this discrete node group
					ctsJointDist = getCtsJointProbDist(nodeGroupNo);

					//update the cts joint prob dist with mean, variance and covariances for the current node
					nd->second->updateCtsJointProbDist(ctsJointDist);

					updated = false;
					for(map<unsigned int, DiscreteNode *>::iterator ndc = discreteNodes.begin(); ndc != discreteNodes.end(); ++ndc)
					{
						if(ndc->second->getLevel() < ndc->second->getNoLevels())
						{
							ndc->second->setLevel(ndc->second->getLevel()+1);
							//set prev nodes to 1
							updated = true;
							break;
						}
						else
						{
							ndc->second->setLevel(1);
						};

					};

					if(!updated) break;

					nodeGroupNo = getNodesGroupNo();
				}while(true); //end of looping thro' discrete node combos

				nd->second->setVisited(true);
				++noVisitedThisLoop;
			};

		};//end of looping thro' nodes

		if(noVisitedThisLoop == 0)
		{
			exitErr("Loop detected in network while calculating the continuous joint probability network!\n");
		};

		totalVisited += noVisitedThisLoop;
	}while(totalVisited < totalNodes);

};

//! Sets up Joint priors names for the cts nodes.
void NetworkDeal::calcCtsJointPriorNames()
{
	
	if(freeClearMemory) map<unsigned int, string>().swap(meanNames); else meanNames.clear();

	//loop thro' the nodes in order of dependency
	//set all nodes to unvisited
	clearVisitedNodes();

	unsigned int noVisitedThisLoop = 0;
	unsigned int totalVisited = 0;
	unsigned int totalNodes = allNetworkNodes.size();

	bool updated = false;
	unsigned int nodeGroupNo;
	map<unsigned int, DiscreteNode *>::iterator aNode;
	unsigned int levelNo;
	CtsMultiDistJoint * ctsJointDist;

	//try to visit every node if the parents are given (if any), if no loops then can visit every loop
	//output even if a string is detected, but try to do in order first
	do{
		noVisitedThisLoop = 0;

		for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
		{
			if(!nd->second->getVisited() && nd->second->allParentsVisited())
			{
				//loop thro' every discrete combo
				//loop thro' all possible combos for the Nodes, like mileometer	
				nodeGroupNo = 1;
				aNode = discreteNodes.begin();
				levelNo = 1;
	
				//set up initial levels
				for(map<unsigned int, DiscreteNode *>::iterator ndd = discreteNodes.begin(); ndd != discreteNodes.end(); ++ndd)
				{
					ndd->second->setLevel(1);
				};

				updated = false;
				if(nd->second->getNodeType() == "c") meanNames[nd->first] = nd->second->getName();

				do{
			
					//get the cts joint name for this discrete node group
					ctsJointDist = getCtsJointProbDist(nodeGroupNo);
					if(ctsJointDist->name == "") ctsJointDist->name = calcDiscreteJointName();
						
					updated = false;
					for(map<unsigned int, DiscreteNode *>::iterator ndc = discreteNodes.begin(); ndc != discreteNodes.end(); ++ndc)
					{
						if(ndc->second->getLevel() < ndc->second->getNoLevels())
						{
							ndc->second->setLevel(ndc->second->getLevel()+1);
							//set prev nodes to 1
							updated = true;
							break;
						}
						else
						{
							ndc->second->setLevel(1);
						};

					};

					if(!updated) break;

					nodeGroupNo = getNodesGroupNo();
				}while(true); //end of looping thro' discrete node combos

				nd->second->setVisited(true);
				++noVisitedThisLoop;
			};

		};//end of looping thro' nodes

		if(noVisitedThisLoop == 0)
		{			
			exitErr("Loop detected in network while calculating the continuous joint probability network!\n");
		};

		totalVisited += noVisitedThisLoop;
	}while(totalVisited < totalNodes);

};

//! Calculates the cts local parameter priors.
void NetworkDeal::calcCtsLocalParaPris()
{
	//loop thro' nodes and calculate the Local Parameter priors
	for(set<CtsDealNode *>::iterator nd = ctsDealNodes.begin(); nd != ctsDealNodes.end(); ++nd)
	{
		(*nd)->calcCtsLocalParameterPrior(this);
	};

};

//! Calculates the cts local parameter posteriors.
void NetworkDeal::calcCtsLocalParaPosts()
{
	//loop thro' nodes and calculate the Local Parameter posteriors
	for(set<CtsDealNode *>::iterator nd = ctsDealNodes.begin(); nd != ctsDealNodes.end(); ++nd)
	{
		(*nd)->calcCtsLocalParameterPost(this);
	};

};

//! Calculates the significance of edges.
void Network::calcEdgeSignifs()
{
	unsigned int oldType = getScoreType();
	setScoreType(0); //log like

	//loop thro' nodes and calculate the edge signifances
	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		nd->second->calcEdgeSignifs(this);
	};

	setScoreType(oldType);
};

//! Gets edge probability.
double Network::getEdgeProb(const string & aParent, const string & nodeStr)
{
	unsigned int nodeNumber1 = allNodeData->getNodeDataNumber(nodeStr);
	unsigned int nodeNumber2 = allNodeData->getNodeDataNumber(aParent);

	Node * aNode1 = getNetworkNode(nodeNumber1);
	Node * aNode2 = getNetworkNode(nodeNumber2);

	return aNode1->getEdgeProb(nodeNumber2);
};

//! Gets number of parents of node.
unsigned int Network::getNumberOfParents(const string & nodeStr)
{
	unsigned int nodeNumber1 = allNodeData->getNodeDataNumber(nodeStr);
	Node * aNode1 = getNetworkNode(nodeNumber1);

	unsigned int numParents = aNode1->getNoParents();

	return numParents;
};

//! Gets number of children of node.
unsigned int Network::getNumberOfChildren(const string & nodeStr)
{
	unsigned int nodeNumber1 = allNodeData->getNodeDataNumber(nodeStr);
	
	map<unsigned int, Node *> parents;
	unsigned int numOfChildren = 0;

	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		parents = nd->second->getParents();
		for(map<unsigned int, Node *>::const_iterator p = parents.begin(); p != parents.end(); ++p)
		{
			if(p->first == nodeNumber1) numOfChildren++;
		};
	};

	return numOfChildren;
};


//! Gets edge signif.
pair<double, unsigned int> Network::getEdgeSign(const string & aParent, const string & nodeStr)
{
	unsigned int nodeNumber1 = allNodeData->getNodeDataNumber(nodeStr);
	unsigned int nodeNumber2 = allNodeData->getNodeDataNumber(aParent);

	Node * aNode1 = getNetworkNode(nodeNumber1);
	Node * aNode2 = getNetworkNode(nodeNumber2);

	return aNode1->getEdgeSign(nodeNumber2);
};

//! Returns number of parameters in network
unsigned int Network::getNumberOfParameters()
{
	unsigned int noParameters = 0;
	
	//loop thro' nodes and calculate the edge signifances
	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		noParameters += nd->second->getNumberOfParameters();
	};

	return noParameters;
};

//! Unfix all nodes.
void Network::setAllNodesUnfixed()
{
	for(map<unsigned int, DiscreteNode *>::iterator nd = discreteNodes.begin(); nd != discreteNodes.end(); ++nd)
	{
		nd->second->setFixed(false);
	};
};

//! Calculates a master prior parameter value given discrete node levels.
double NetworkDeal::getDiscreteMasterPriorPara()
{
	double masterPriorPara = 0;

	//loop thro' all possible combos for the Nodes, like mileometer, keeping values fixed in nodeValues unless equal to 0
	unsigned int nodeGroupNo;
	
	//set up initial levels for unfixed nodes
	for(map<unsigned int, DiscreteNode *>::iterator nd = discreteNodes.begin(); nd != discreteNodes.end(); ++nd)
	{
		if(!nd->second->getFixed()) nd->second->setLevel(1);
	};

	do{
		nodeGroupNo = getNodesGroupNo();
		
		//get prob
		masterPriorPara += getDiscreteJointPriProb(nodeGroupNo);

		bool updated = false;
		for(map<unsigned int, DiscreteNode *>::iterator ndc = discreteNodes.begin(); ndc != discreteNodes.end(); ++ndc)
		{
			//only update nodes that are not fixed
			if(!ndc->second->getFixed())
			{
				if(ndc->second->getLevel() < ndc->second->getNoLevels())
				{
					ndc->second->setLevel(ndc->second->getLevel()+1);
					//set prev nodes to 1
					updated = true;
					break;
				}
				else
				{
					ndc->second->setLevel(1);
				};
			};
		};

		if(!updated) break;
	
	}while(true);

	//scale by imaginary sample size to weight prior belief
	return masterPriorPara*imaginarySampleSize;
};

//! Calculates a master prior parameter value given discrete node levels.
CtsMultiDistMaster * NetworkDeal::getCtsMasterPriorDist(CtsNode * ctsNode)
{
	double jointPriorProbSum = 0; //discrete probs summed
	double jointPriorProb;
	CtsMultiDistJoint * jointCtsMultiDist;
	CtsMultiDistMaster * masterCtsMultiDist = new CtsMultiDistMaster();
	list<unsigned int> nodeGroupsForCtsNodeFixed;

	//set up means, node+ cts parents	
	unsigned int nodeID = ctsNode->getNodeID();
	list<unsigned int> ctsNodeAndParents;
	ctsNodeAndParents.push_back(nodeID);

	masterCtsMultiDist->means.push_back(0);
	list<unsigned int> ctsParents = ctsNode->getCtsParents();
	for(list<unsigned int>::const_iterator cnp = ctsParents.begin(); cnp != ctsParents.end(); ++cnp)
	{
		ctsNodeAndParents.push_back(*cnp);
		masterCtsMultiDist->means.push_back(0);
	};
	
	//setup zero covaraince matrix
	list<double> zeroCovariances;
	for(list<unsigned int>::const_iterator cp = ctsNodeAndParents.begin(); cp != ctsNodeAndParents.end(); ++cp)
	{
		zeroCovariances.push_back(0);
		masterCtsMultiDist->covariances.push_front(zeroCovariances);
	};

	map<unsigned int, double>::const_iterator jm; 
	list<double>::iterator mm, mm2;
	
	//loop thro' all possible combos for the Nodes, like mileometer, keeping values fixed in nodeValues unless equal to 0
	unsigned int nodeGroupNo;
	
	//set up initial levels for unfixed nodes
	for(map<unsigned int, DiscreteNode *>::iterator nd = discreteNodes.begin(); nd != discreteNodes.end(); ++nd)
	{
		if(!nd->second->getFixed()) nd->second->setLevel(1);
	};

	
	do{
		nodeGroupNo = getNodesGroupNo();
		nodeGroupsForCtsNodeFixed.push_back(nodeGroupNo); //keep group nos for calc'ing variance

		//get prob
		jointPriorProb = getDiscreteJointPriProb(nodeGroupNo); //*imaginarySampleSize; //sample size cancels out for means //alpha_i
		jointPriorProbSum += jointPriorProb;
		
		//get cts joint dist and add subset of this node to total
		jointCtsMultiDist = getCtsJointProbDist(nodeGroupNo);
	
		//loop thro' mater means and update
		mm = masterCtsMultiDist->means.begin(); 
		for(list<unsigned int>::const_iterator cp = ctsNodeAndParents.begin(); cp != ctsNodeAndParents.end(); ++cp, ++mm)
		{

			jm = jointCtsMultiDist->means.find(*cp); //find corresponding joint mean			

			if(jm != jointCtsMultiDist->means.end())
			{
				// joint prior *  weight given by discrete joint prior(alpha_j)
				*mm += jm->second*jointPriorProb;
			}
			else
			{
				exitErr("Problem setting up master prior means for continuous node!");
			};
		};

		//start updating group number
		bool updated = false;
		for(map<unsigned int, DiscreteNode *>::iterator ndc = discreteNodes.begin(); ndc != discreteNodes.end(); ++ndc)
		{
			//only update nodes that are not fixed
			if(!ndc->second->getFixed())
			{
				if(ndc->second->getLevel() < ndc->second->getNoLevels())
				{
					ndc->second->setLevel(ndc->second->getLevel()+1);
					//set prev nodes to 1
					updated = true;
					break;
				}
				else
				{
					ndc->second->setLevel(1);
				};
			};
		};

		if(!updated) break;
		//end updating group number

	}while(true);

	//divide by total
	for(list<double>::iterator mi = masterCtsMultiDist->means.begin(); mi != masterCtsMultiDist->means.end(); ++mi)
	{
		*mi /= jointPriorProbSum;
	};
	
	map<unsigned int, double>::const_iterator jv, jm2, jcov2; 
	map<unsigned int, map<unsigned int, double> >::iterator jcov;
	double diff;
	double jointPriorParaSum = 0;
	double jointPriorPara;
	double aJointCovariance = 0;
	list< list<double> >::iterator mcov;
	list<double>::iterator mcov2;
	unsigned int row = 1;
	unsigned int col = 1;

	//do variances, loop thro' all discrete nodes group nos, where node+parents are fixed
	for(list<unsigned int>::const_iterator grp = nodeGroupsForCtsNodeFixed.begin(); grp != nodeGroupsForCtsNodeFixed.end(); ++grp)
	{
		jointPriorPara = getDiscreteJointPriProb(nodeGroupNo)*imaginarySampleSize; //alpha_i
		jointPriorParaSum += jointPriorPara;

		//get cts joint dist 
		jointCtsMultiDist = getCtsJointProbDist(nodeGroupNo);
		
		row = 1;
		mm = masterCtsMultiDist->means.begin();		
		mcov = masterCtsMultiDist->covariances.begin();

		for(list<unsigned int>::const_iterator cnp1 = ctsNodeAndParents.begin(); cnp1 != ctsNodeAndParents.end(); ++cnp1, ++mm, ++mcov, ++row)
		{
			//find corresponding row in joint prior matrix			
			jm = jointCtsMultiDist->means.find(*cnp1);
			if(jm != jointCtsMultiDist->means.end())
			{
				mm2 = masterCtsMultiDist->means.begin();
				mcov2 = (*mcov).begin();
				col = 1;

				for(list<unsigned int>::const_iterator cnp2 = ctsNodeAndParents.begin(); cnp2 != ctsNodeAndParents.end(); ++cnp2, ++mm2, ++col)
				{
				
					if(col >= row) //only calc "half" of covariance matrix as symetric and fill in later
					{
						jm2 = jointCtsMultiDist->means.find(*cnp2);
						aJointCovariance = jointCtsMultiDist->getCovariance(*cnp1, *cnp2);

						if(jm2 != jointCtsMultiDist->means.end() && mcov2 != (*mcov).end())
						{
							diff = (jm->second - *mm)*(jm2->second - *mm2);
							*mcov2 += aJointCovariance*(jointPriorPara - 1) + jointPriorPara*diff*diff;
						}
						else
						{
							exitErr("Problem setting up master prior covariance for a continuous node!");
						};

						++mcov2;
					};
				};

			}
			else
			{
				exitErr("Problem setting up master prior covariance for continuous node!");
			};

		}; 
	

	};//end of discrete node group loop

	masterCtsMultiDist->jointDiscretePriorParaSum = jointPriorParaSum;

	list< list<double> >::reverse_iterator rmcovRowToAdd;
	list<double>::reverse_iterator rmcov2;
	list<double>::reverse_iterator rmcov2prev;

	//fill in other half of the covariance matrix and scale
	for(list< list<double> >::reverse_iterator rmcov = masterCtsMultiDist->covariances.rbegin(); rmcov != masterCtsMultiDist->covariances.rend(); ++rmcov)
	{
		rmcov2 = rmcov->rbegin();
		if(rmcov2 == rmcov->rend()) exitErr("Problem constructing master prior covariance matrix!");
		rmcov2prev = rmcov2; ++rmcov2;
		
		rmcovRowToAdd = masterCtsMultiDist->covariances.rbegin();

		while(rmcov2 != rmcov->rend())
		{
			if(rmcovRowToAdd == masterCtsMultiDist->covariances.rend()) exitErr("Problem constructing a master prior covariance matrix!");

			rmcovRowToAdd->push_front(*rmcov2prev);

			++rmcovRowToAdd;
			rmcov2prev = rmcov2; ++rmcov2;
		};

	};

	return masterCtsMultiDist;
};

//! Calculate score for this network (DAG).
double Network::calcScore()
{
	if(!checkNetworkIsValid())
	{
		string mess = "Attempt to evaluate an invalid network!\n" + getNetworkString() + "\n";	
		exitErr(mess);
	};

	double score = 0;
	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{	
		score += nd->second->calcUpdatedScoreBit(this);
	};

	return score;
};

//! Calculate score for this network (DAG).
double NetworkDeal::calcScore()
{
	if(!checkNetworkIsValid())
	{
		string mess = "Attempt to evaluate an invalid network!\n" + getNetworkString() + "\n";	
		exitErr(mess);
	};

	calculatePriorsAndPosteriors(); //seems to be necessary for deal to always recalculate all priors and posteriors...

	double score = 0;
	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{	
		score += nd->second->calcUpdatedScoreBit(this);
	};
	return score;
};

//! Determine whether there is a loop or not in the network.
bool Network::hasLoop()
{
	//set all nodes to unvisited
	clearVisitedNodes();

	unsigned int noVisitedThisLoop = 0;
	unsigned int totalVisited = 0;
	unsigned int totalNodes = allNetworkNodes.size();

	//try to visit every node if the parents are given (if any), if no loops then can visit every loop
	do{
		noVisitedThisLoop = 0;

		for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
		{
			if(!nd->second->getVisited() && nd->second->allParentsVisited())
			{
				nd->second->setVisited(true);
				++noVisitedThisLoop;
			};
		};

		totalVisited += noVisitedThisLoop;
	}while(totalVisited < totalNodes && noVisitedThisLoop > 0);

	
	//should has visited all nodes if no loops
	return totalVisited < totalNodes;
};

//Display loop for debugging purposes
void Network::displayLoop()
{	
	out("Loop:\n");
	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		if(!nd->second->getVisited())
		{
			out(nd->second->getCondName()); out("\n");
		};
	};	
};

//! Sets all nodes to unvisited to loop thro' in order of dependancy.
void Network::clearVisitedNodes()
{
	//set all nodes to unvisited
	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		nd->second->setVisited(false);
	};
};

//! Sets initial data iterators for looping thro' data.
void Network::setInitialDataIterators2()
{
	//set all node iterators
	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		nd->second->setInitialDataValue2();
	};

};

//! Moves data iterators to next pt.
void Network::advanceDataIterators2()
{
	//move iterators on
	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		nd->second->nextDataValueIter2();
	};

};

//! Sets initial data iterators for looping thro' data.
void Network::setInitialDataIterators3()
{
	//set all node iterators
	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		nd->second->setInitialDataValue3();
	};

};

//! Moves data iterators to next pt.
void Network::advanceDataIterators3()
{
	//move iterators on
	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		nd->second->nextDataValueIter3();
	};

};

//! Sets estimates of the standard devs of node data.
void Network::setStandardDevs()
{
	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		nd->second->setStDev();
	};
};

//! Returns other nodes connect to this one, either children or parents.
list<unsigned int> Network::getConnectedNodes(const unsigned int & nodeID)
{
	list<unsigned int> connectedNodes;
	Node * aNode = getNetworkNode(nodeID);

	//Add children
	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		if(nd->second->parentExists(nodeID)) connectedNodes.push_back(nd->first);
	};

	//Add parents
	map<unsigned int, Node *> theParents = aNode->getParents();
	
	for(map<unsigned int, Node *>::const_iterator pa = theParents.begin(); pa != theParents.end(); ++pa)
	{
		connectedNodes.push_back(pa->first);
	};

	return connectedNodes;
};

//! Returns the number of edges not connected to a missing node of the network.
void Network::getNoEdgesNotMissing(const map<unsigned int, set<unsigned int> > & groupNodesWithMissingData, Network * network, unsigned int & noEdgesNotMissing, unsigned int & totalEdges, unsigned int & noNonMissingSingltonNodes, unsigned int & totalSingltonNodes)
{
	set<unsigned int> missingNodes;
	set<unsigned int>::const_iterator mn, mnpa;
	unsigned int noFactorChildNodes = 0; 

	map<unsigned int, Node *> someParents;

	for(map<unsigned int, set<unsigned int> >::const_iterator i = groupNodesWithMissingData.begin(); i != groupNodesWithMissingData.end(); ++i)
	{
		for(set<unsigned int>::const_iterator j = i->second.begin(); j != i->second.end(); ++j)
		{
			missingNodes.insert(convertID(*j)); //convert to bootstrap network ID
		};
	};

	set<unsigned int> nodesWithEdges;
	set<unsigned int>::const_iterator nwe;
	noNonMissingSingltonNodes = 0;
	totalSingltonNodes = 0;
	noEdgesNotMissing = 0;
	totalEdges = 0;
	for(map<unsigned int, Node *>::const_iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{	
		someParents = nd->second->getParents();
		totalEdges += someParents.size();

		
		if(someParents.size() != 0) 
		{
			nodesWithEdges.insert(nd->first);
			for(map<unsigned int, Node *>::const_iterator pa0 = someParents.begin(); pa0 != someParents.end(); ++pa0)
			{
				nodesWithEdges.insert(pa0->first);
			};
		};
		
		mn = missingNodes.find(nd->first);
		if(mn == missingNodes.end()) //not missing
		{
			for(map<unsigned int, Node *>::const_iterator pa = someParents.begin(); pa != someParents.end(); ++pa)
			{
				 mnpa = missingNodes.find(pa->first);
				 if(mnpa == missingNodes.end()) //not missing
				 {
					 noEdgesNotMissing++;
				 };
			};
		}
		
	};

	//count singleton nodes
	for(map<unsigned int, Node *>::const_iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		nwe = nodesWithEdges.find(nd->first);
		if(nwe == nodesWithEdges.end()) //singlton node
		{
			totalSingltonNodes++;
			mnpa = missingNodes.find(nd->first);
			if(mnpa == missingNodes.end()) //not missing
			{
				noNonMissingSingltonNodes++;
			};
		};
	};
};

//! Sets ID map.
void Network::setMapID(const unsigned int & otherNetID, const unsigned int & thisNetID)
{
	nodeIDmap[otherNetID] = thisNetID;
};

//! Sets up 1-to-1 map between nodes in network and bootstrap network nodes.
void Network::setNodeMap(const string & preNodeName, Network * bootstrapNetwork)
{
	nodeIDmap.clear();
	unsigned int nodeID;
	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		nodeID = allNodeData->getNodeDataNumberInit(preNodeName + nd->second->getName());
	
		setMapID(nodeID, nd->first);
		bootstrapNetwork->setMapID(nd->first, nodeID);
	};

};

//! Returns ID of corresponding node in another network.
unsigned int Network::convertID(const unsigned int & nodeID)
{
	map<unsigned int, unsigned int>::const_iterator i = nodeIDmap.find(nodeID);

	if(i != nodeIDmap.end())
	{
		return i->second;
	};

	string msg = "Cannot convert ID for node " + toString(nodeID) + "!";
	exitErr(msg);

	return 0;
};

//! Randomly replaces missing values with randomly sampled data.
void Network::randomlyFillMissingValues()
{
	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		nd->second->randomlyFillMissingValues();
	};
};

//! Sets up was imputed list.
void Network::setupWasImputed()
{
	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		nd->second->setupWasImputed();
	};
};

//! Sets data as non missing if it was imputed.
void Network::updateImputedDataAsNonMissing()
{
	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		nd->second->updateImputedDataAsNonMissing();
	};
};

//! Returns node that have missing data for the current indiv given by iterators.
list<unsigned int> Network::getNodesWithMissingData()
{
	list<unsigned int> nodes;

	//test if nodes are missing here
	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		if(nd->second->nodeDataIsMissing2()) nodes.push_back(nd->first);
	};

	return nodes;
}

//! Imputes missing data using NN for one indiv, called from original data network.
unsigned int Network::imputeDataNN(const map<unsigned int, set<unsigned int> > & groupNodesWithMissingData, const map<unsigned int, set<unsigned int> > & groupNodesWithCompleteData, Network * bootstrapNetwork, const bool & doNotDoAdjust, const bool & updateImpImmed)
{

	map<unsigned int, double> nodeDataCtsValues; //node ID, value at indiv i
	map<unsigned int, unsigned int> nodeDataDisValues; //node ID, value at indiv i

	map<unsigned int, double> groupNodesBestDistance; //group ID, shortest distance to NN
	map<unsigned int, int> groupNodesBestSplitter; //group ID, best random no used to choose between pts with the same distance
	map<unsigned int, double> nodeDataBestCtsValue; //(missing data) node ID, best found cts value for this node
	map<unsigned int, unsigned int> nodeDataBestDisValue; //(missing data) node ID, best found dis value for this node
	map<unsigned int, double> partialVarStDevs;
	unsigned int groupIDNodeId;

	map<unsigned int, set<unsigned int> >::const_iterator gnwcd;
	map<unsigned int, set<unsigned int> >::const_iterator gnwmd;
	map<unsigned int, double>::iterator bestDis;
	map<unsigned int, int>::iterator bestDisSplit; 
	map<unsigned int, double>::iterator nodeBestCtsVal; 
	map<unsigned int, unsigned int>::iterator nodeBestDisVal;
	map<unsigned int, double>::const_iterator psd;
	set<unsigned int> usedNums; // if distance is the same draw a random no 1 to noIndivs, these are the used nums. Smallest no wins
	double distance, diff, stDevPartial;
	bool doPartialCalc, doLoopDisGroup, vetoDoPartialCalc2, vetoDoPartialCalc3;
	map<unsigned int, Node *> theParents;
	list<unsigned int> bootNodeIDsAdjVars;
	list<unsigned int> adjMethods;
	double adjValueCurrentIndiv;
	map<unsigned int, double> adjValuesCurrentIndiv; //node ID for complete data node surrounding imputed nodes adjusted value
	map<unsigned int, double>::const_iterator avci;
	bool currentAdjCalced = false;
	unsigned int thisParentBootID, otherParentOrigID;
	int randNo;
	bool validNN;
	bool setAsNN;
	unsigned int bNodeID;
	Node * aNode;
	Node * aNodeBoot;
	Node * bNode;
	//bool noProbAdj = false;//true;

	unsigned int noLevels;
	double diff0;


	//loop thro' data using a different iterator to where the NN is up to, to find nearest neighbours
	setInitialDataIterators3();

	unsigned int amountData = allNodeData->getAmountOfData();

//	loop thro' indivs iter 3
	for(unsigned int indiv = 1; indiv <= amountData; ++indiv)
	{

		validNN = true;

		//loop thro' nodes(group of nodes) with missing data, and corresponding connected complete data
		gnwmd = groupNodesWithMissingData.begin();
		for(gnwcd = groupNodesWithCompleteData.begin(); gnwcd != groupNodesWithCompleteData.end() && validNN; ++gnwcd, ++gnwmd)
		{
			
			distance = 0;		

			//check data is complete for the missing data nodes, iter3
			for(set<unsigned int>::const_iterator nd = gnwmd->second.begin(); nd != gnwmd->second.end() && validNN; ++nd)
			{		
				aNode = getNetworkNode(*nd);
				if(aNode->nodeDataIsMissing3()) validNN = false;				
			};

			if(validNN && gnwcd->second.size() != 0)
			{

				//calculate distance to current indiv, if all data is complete that we need
				for(set<unsigned int>::const_iterator nd = gnwcd->second.begin(); nd != gnwcd->second.end() && validNN; ++nd)
				{

					aNode = getNetworkNode(*nd);
					if(!aNode->nodeDataIsMissing3())
					{

						if(aNode->getIsDiscreteNode())
						{						
							if(aNode->getNodeDataDisValue2() == aNode->getNodeDataDisValue3()) diff = 0;
							else diff = 1;					
						}
						else
						{													
							//check if exists other parents other than missing nodes being imputed together
							doPartialCalc = false;
							vetoDoPartialCalc2 = false;
							vetoDoPartialCalc3 = false;
							doLoopDisGroup = false;

							if(freeClearMemory)
							{
								list<unsigned int>().swap(bootNodeIDsAdjVars);
								list<unsigned int>().swap(adjMethods);
							}
							else
							{
								bootNodeIDsAdjVars.clear();
								adjMethods.clear();
							}
							theParents = bootstrapNetwork->getNetworkNode(bootstrapNetwork->convertID(*nd))->getParents();

							//check that the node is a parent of a node in the missing node group firstly, then check for other parents - do we need to adjust with partial calc?					
							for(set<unsigned int>::const_iterator ndmiss = gnwmd->second.begin(); ndmiss != gnwmd->second.end() && validNN; ++ndmiss)
							{
								thisParentBootID = bootstrapNetwork->convertID(*ndmiss);
								aNodeBoot = bootstrapNetwork->getNetworkNode(thisParentBootID);								
								if(theParents.find(thisParentBootID) != theParents.end()) //missing node is a parent and not discrete
								{

									if(aNodeBoot->getIsDiscreteNode()) //if a missing node is a parent and discrete then cannot adjust as we don't know group
									{
										doPartialCalc = false;
										break;
										doLoopDisGroup = true;
									};

									for(map<unsigned int, Node *>::const_iterator pa = theParents.begin(); pa != theParents.end(); ++pa)
									{
							
										//convert to orig nodes
										otherParentOrigID = convertID(pa->first);										
										
										if(thisParentBootID != pa->first)
										{
											if(pa->second->getIsDiscreteNode())
											{		
												if(getNetworkNode(otherParentOrigID)->nodeDataIsMissing2()) vetoDoPartialCalc2 = true;
												else if(getNetworkNode(otherParentOrigID)->nodeDataIsMissing3()) vetoDoPartialCalc3 = true; //if parent is discrete and missing of variable to be adjusted then cannot calc the adj - so skip this indiv if adj req'd
											}
											else if(gnwmd->second.find(otherParentOrigID) == gnwmd->second.end())  //check is not in this missing group node data	
											{
												if(!pa->second->parentExists(thisParentBootID)) //check parents are not connected to one another
												{												
													if(!aNodeBoot->parentExists(pa->first))
													{	
														if(!getNetworkNode(otherParentOrigID)->nodeDataIsMissing2()) //only add adjust variable if non missing for impute indiv
														{
															bootNodeIDsAdjVars.push_back(pa->first);
															doPartialCalc = true;
														};
													};
												};
											};
										};
									};
								};			
							};
									
							
							if(doPartialCalc && !vetoDoPartialCalc2)
							{
								if(vetoDoPartialCalc3)
								{
									validNN = false;
								}
								else
								{
									//check other parent nodes have data, determine adjust methods
									for(list<unsigned int>::const_iterator bv = bootNodeIDsAdjVars.begin(); bv != bootNodeIDsAdjVars.end(); ++bv)
									{
										bNodeID = convertID(*bv);  //convert to orig nodes		
										bNode = getNetworkNode(bNodeID); //get node in original data
									
										if(!bNode->nodeDataIsMissing2() && !doNotDoAdjust)
										{
											adjMethods.push_back(0); //adjust as A - beta_B * B
											if(bNode->nodeDataIsMissing3()) validNN = false; // for this potential NN cannot adjust as so. Will be some NNs, as there must be some complete data Indivs
											
											if(bNode->nodeDataIsMissing2() || bNode->nodeDataIsMissing3()) validNN = false;
										}
										else
										{
											adjMethods.push_back(1); //adjust as A, i.e. do not adjust
											//adjMethods.push_back(2); //adjust as A - beta_B * mean(B)
											//adjMethods.push_back(3); //adjust as A - beta_B * random(B)
										
											//if(noProbAdj) validNN = false; //skip this individual
										};

										//if(bNode->nodeDataIsMissing2() || bNode->nodeDataIsMissing3()) validNN = false;								
									};
								};							

								if(validNN)
								{
																		
									aNodeBoot = bootstrapNetwork->getNetworkNode(bootstrapNetwork->convertID(*nd)); //complete data node
												
									groupIDNodeId = aNodeBoot->getNodeAndParentsGroupNo2(this); //use orig network for data, boot net for structure
									
									if(!doLoopDisGroup)
									{
										psd = partialVarStDevs.find(*nd);
										avci = adjValuesCurrentIndiv.find(*nd);
										if(avci != adjValuesCurrentIndiv.end())
										{
											stDevPartial = psd->second;
											adjValueCurrentIndiv = avci->second;
										}
										else
										{
										
											stDevPartial = aNode->getNodeDataStDevAdj(gnwmd->second, bootNodeIDsAdjVars, adjMethods, this, bootstrapNetwork);
											partialVarStDevs[*nd] = stDevPartial;
		
											adjValueCurrentIndiv = aNodeBoot->getNodeDataCtsValueAdj(true, gnwmd->second, bootNodeIDsAdjVars, adjMethods, this, bootstrapNetwork, doLoopDisGroup, 0);
											adjValuesCurrentIndiv[*nd] = adjValueCurrentIndiv;										
										};
				
										diff = (adjValueCurrentIndiv - aNodeBoot->getNodeDataCtsValueAdj(false, gnwmd->second, bootNodeIDsAdjVars, adjMethods, this, bootstrapNetwork, doLoopDisGroup, 0))/stDevPartial; //adjusted for partials
									}
									else
									{
										//try each group for imputing node
										psd = partialVarStDevs.find(*nd);
										avci = adjValuesCurrentIndiv.find(*nd);
										if(psd != partialVarStDevs.end())
										{
											stDevPartial = psd->second;											
										}
										else
										{										
											stDevPartial = aNode->getNodeDataStDevAdj(gnwmd->second, bootNodeIDsAdjVars, adjMethods, this, bootstrapNetwork);
											partialVarStDevs[*nd] = stDevPartial;																					
										};
										
										noLevels = 2;//aNodeBoot->getNoLevels(); //no groups not levs...
										
										for(unsigned int lev = 1; lev <= noLevels; ++lev)
										{
										
											adjValueCurrentIndiv = aNodeBoot->getNodeDataCtsValueAdj(true, gnwmd->second, bootNodeIDsAdjVars, adjMethods, this, bootstrapNetwork, doLoopDisGroup, lev);

											diff0 = (adjValueCurrentIndiv - aNodeBoot->getNodeDataCtsValueAdj(false, gnwmd->second, bootNodeIDsAdjVars, adjMethods, this, bootstrapNetwork, doLoopDisGroup, 0))/stDevPartial; //adjusted for partials

											if(lev == 1 || fabs(diff0) < fabs(diff)) diff = diff0;
											
										};
																			
										adjValuesCurrentIndiv[*nd] = adjValueCurrentIndiv;
									};
							
								};
								
							}
							else
							{
								diff = (aNode->getNodeDataCtsValue2() - aNode->getNodeDataCtsValue3())/aNode->getStDev(); 						
							};
						};//end of discrete or cts complete data node

						distance += diff*diff;
					}
					else
					{
						validNN = false;  //complete data node is missing for this potential NN (iter 3)					
					};
				};//end of complete data nodes connected to missing nodes to be imputed
			};//end of if there are any complete data nodes

			//update nearest neighbour if better than best previous
			if(validNN)
			{				
				setAsNN = false;

				bestDis = groupNodesBestDistance.find(gnwmd->first);

				//set NN if not set previously
				if(bestDis == groupNodesBestDistance.end())
				{
					setAsNN = true;
				}
				else if(distance == 0 || distance == bestDis->second)
				{
					bestDisSplit = groupNodesBestSplitter.find(gnwmd->first);
					if(bestDisSplit == groupNodesBestSplitter.end())
					{
						groupNodesBestSplitter[gnwmd->first] = rand();
						bestDisSplit = groupNodesBestSplitter.find(gnwmd->first);
					};
						
					randNo = rand();
					if(randNo < bestDisSplit->second)
					{
						groupNodesBestSplitter[gnwmd->first] = randNo;
						setAsNN = true;
					};

				}
				else if(distance < bestDis->second)
				{
					setAsNN = true;
				};
					
				if(setAsNN)
				{					
					groupNodesBestDistance[gnwmd->first] = distance;
					for(set<unsigned int>::const_iterator nd = gnwmd->second.begin(); nd != gnwmd->second.end() && validNN; ++nd)
					{
						aNode = getNetworkNode(*nd);
						
						if(aNode->getIsDiscreteNode())
						{
							nodeDataBestDisValue[*nd] = aNode->getNodeDataDisValue3(); 
						}
						else
						{						
							nodeDataBestCtsValue[*nd] = aNode->getNodeDataCtsValue3(); 
						};
					};
				};

			};//end of checking if this indiv may be NN
		
		}; // end of missing data node (or nodes) loop, nodes that are imputed together

		advanceDataIterators3();

	};//end of indivs loop

	//set missing data with imputed data
	for(map<unsigned int, double>::const_iterator ndCts = nodeDataBestCtsValue.begin(); ndCts != nodeDataBestCtsValue.end(); ++ndCts)
	{		
		aNode = getNetworkNode(ndCts->first);	
		aNode->setImputedDataCts(ndCts->second, updateImpImmed);
	};

	for(map<unsigned int, unsigned int>::const_iterator ndDis = nodeDataBestDisValue.begin(); ndDis != nodeDataBestDisValue.end(); ++ndDis)
	{		
		aNode = getNetworkNode(ndDis->first);
		aNode->setImputedDataDis(ndDis->second, updateImpImmed);
	};
	
	//calculate how many points failed to impute
	unsigned int noToImpute = 0;
	for(map<unsigned int, set<unsigned int> >::const_iterator i = groupNodesWithMissingData.begin(); i != groupNodesWithMissingData.end(); ++i)
	{
		noToImpute += i->second.size();

		/*for(set<unsigned int>::const_iterator j = i->second.begin(); j != i->second.end(); ++j)
		{
			aNode = getNetworkNode(*j);
		};*/
	};

	return noToImpute - (nodeDataBestCtsValue.size() + nodeDataBestDisValue.size());
};

//! Replaces missing data with the mean for each variable (terrible method).  
unsigned int Network::imputeDataMean(const map<unsigned int, set<unsigned int> > & groupNodesWithMissingData, Network * bootstrapNetwork, const bool & updateImpImmed)
{
	map<unsigned int, double> nodeDataBestCtsValue; //(missing data) node ID, best found cts value for this node
	map<unsigned int, unsigned int> nodeDataBestDisValue; //(missing data) node ID, best found dis value for this node

	map<unsigned int, set<unsigned int> >::const_iterator gnwmd;
	
	Node * aNode;
	
	//loop thro' nodes(group of nodes) with missing data, and corresponding connected complete data		
	for(gnwmd = groupNodesWithMissingData.begin(); gnwmd != groupNodesWithMissingData.end(); ++gnwmd)
	{	
					
		for(set<unsigned int>::const_iterator nd = gnwmd->second.begin(); nd != gnwmd->second.end(); ++nd)
		{
			aNode = getNetworkNode(*nd);

			if(aNode->getIsDiscreteNode())
			{
				//nodeDataBestDisValue[*nd] = aNode->getNodeMode(); //maybe later
			}
			else
			{
				nodeDataBestCtsValue[*nd] = aNode->getMean(); // mean will be set after the first use, so mean will remain the same for all missing data for a given variable
			};
		};
				
	}; // end of missing data node (or nodes) loop, nodes that are imputed together

	//set missing data with imputed data
	for(map<unsigned int, double>::const_iterator ndCts = nodeDataBestCtsValue.begin(); ndCts != nodeDataBestCtsValue.end(); ++ndCts)
	{
		aNode = getNetworkNode(ndCts->first);
		aNode->setImputedDataCts(ndCts->second, updateImpImmed);
	};

	for(map<unsigned int, unsigned int>::const_iterator ndDis = nodeDataBestDisValue.begin(); ndDis != nodeDataBestDisValue.end(); ++ndDis)
	{
		aNode = getNetworkNode(ndDis->first);
		aNode->setImputedDataDis(ndDis->second, updateImpImmed);
	};

	//calculate how many points failed to impute
	unsigned int noToImpute = 0;
	for(map<unsigned int, set<unsigned int> >::const_iterator i = groupNodesWithMissingData.begin(); i != groupNodesWithMissingData.end(); ++i)
	{
		noToImpute += i->second.size();

		/*for(set<unsigned int>::const_iterator j = i->second.begin(); j != i->second.end(); ++j)
		{
			aNode = getNetworkNode(*j);
		};*/
	};

	return noToImpute - (nodeDataBestCtsValue.size() + nodeDataBestDisValue.size());
};



//! Returns a string for the network.
string Network::getNetworkString(const unsigned int & charLimit)
{
	string networkString = "";

	//set all nodes to unvisited
	clearVisitedNodes();

	unsigned int noVisitedThisLoop = 0;
	unsigned int totalVisited = 0;
	unsigned int totalNodes = allNetworkNodes.size();
	bool loopDetected = false;
	set<string> nodeCondNames;
	list<Node *> visitedNodesThisLoop;

	//try to visit every node if the parents are given (if any), if no loops then can visit every loop
	//output even if a string is detected, but try to do in order first
	do{
		noVisitedThisLoop = 0;
		nodeCondNames.clear();
		visitedNodesThisLoop.clear();

		for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
		{
			if(!nd->second->getVisited() && (nd->second->allParentsVisited() || loopDetected))
			{
				if(!nd->second->getIsFactorChildNode()) nodeCondNames.insert(nd->second->getCondName()); //make alpha

				//nd->second->setVisited(true);
				visitedNodesThisLoop.push_back(nd->second);
				++noVisitedThisLoop;
			};
		};

		if(noVisitedThisLoop == 0) loopDetected = true;

		for(set<string>::const_iterator ncn = nodeCondNames.begin(); ncn != nodeCondNames.end(); ++ncn)
		{
			networkString += *ncn;			
		};
		
		for(list<Node *>::iterator vntl = visitedNodesThisLoop.begin(); vntl != visitedNodesThisLoop.end(); ++vntl)
		{
			(*vntl)->setVisited(true);
		};

		totalVisited += noVisitedThisLoop;
	}while(totalVisited < totalNodes);

	if(charLimit != 0)
	{
		unsigned int length = networkString.length();
		if(length > charLimit)
		{
			networkString = networkString.substr(0, charLimit) + "...";
		};
	};

	//should has visited all nodes if no loops
	return networkString;
};

//! Returns a CPDAG string for the network
string Network::getPDAGNetworkString()
{
	if(hasLoop())
	{
		displayLoop();
		out("\n");
		out(getNetworkString(0));
		out("\n");
		
		exitErr("Unable to get string for PDAG as the network has a loop in it!");
	};
	 
	set<string> vStructures; // collider node + node1 pa + node2 pa
	set<string> vStrConnections; // collider node + parent node, used to check if already exists when setting up nodeConnections
	map<string, set<string> > nodeConnectUndir;
	map<string, set<string> >::iterator ndc;
	map<string, set<string> > nodeConnectDir;
	set<string> isolatedNodes;
	map<string, string>::iterator itn;
	set<string>::const_iterator v;
	set<string>::const_iterator v2;
	set<string> someNodes;
	
	map<unsigned int, Node *> someParents;
	string nodeName1, nodeName2, nodeNameTemp;
	set<string> noParentNodes;
	map<unsigned int, Node *>::const_iterator pa1, pa2;
	string vStructString, connectedNodes;

	//set up structure of the network firstly and then construct the string
	//set up v-structures
	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		someParents = nd->second->getParents();
		
		if(someParents.size() > 1)
		{
			//loop thro' pairs of parents
			pa1 = someParents.begin();
			pa2 = someParents.begin(); ++pa2;

			while(pa1 != someParents.end() && pa2 != someParents.end())
			{
				if(!edgeExists(pa1->first, pa2->first) && !edgeExists(pa2->first, pa1->first))
				{
					nodeName1 = pa1->second->getDisplayName();
					nodeName2 = pa2->second->getDisplayName();

					if(nodeName1 > nodeName2) //order alphabetically
					{
						nodeNameTemp = nodeName1;
						nodeName1 = nodeName2;
						nodeName2 = nodeNameTemp;
					};

					vStructString = nd->second->getDisplayName() + "<" + nodeName1 + ":" + nodeName2;
					vStructures.insert(vStructString);

					if(nodeName1 < nd->second->getDisplayName()) connectedNodes = nodeName1 + " " + nd->second->getDisplayName();
					else connectedNodes = nd->second->getDisplayName() + " " + nodeName1;
					vStrConnections.insert(connectedNodes); 

					if(nodeName2 < nd->second->getDisplayName()) connectedNodes = nodeName2 + " " + nd->second->getDisplayName();
					else connectedNodes = nd->second->getDisplayName() + " " + nodeName2;
					vStrConnections.insert(connectedNodes);
				};

				++pa2;
				if(pa2 == someParents.end())
				{
					++pa1; 
					pa2 = someParents.begin();
					//start pa2 one after pa1
					if(pa1 != someParents.end())
					{
						while(pa2->first != pa1->first) ++pa2;
						++pa2;
					};
				};
			};
		}
		else if(someParents.size() == 0) noParentNodes.insert(nd->second->getDisplayName());

	};

	bool connectionFound;
	map<unsigned int, Node *> parentsOfParent;

	//set up non-directed connections and directed connections that are not v-strs
	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		someParents = nd->second->getParents();
		
		for(map<unsigned int, Node *>::const_iterator pa = someParents.begin(); pa != someParents.end(); ++pa)
		{
			nodeName1 = nd->second->getDisplayName();
			nodeName2 = pa->second->getDisplayName();

			if(nodeName1 > nodeName2) //order alphabetically
			{
				nodeNameTemp = nodeName1;
				nodeName1 = nodeName2;
				nodeName2 = nodeNameTemp;
			};

			//check if already in v-strs
			connectionFound = false;
			connectedNodes = nodeName1 + " " + nodeName2;

			v = vStrConnections.find(connectedNodes);
			if(v != vStrConnections.end())
			{
				connectionFound = true;
			}
			else
			{
				connectedNodes = nodeName1 + " " + nodeName2;
				v = vStrConnections.find(connectedNodes);
				if(v != vStrConnections.end())
				{
					connectionFound = true;
				};
			};

			if(!connectionFound)
			{
				ndc = nodeConnectUndir.find(nodeName1);
				if(ndc != nodeConnectUndir.end())
				{
					ndc->second.insert(nodeName2);
				}
				else
				{
					someNodes.clear();
					someNodes.insert(nodeName2);
					nodeConnectUndir[nodeName1] = someNodes;	
				};								
			};

		};
		
	};

	bool nodeFound;
	unsigned int pos1;

	//check if some nodes have no connections
	for(set<string>::const_iterator npn = noParentNodes.begin(); npn != noParentNodes.end(); ++npn)
	{
		nodeFound = false;
		for(map<string, set<string> >::const_iterator nc = nodeConnectUndir.begin(); nc != nodeConnectUndir.end() && !nodeFound; ++nc)
		{
			if(*npn == nc->first) nodeFound = true;

			for(set<string>::const_iterator otn = nc->second.begin(); otn != nc->second.end() && !nodeFound; ++otn)
			{
				if(*npn == *otn) nodeFound = true;
			};	
		};

		for(set<string>::const_iterator vs = vStrConnections.begin(); vs != vStrConnections.end() && !nodeFound; ++vs)
		{
			pos1 = (unsigned int)(*vs).find_first_of(' ');

			nodeName1 = (*vs).substr(0, pos1);
			nodeName2 = (*vs).substr(pos1 + 1);

			if(*npn == nodeName1 || *npn == nodeName2) nodeFound = true;
		};

		if(!nodeFound) isolatedNodes.insert(*npn);
	};

	string CPDAGString = "";

	//isolated nodes
	for(set<string>::const_iterator isn = isolatedNodes.begin(); isn != isolatedNodes.end(); ++isn)
	{
		CPDAGString += "(" + *isn + ")";
	};	

	//undirected connections
	for(map<string, set<string> >::const_iterator nc = nodeConnectUndir.begin(); nc != nodeConnectUndir.end(); ++nc)
	{
		CPDAGString += "(" + nc->first + "|";
		for(set<string>::const_iterator otn = nc->second.begin(); otn != nc->second.end(); )
		{
			CPDAGString += *otn;
			++otn;
			if(otn != nc->second.end()) CPDAGString += ":";
		};
		CPDAGString += ")";
	};

	//directed connections
	for(map<string, set<string> >::const_iterator ncu = nodeConnectDir.begin(); ncu != nodeConnectDir.end(); ++ncu)
	{
		CPDAGString += "(" + ncu->first + "^";
		for(set<string>::const_iterator otnu = ncu->second.begin(); otnu != ncu->second.end(); )
		{
			CPDAGString += *otnu;
			++otnu;
			if(otnu != ncu->second.end()) CPDAGString += ":";
		};
		CPDAGString += ")";
	};

	//directed v-str connections
	for(set<string>::const_iterator vs = vStructures.begin(); vs != vStructures.end(); ++vs)
	{
		CPDAGString += "(" + *vs + ")";
	};

	return CPDAGString;
};

//! Calculates all network scores given by the nodes in the initial network.
unsigned int Network::calculateAllScores(const string & allScoresFilename, double & bestNetworkScore, string & bestNetwork)
{
	calculateAllScoresSetup();

	cacheAllNodes();

	ofstream allScoresFile(allScoresFilename.c_str());

	bool updated = true;	
	double score;
	bool firstScore = true;
	unsigned int noNetworksEval = 0;

	while(updated)
	{
		
		//evaluate network
		if(!hasLoop())
		{
			noNetworksEval++;

			score = getScoreFromCache();

			if((firstScore || score > bestNetworkScore) && score*0 == 0)
			{
				bestNetworkScore = score;
				bestNetwork = getNetworkString();
				firstScore = false;
			};

			allScoresFile << noNetworksEval << " " << score << " " << getNetworkString() << "\n";	

		};

		updateNextNetwork(updated);
	};
	
	allScoresFile.close();

	return noNetworksEval;
};

//! Calculates all network scores given by the nodes in the initial network.
void Network::calculateAllScoresSetup()
{
	
	//remove all edges firstly
	removeAllEdges();

	//add white edges to network
	for(set<pair<unsigned int, unsigned int> >::const_iterator we = whiteEdges.begin(); we != whiteEdges.end(); ++we)
	{
		if(!edgeAllowed(we->first, we->second))
		{
			outErr("Edge "); outErr(getNetworkNode(we->first)->getDisplayName()); outErr(" --> "); outErr(getNetworkNode(we->second)->getDisplayName()); outErr(" is not allowed, but the white edge list states otherwise!\n");
			exitErr("");
		}
		else if(isNodeCts(we->first) && isNodeDiscrete(we->second))
		{
			outErr("Edge "); outErr(getNetworkNode(we->first)->getDisplayName()); out(" --> "); outErr(this->getNetworkNode(we->second)->getDisplayName()); outErr(" is not allowed as it is discrete --> cts, but is stated in the list of white edges!\n");
			exitErr("");
		};

		addEdge(we->first, we->second);
	};

	
	for(map<unsigned int, Node *>::const_iterator nd1 = allNetworkNodes.begin(); nd1 != allNetworkNodes.end(); ++nd1)
	{
		for(map<unsigned int, Node *>::const_iterator nd2 = allNetworkNodes.begin(); nd2 != allNetworkNodes.end(); ++nd2)
		{
			if(nd1->first < nd2->first 
				&& isWhiteEdge(nd1->first, nd2->first) && isWhiteEdge(nd2->first, nd1->first) //edges that must be present but in any direction 
				&& !nd1->second->getIsFactorChildNode() && !nd2->second->getIsFactorChildNode()
				&& (edgeAllowed(nd1->first, nd2->first) || edgeAllowed(nd2->first, nd1->first)))				
			{
				nodePairs.push_back(make_pair(nd1->first, nd2->first));
				nodePairConnections.push_back(2);
			}
			else if(nd1->first < nd2->first 
				&& !isWhiteEdge(nd1->first, nd2->first) && !isWhiteEdge(nd2->first, nd1->first)
				&& !nd1->second->getIsFactorChildNode() && !nd2->second->getIsFactorChildNode()
				&& (edgeAllowed(nd1->first, nd2->first) || edgeAllowed(nd2->first, nd1->first)))			
			{
				nodePairs.push_back(make_pair(nd1->first, nd2->first));
				nodePairConnections.push_back(1);
			};
		};
	};

};

//! Calculates all network scores given by the nodes in the initial network.
void Network::updateNextNetwork(bool & updated)
{

	//update network	
	updated = false;
	list<pair<unsigned int, unsigned int> >::const_iterator ndPrs = nodePairs.begin();
	for(list<unsigned int>::iterator ndPairCh = nodePairConnections.begin(); ndPairCh != nodePairConnections.end(); ++ndPairCh, ++ndPrs)
	{
		if(*ndPairCh < 3)
		{
			(*ndPairCh)++;

			//check is ok, otherwise update again
			if(*ndPairCh == 2 && !edgeAllowed(ndPrs->first, ndPrs->second))				
			{
				(*ndPairCh)++; //i.e. set to 3 because 2 not good
			};
				
			if(*ndPairCh == 3 && !edgeAllowed(ndPrs->second, ndPrs->first))				
			{
				if(isWhiteEdge(ndPrs->first, ndPrs->second) && isWhiteEdge(ndPrs->second, ndPrs->first)) *ndPairCh = 2;
				else *ndPairCh = 1; //i.e. set to 1 because 3 not good
			}
			else
			{
				//set prev nodes to 1
				updated = true;
				break;
			};
		}
		else
		{
			if(isWhiteEdge(ndPrs->first, ndPrs->second) && isWhiteEdge(ndPrs->second, ndPrs->first)) *ndPairCh = 2;
			else *ndPairCh = 1;
		};

	};

	if(!updated) return;

	//update network to match what is given in nodePairConnections
	//check for black and white nodes and discrete, parent allowed nodes
	set<unsigned int> updatedNodes;
	
	list<pair<unsigned int, unsigned int> >::const_iterator ndPairs = nodePairs.begin();
	for(list<unsigned int>::const_iterator ndPairUp = nodePairConnections.begin(); ndPairUp != nodePairConnections.end(); ++ndPairUp, ++ndPairs)
	{
		
		if(*ndPairUp == 1) //remove edge if exists
		{
			if(edgeExists(ndPairs->first, ndPairs->second)) {removeEdge(ndPairs->first, ndPairs->second); updatedNodes.insert(ndPairs->second);}
			else if(edgeExists(ndPairs->second, ndPairs->first)) {removeEdge(ndPairs->second, ndPairs->first); updatedNodes.insert(ndPairs->first);};
		}
		else if(*ndPairUp == 2) //add edge 1 --> 2
		{	
			if(edgeExists(ndPairs->second, ndPairs->first))
			{
				removeEdge(ndPairs->second, ndPairs->first);
				addEdge(ndPairs->first, ndPairs->second);	
				updatedNodes.insert(ndPairs->first);
				updatedNodes.insert(ndPairs->second);
			}
			else if(!edgeExists(ndPairs->first, ndPairs->second)) {addEdge(ndPairs->first, ndPairs->second); updatedNodes.insert(ndPairs->second);};
		}
		else if(*ndPairUp == 3) //add edge 2 --> 1
		{	
			if(edgeExists(ndPairs->first, ndPairs->second))
			{
				removeEdge(ndPairs->first, ndPairs->second);
				addEdge(ndPairs->second, ndPairs->first);
				updatedNodes.insert(ndPairs->first);
				updatedNodes.insert(ndPairs->second);
			}
			else if(!edgeExists(ndPairs->second, ndPairs->first)) {addEdge(ndPairs->second, ndPairs->first); updatedNodes.insert(ndPairs->first);};
		};

	}; //end of updating network

	for(set<unsigned int>::const_iterator un = updatedNodes.begin(); un != updatedNodes.end(); ++un)
	{
		setScoreBitInCache(*un);
	};

};

//! Check the current network is valid
bool Network::checkNetworkIsValid()
{
	map<unsigned int, Node *> parents;
	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		parents = nd->second->getParents();
		for(map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end(); ++pa)
		{
			if(!edgeAllowed(pa->first, nd->first) && !pa->second->getIsFactorChildNode())
			{
				return false;
			};
		};
		
	};

	return true;
};


//! Clears cache for nodes, cts node linear reg for initial prior
void Network::clearCache()
{
	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		nd->second->clearCache();		
	};
};

//! Calculates score differences in changing the edges.
void Network::evaluateAllEdgeChanges()
{
	//clear list of previos evaluations
	
	if(freeClearMemory) multimap<double, pair<unsigned int, unsigned int> >().swap(evaluations); else evaluations.clear();

	double score;
	
	//loop thro' pairs of nodes 
	for(map<unsigned int, Node *>::iterator nd1 = allNetworkNodes.begin(); nd1 != allNetworkNodes.end(); ++nd1)
	{
		for(map<unsigned int, Node *>::iterator nd2 = allNetworkNodes.begin(); nd2 != allNetworkNodes.end(); ++nd2)
		{
			if(nd1->first != nd2->first)
			{ 
	
			if(nd2->second->parentExists(nd1->first)) // nd1 --> nd2 exists
			{
				if((!isWhiteEdge(nd1->first, nd2->first) || (isWhiteEdge(nd1->first, nd2->first) && isWhiteEdge(nd2->first, nd1->first))) //can switch white edge if white in both directions					
					&& edgeAllowed(nd2->first, nd1->first)
					)
				{
					//try REVERSING edge, nd1 <-- nd2	
					removeEdge(nd1->first, nd2->second);
					addEdge(nd2->second, nd1->second);

					if(!hasLoop())
					{						
						score = calcScore();
						evaluations.insert(make_pair(score, make_pair(nd2->first, nd1->first)));						
					};

					//change back: now put back the edge, nd1 --> nd2 
					removeEdge(nd2->first, nd1->second);
					addEdge(nd1->second, nd2->second);
				};
			}
			else if(nd1->second->parentExists(nd2->first))  // nd1 <-- nd2 exists
			{
				if(!isWhiteEdge(nd2->first, nd1->first))
				{
					//try DELETING edge, nd2 --> nd1
					removeEdge(nd2->first, nd1->second);
					score = calcScore();
					evaluations.insert(make_pair(score, make_pair(nd2->first, nd1->first)));
					
					//change back: now put back the edge, nd2 --> nd1 	
					addEdge(nd2->second, nd1->second);		
				};
			}
			else
			{			
				//try ADDING edge nd1 <-- nd2
				if(edgeAllowed(nd2->first, nd1->first))
				{							
					addEdge(nd2->second, nd1->second);
	
					if(!hasLoop())
					{					
						score = calcScore();
						evaluations.insert(make_pair(score, make_pair(nd2->first, nd1->first)));				
					};

					//change back: remove edge				
					removeEdge(nd2->first, nd1->second);
				};
		
			};			
		
			};//end of nodes not the same
		};		
				
	};

};

//! Calculates and caches score bit for every node and returns network score.
double NetworkBNLearn::cacheAllNodes()
{
	if(!checkNetworkIsValid())
	{
		string mess = "Attempt to evaluate an invalid network!\n" + getNetworkString() + "\n";	
		exitErr(mess);
	};

	double score = 0;
	double scoreBit;

	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{	
		scoreBit = nd->second->calcUpdatedScoreBit(this);
		nodeScoreCache[nd->first] = scoreBit;
		score += scoreBit;
	};

	return score;
};

//! Calculates and caches score bit for every node and returns network score.
double NetworkBNLearn::getScoreFromCache()
{
	double score = 0;

	for(map<unsigned int, double>::const_iterator nsc = nodeScoreCache.begin(); nsc != nodeScoreCache.end(); ++nsc)
	{	
		score += nsc->second; //add score bit
	};

	return score;
};

//! Sets score bit for one node in the cache.
void NetworkBNLearn::setScoreBitInCache(const unsigned int & nodeID, const double & value)
{
	nodeScoreCache[nodeID] = value;
};

//! Calculates and sets score bit
void NetworkBNLearn::setScoreBitInCache(const unsigned int & nodeID)
{
	nodeScoreCache[nodeID] = getNetworkNode(nodeID)->calcUpdatedScoreBit(this);
};

double NetworkBNLearn::getNodeScoreFromCache(const unsigned int & nodeID)
{
	map<unsigned int, double>::const_iterator nsc;
	nsc = nodeScoreCache.find(nodeID);

	if(nsc != nodeScoreCache.end()) return nsc->second;

	return 0;
};

//! Calculates the updated network score when only 2 nodes have changed.
double NetworkBNLearn::calcScore2DifferentNodes(const double & origScore, const unsigned int & nodeID1, const unsigned int & nodeID2)
{
	double score = origScore;

	for(map<unsigned int, double>::const_iterator nsc = nodeScoreCache.begin(); nsc != nodeScoreCache.end(); ++nsc)
	{
		if(nsc->first == nodeID1 || nsc->first == nodeID2)
		{
			score -= nsc->second; //take away original score
			score += getNetworkNode(nsc->first)->calcUpdatedScoreBit(this); //add new score
		};

	};

	return score;
};

//! Calculates score differences in changing the edges, in a more efficient way for bnlearn network.
void NetworkBNLearn::evaluateAllEdgeChangesDiffs()
{
	//cache all network node evals
	double origScore = getScoreFromCache();

	double scoreDiff;
	
	 //clear and free memory
	if(freeClearMemory) multimap<double, pair<unsigned int, unsigned int> >().swap(evaluations); else evaluations.clear();

	
	//loop thro' pairs of nodes 
	for(map<unsigned int, Node *>::iterator nd1 = allNetworkNodes.begin(); nd1 != allNetworkNodes.end(); ++nd1)
	{
		for(map<unsigned int, Node *>::iterator nd2 = allNetworkNodes.begin(); nd2 != allNetworkNodes.end(); ++nd2)
		{
			if(nd1->first != nd2->first)
			{ 
	
				if(nd2->second->parentExists(nd1->first)) // nd1 --> nd2 exists
				{
					if((!isWhiteEdge(nd1->first, nd2->first) || (isWhiteEdge(nd1->first, nd2->first) && isWhiteEdge(nd2->first, nd1->first))) //can switch white edge if white in both directions					
						&& edgeAllowed(nd2->first, nd1->first)
						)
					{
						//try REVERSING edge, nd1 <-- nd2	
						removeEdge(nd1->first, nd2->second);
						addEdge(nd2->second, nd1->second);

						//if(!hasLoop())
						//{						
							scoreDiff = calcScore2DifferentNodes(origScore, nd1->first, nd2->first) - origScore;
							if(scoreDiff > 0) evaluations.insert(make_pair(scoreDiff, make_pair(nd2->first, nd1->first)));						
						//};

						//change back: now put back the edge, nd1 --> nd2 
						removeEdge(nd2->first, nd1->second);
						addEdge(nd1->second, nd2->second);
					};
				}
				else if(nd1->second->parentExists(nd2->first))  // nd1 <-- nd2 exists
				{
					if(!isWhiteEdge(nd2->first, nd1->first))
					{
						//try DELETING edge, nd2 --> nd1
						removeEdge(nd2->first, nd1->second);
						scoreDiff = calcScore2DifferentNodes(origScore, nd1->first, 0) - origScore;					
						if(scoreDiff > 0) evaluations.insert(make_pair(scoreDiff, make_pair(nd2->first, nd1->first)));
					
						//change back: now put back the edge, nd2 --> nd1 	
						addEdge(nd2->second, nd1->second);		
					};
				}
				else
				{			
					//try ADDING edge nd1 <-- nd2
					if(edgeAllowed(nd2->first, nd1->first))
					{							
						addEdge(nd2->second, nd1->second);
	
						//if(!hasLoop())
						//{					
							scoreDiff = calcScore2DifferentNodes(origScore, nd1->first, 0) - origScore;					
							if(scoreDiff > 0) evaluations.insert(make_pair(scoreDiff, make_pair(nd2->first, nd1->first)));				
						//};

						//change back: remove edge				
						removeEdge(nd2->first, nd1->second);
					};
		
				};
			
		
			};//end of nodes not the same
		};		
				
	};

	//clear cache
	//nodeScoreCache.clear();
};

#ifdef USING_OPEN_MPI

bool NetworkBNLearn::updateNetworkParallel(double & networkScore)
{
	//if(evaluations.size() == 0) return false; //exitErr("Failed to evaluate any network changes!\n");

	multimap<double, pair<unsigned int, unsigned int> >::const_reverse_iterator eval;
	eval = evaluations.rbegin();

	//skip any NaNs
	while(eval != evaluations.rend() && !(eval->first*0 == 0))	eval++;
	
	double initEval = 0;
	unsigned int updateType;

	//pick between equivalent edge changes randomly or will create a bias in the search
	set<pair<unsigned int, unsigned int> > possibleChanges; //parentID, childID
	unsigned int possChangesSize = possibleChanges.size();
	unsigned int bestParentID = 0;
	unsigned int bestChildID = 0;
	set<pair<unsigned int, unsigned int> >::const_iterator pc = possibleChanges.begin();

	double bestEval = 0;
	bool bestFound = false;
	bool noMoreValidEvals = false;

	//if(eval != evaluations.rend() && eval->first > networkScore + 1e-8) //do not want to flip to an equivalent network possible already searched so add numerical accuarcy constant
	//{
	while(!noMoreValidEvals && !bestFound && eval != evaluations.rend())
	{
		//set initial possible edge change
		if(possChangesSize == 0)
		{
			if(eval->first > 1e-8)
			{
				initEval = eval->first;
				bestEval = initEval;
				possibleChanges.insert(make_pair(eval->second.first, eval->second.second));

				eval++;

				while(eval != evaluations.rend() && fabs(eval->first - initEval) < 1e-8)
				{
					possibleChanges.insert(make_pair(eval->second.first, eval->second.second));			
					eval++;
				};
			}
			else
			{
				//return false; //all subsequent edge changes are worse
				bestEval = 0;
				bestParentID = 0;
				bestChildID = 0;				
				noMoreValidEvals = true;
			};		
		};

		if(!noMoreValidEvals)
		{
			possChangesSize = possibleChanges.size();
			pc = possibleChanges.begin();		

			if(possChangesSize == 1) 
			{				
				bestParentID = pc->first;
				bestChildID = pc->second;
			}
			else
			{
				unsigned int ranNo = rand() % possChangesSize;
				for(unsigned int i = 1; i <= ranNo; ++i)
				{
					pc++;
				};
				bestParentID = pc->first;
				bestChildID = pc->second;
			};

			//if(updateNetworkEdgeLoopCheck(bestParentID, bestChildID, updateType))
			if(checkLoopForUpdateNetworkEdgeParallel(bestParentID, bestChildID, updateType))	
			{
				bestFound = true;
				//networkScore += initEval; //update score
				//evaluateAllEdgeChangesDiffsUpdate(bestParentID, bestChildID, updateType);			
				//return true;
			};
		
			if(possChangesSize == 1)
			{
				possibleChanges.clear();				
			}
			else
			{			
				possibleChanges.erase(pc);
			};

			possChangesSize = possibleChanges.size();
		};
	};
	
	if(!bestFound)
	{
		bestEval = 0;
		bestParentID = 0;
		bestChildID = 0;	
	};

	updateBestEvalParallel(bestEval, bestParentID, bestChildID, updateType);

	if(updateType == 0) return false;

	networkScore += bestEval;//initEval;
	evaluateAllEdgeChangesDiffsUpdateParallel(bestParentID, bestChildID, updateType);
	
	return true;
};

//! Calculates score differences in changing the edges, in a more efficient way for bnlearn network.
void NetworkBNLearn::evaluateAllEdgeChangesDiffsParallel()
{
	//cache all network node evals
	double origScore = getScoreFromCache();

	double scoreDiff;
	
	if(freeClearMemory) multimap<double, pair<unsigned int, unsigned int> >().swap(evaluations); else evaluations.clear();

	int noCalcs = allNetworkNodes.size();
	noCalcs = noCalcs*noCalcs - allNetworkNodes.size();
                                                                                                                               
	int numToProcDiv = noCalcs/processSize;
	int numToProcRemainer = noCalcs - numToProcDiv*processSize;

	//startNo = 0;
	//endNo = 0;
	int arraySize = numToProcDiv; //must be the same for all processes
	if(numToProcRemainer > 0) arraySize++;
	
	/*if(processRank < numToProcRemainer)
	{
		startNo = (numToProcDiv+1)*processRank + 1;
		endNo = (numToProcDiv+1)*(processRank+1);
	}
	else
	{
		startNo = (numToProcDiv+1)*numToProcRemainer + numToProcDiv*(processRank - numToProcRemainer) + 1;
		endNo = (numToProcDiv+1)*numToProcRemainer + numToProcDiv*(processRank - numToProcRemainer + 1);
	};*/

	double* evaluationsScoreDiff = new double[arraySize];
	unsigned int* evaluationsFromNode = new unsigned int[arraySize];
	unsigned int* evaluationsToNode = new unsigned int[arraySize];

	currentEvalNumber = 0;
	int recordedEvals = 0;

	//loop thro' pairs of nodes 
	for(map<unsigned int, Node *>::iterator nd1 = allNetworkNodes.begin(); nd1 != allNetworkNodes.end(); ++nd1)
	{
		for(map<unsigned int, Node *>::iterator nd2 = allNetworkNodes.begin(); nd2 != allNetworkNodes.end(); ++nd2)
		{
			if(nd1->first != nd2->first)
			{ 
				currentEvalNumber++;

				//if(currentEvalNumber >= startNo && currentEvalNumber <= endNo)
				if(currentEvalNumber % processSize == processRank)
				{
					if(nd2->second->parentExists(nd1->first)) // nd1 --> nd2 exists
					{
						if((!isWhiteEdge(nd1->first, nd2->first) || (isWhiteEdge(nd1->first, nd2->first) && isWhiteEdge(nd2->first, nd1->first))) //can switch white edge if white in both directions					
							&& edgeAllowed(nd2->first, nd1->first)
							)
						{
							//try REVERSING edge, nd1 <-- nd2	
							removeEdge(nd1->first, nd2->second);
							addEdge(nd2->second, nd1->second);
					
							scoreDiff = calcScore2DifferentNodes(origScore, nd1->first, nd2->first) - origScore;
							if(scoreDiff > 0)
							{
								//evaluations.insert(make_pair(scoreDiff, make_pair(nd2->first, nd1->first)));						
								evaluationsScoreDiff[recordedEvals] = scoreDiff;
								evaluationsFromNode[recordedEvals] = nd2->first;
								evaluationsToNode[recordedEvals] = nd1->first;	
							}
							else
							{
								evaluationsScoreDiff[recordedEvals] = 0;
							};				

							recordedEvals++;
							//change back: now put back the edge, nd1 --> nd2 
							removeEdge(nd2->first, nd1->second);
							addEdge(nd1->second, nd2->second);
						};
					}
					else if(nd1->second->parentExists(nd2->first))  // nd1 <-- nd2 exists
					{
						if(!isWhiteEdge(nd2->first, nd1->first))
						{
							//try DELETING edge, nd2 --> nd1
							removeEdge(nd2->first, nd1->second);
							scoreDiff = calcScore2DifferentNodes(origScore, nd1->first, 0) - origScore;					
							if(scoreDiff > 0)
							{
								//evaluations.insert(make_pair(scoreDiff, make_pair(nd2->first, nd1->first)));
								evaluationsScoreDiff[recordedEvals] = scoreDiff;
								evaluationsFromNode[recordedEvals] = nd2->first;
								evaluationsToNode[recordedEvals] = nd1->first;								
							}
							else
							{
								evaluationsScoreDiff[recordedEvals] = 0;
							};	
					
							recordedEvals++;
							//change back: now put back the edge, nd2 --> nd1 	
							addEdge(nd2->second, nd1->second);		
						};
					}
					else
					{			
						//try ADDING edge nd1 <-- nd2
						if(edgeAllowed(nd2->first, nd1->first))
						{							
							addEdge(nd2->second, nd1->second);
				
							scoreDiff = calcScore2DifferentNodes(origScore, nd1->first, 0) - origScore;					
							if(scoreDiff > 0)
							{
								//evaluations.insert(make_pair(scoreDiff, make_pair(nd2->first, nd1->first)));				
								evaluationsScoreDiff[recordedEvals] = scoreDiff;
								evaluationsFromNode[recordedEvals] = nd2->first;
								evaluationsToNode[recordedEvals] = nd1->first;								
							}
							else
							{
								evaluationsScoreDiff[recordedEvals] = 0;
							};	
						
							recordedEvals++;
							//change back: remove edge				
							removeEdge(nd2->first, nd1->second);
						};
					};
			
				};//end check if calculating

			};//end of nodes not the same
		};		
				
	};
		
	while(recordedEvals < arraySize)
	{
		evaluationsScoreDiff[recordedEvals] = 0;
		recordedEvals++;
	};

	unsigned int totalArraySize = processSize*arraySize;

	double* evaluationsScoreDiffAll = new double[totalArraySize];
	unsigned int* evaluationsFromNodeAll = new unsigned int[totalArraySize];
	unsigned int* evaluationsToNodeAll = new unsigned int[totalArraySize];

	//gather and broadcast all results to all processes                                                                                                                                                                  
	MPI_Allgather(&evaluationsScoreDiff[0], arraySize, MPI_DOUBLE, &evaluationsScoreDiffAll[0], arraySize, MPI_DOUBLE, MPI_COMM_WORLD);
   
	MPI_Allgather(&evaluationsFromNode[0], arraySize, MPI_UNSIGNED, &evaluationsFromNodeAll[0], arraySize, MPI_UNSIGNED, MPI_COMM_WORLD);

	MPI_Allgather(&evaluationsToNode[0], arraySize, MPI_UNSIGNED, &evaluationsToNodeAll[0], arraySize, MPI_UNSIGNED, MPI_COMM_WORLD);

	//now put in usual structure in all processes
	for(int res = 0; res < totalArraySize; ++res)
	{
		if(evaluationsScoreDiffAll[res] > 0)
		{
			evaluations.insert(make_pair(evaluationsScoreDiffAll[res], make_pair(evaluationsFromNodeAll[res], evaluationsToNodeAll[res])));
		};
	};

	delete evaluationsScoreDiff;
	delete evaluationsFromNode;
	delete evaluationsToNode;
	delete evaluationsScoreDiffAll;
	delete evaluationsFromNodeAll;
	delete evaluationsToNodeAll;
};

//! Updates for the best eval across all the different processes for Open MPI parallel processing.
void NetworkBNLearn::updateBestEvalParallel(double & bestEval, unsigned int & bestParentID, unsigned int & bestChildID, unsigned int & updateType)
{
	BestEval bestEvalType;
	bestEvalType.eval = bestEval;
	bestEvalType.parentID = bestParentID;
	bestEvalType.childID = bestChildID;

	//MPI_Aint 
	MPI_Datatype oldtypes[3];
	oldtypes[0] = MPI_DOUBLE;
	oldtypes[1] = MPI_UNSIGNED;
	oldtypes[2] = MPI_UNSIGNED;

	MPI_Datatype ntype1;
	int blocklens[3];
	blocklens[0] = 1;
	blocklens[1] = 1;
	blocklens[2] = 1;

	MPI_Aint indices[3];
	MPI_Aint doubLoc, unintLoc1, unintLoc2;
	MPI_Get_address(&bestEvalType.eval, &doubLoc);
	MPI_Get_address(&bestEvalType.parentID, &unintLoc1);
	MPI_Get_address(&bestEvalType.childID, &unintLoc2);

	indices[0] = 0;
	indices[1] = unintLoc1 - doubLoc;
	indices[2] = unintLoc2 - doubLoc;

	MPI_Type_create_struct(3, blocklens, indices, oldtypes, &ntype1);

	MPI_Type_commit(&ntype1);

	BestEval bestEvalTypeRecv;
	MPI_Status status;

	//send and receive the best evals
	if(processRank != 0)
	{
		MPI_Ssend(&bestEvalType, 1, ntype1, 0, 123, MPI_COMM_WORLD);
	}
	else
	{
		double bestSplitter = 0;
		double splitter;
		for(int i = 1; i < processSize; i++)
		{
			MPI_Recv(&bestEvalTypeRecv, 1, ntype1, i, 123, MPI_COMM_WORLD, &status);
			
			if(bestEvalTypeRecv.eval > 1e-8)
			{
				if(bestEvalTypeRecv.eval > bestEvalType.eval)
				{
					bestEvalType.eval = bestEvalTypeRecv.eval;
					bestEvalType.parentID = bestEvalTypeRecv.parentID;
					bestEvalType.childID = bestEvalTypeRecv.childID;
				}
				else if(bestEvalTypeRecv.eval == bestEvalType.eval)
				{
					//if a draw for the best edge change choose one at random, so each is equally likely to be chosen
					if(bestSplitter==0) bestSplitter = rand()/(double)RAND_MAX;
					splitter = rand()/(double)RAND_MAX;;
					if(splitter > bestSplitter)
					{
						bestSplitter = splitter;
						bestEvalType.eval = bestEvalTypeRecv.eval;
						bestEvalType.parentID = bestEvalTypeRecv.parentID;
						bestEvalType.childID = bestEvalTypeRecv.childID;
					};
				};
			};

		};

	};

	//broadcast the best update to all processes
	MPI_Bcast(&bestEvalType, 1, ntype1, 0, MPI_COMM_WORLD);

	bestEval = bestEvalType.eval;
	bestParentID = bestEvalType.parentID;
	bestChildID = bestEvalType.childID;

	if(bestEvalType.eval > 1e-8)
	{
		//now do the edge update
		//setup nodes	
		Node * parentNode;	
		Node * childNode;
		bool updateParent = false;

		parentNode = getNetworkNode(bestEvalType.parentID);
		childNode = getNetworkNode(bestEvalType.childID);

		//make change, first node is parent and the second the child
		if(edgeExists(bestEvalType.parentID, bestEvalType.childID)) //parent --> child
		{
			//if the edges exists then delete it
			removeEdge(bestEvalType.parentID, childNode);
			updateType = 1;		
		}
		else if(edgeExists(bestEvalType.childID, bestEvalType.parentID)) // parent <-- child
		{
			//if edge the other way exists then reverse the edge
			reverseEdge(childNode, parentNode);
			updateParent = true;
			updateType = 2;		
		}
		else
		{
			//no edge exists so add a new one
			addEdge(parentNode, childNode);
			updateType = 3;		
		};

		setScoreBitInCache(bestEvalType.childID, childNode->calcUpdatedScoreBit(this));//calcScoreBit(this));
		if(updateParent)
		{		
			setScoreBitInCache(bestEvalType.parentID, parentNode->calcUpdatedScoreBit(this));//calcScoreBit(this));
		};
	}
	else
	{
		updateType = 0;
	};

	//broadcast seed, ensure all processes have same random seed, may be different in eval splitting 
	unsigned int randomSeed = rand();
	MPI_Bcast(&randomSeed, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	srand(randomSeed);
};

//! Check if updating this edge would result in a loop.
bool NetworkBNLearn::checkLoopForUpdateNetworkEdgeParallel(const unsigned int & parentNo, const unsigned int & childNo, unsigned int & updateType)
{
	//setup nodes	
	Node * parentNode;	
	Node * childNode;
	bool updateParent = false;

	parentNode = getNetworkNode(parentNo);
	childNode = getNetworkNode(childNo);

	//make change, first node is parent and the second the child
	if(edgeExists(parentNo, childNo)) //parent --> child
	{
		//if the edges exists then delete it
		removeEdge(parentNo, childNode);
		updateType = 1;	
	}
	else if(edgeExists(childNo, parentNo)) // parent <-- child
	{
		//if edge the other way exists then reverse the edge
		reverseEdge(childNode, parentNode);
		updateParent = true;
		updateType = 2;	
	}
	else
	{
		//no edge exists so add a new one
		addEdge(parentNode, childNode);
		updateType = 3;	
	};

	bool wouldHaveLoop = hasLoop();

	//undo change
	if(updateType==1) addEdge(parentNode, childNode);
	else if(updateType==2) reverseEdge(parentNode, childNode);
	else removeEdge(parentNo, childNode);

	return !wouldHaveLoop; //is ok to update
};

#endif

double NetworkBNLearn::evaluateEdgeChangeDiff(const double & networkScore, const unsigned int & nodeID1, const unsigned int & nodeID2)
{
	double scoreDiff = 0;
	Node * node1 = getNetworkNode(nodeID1);
	Node * node2 = getNetworkNode(nodeID2);

	if(node2->parentExists(nodeID1)) // nd1 --> nd2 exists
	{
		if((!isWhiteEdge(nodeID1, nodeID2) || (isWhiteEdge(nodeID1, nodeID2) && isWhiteEdge(nodeID2, nodeID1))) //can switch white edge if white in both directions					
			&& edgeAllowed(nodeID2, nodeID1)
			)
		{
			//try REVERSING edge, nd1 <-- nd2	
			removeEdge(nodeID1, node2);
			addEdge(node2, node1);
			
			scoreDiff = calcScore2DifferentNodes(networkScore, nodeID1, nodeID2) - networkScore;
			
			//change back: now put back the edge, nd1 --> nd2 
			removeEdge(nodeID2, node1);
			addEdge(node1, node2);
		};
	}
	else if(getNetworkNode(nodeID1)->parentExists(nodeID2))  // nd1 <-- nd2 exists
	{
		if(!isWhiteEdge(nodeID2, nodeID1))
		{
			//try DELETING edge, nd2 --> nd1
			removeEdge(nodeID2, node1);

			scoreDiff = calcScore2DifferentNodes(networkScore, nodeID1, 0) - networkScore;					
					
			//change back: now put back the edge, nd2 --> nd1 	
			addEdge(node2, node1);		
		};
	}
	else
	{			
		//try ADDING edge nd1 <-- nd2
		if(edgeAllowed(nodeID2, nodeID1))
		{							
			addEdge(node2, node1);
				
			scoreDiff = calcScore2DifferentNodes(networkScore, nodeID1, 0) - networkScore;					
			
			//change back: remove edge				
			removeEdge(nodeID2, node1);
		};
		
	};

	return scoreDiff;
};

//! Adds one edge evaluation
void NetworkBNLearn::addEdgeChangeDiffOneEdge(const double & networkScore, const unsigned int & nodeFromID, const unsigned int & nodeToID)
{
	double scoreDiff = evaluateEdgeChangeDiff(networkScore, nodeFromID, nodeToID);
	
	if(scoreDiff > 0)
	{
		evaluations.insert(make_pair(scoreDiff, make_pair(nodeFromID, nodeToID)));		
	};
};

//! Adds new evaluations based one one node changing a parent.
void NetworkBNLearn::addEdgeChangeDiffOneNode(const double & networkScore, const unsigned int & nodeID)
{
	double scoreDiff = 0;

	for(map<unsigned int, Node *>::iterator nd1 = allNetworkNodes.begin(); nd1 != allNetworkNodes.end(); ++nd1)
	{
		if(nd1->first != nodeID) scoreDiff = evaluateEdgeChangeDiff(networkScore, nodeID, nd1->first);
		else scoreDiff = 0;

		if(scoreDiff > 0)
		{
			evaluations.insert(make_pair(scoreDiff, make_pair(nd1->first, nodeID)));			
		};
	
	};

	map<unsigned int, Node *> parents = getNetworkNode(nodeID)->getParents();

	for(map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end(); ++pa)
	{
		scoreDiff = evaluateEdgeChangeDiff(networkScore, pa->first, nodeID);
		if(scoreDiff > 0)
		{
			evaluations.insert(make_pair(scoreDiff, make_pair(nodeID, pa->first)));			
		}
	
	};
};

//! Updates evaluations for some node pairings are keeps ones that have not changed.
void NetworkBNLearn::evaluateAllEdgeChangesDiffsUpdate(const unsigned int & parentID, const unsigned int & childID, unsigned int & updateType)
{
	//updateType, 1 = delete, 2 = reverse, 3 = add
	double origScore = getScoreFromCache();

	//remove evaluations that need to be updated firstly
	unsigned int nodeUpdate1 = childID;
	unsigned int nodeUpdate2 = 0;

	if(updateType==2) nodeUpdate2 = parentID;

	for(multimap<double, pair<unsigned int, unsigned int> >::iterator eval = evaluations.begin(); eval != evaluations.end(); )
	{
		if(eval->second.second == nodeUpdate1 || 
		   eval->second.second == nodeUpdate2 ||
		   (getNetworkNode(childID)->parentExists(eval->second.second) && (eval->second.first == nodeUpdate1 || eval->second.first == nodeUpdate2)) ||
		   (updateType == 2 && getNetworkNode(parentID)->parentExists(eval->second.second) && (eval->second.first == nodeUpdate1 || eval->second.first == nodeUpdate2)) ||
		   (updateType == 1 && eval->second.first == childID && eval->second.second == parentID)		   
		   )
		{			
			evaluations.erase(eval++);
		}		
		else
		{
			++eval;
		};

	};


	//now add in updated evaluations
	addEdgeChangeDiffOneNode(origScore, nodeUpdate1);
	if(nodeUpdate2 != 0) addEdgeChangeDiffOneNode(origScore, nodeUpdate2); 
	
	//If an edge was deleted then need to reevaluate the edge also in the other direction (as would have been taken previously from a reverse edge and childNode may have changed parents since then)
	if(updateType == 1) addEdgeChangeDiffOneEdge(origScore, childID, parentID); 
};

#ifdef USING_OPEN_MPI

//! Adds one edge evaluation in parallel using Open MPI.
void NetworkBNLearn::addEdgeChangeDiffOneEdgeParallel(const double & networkScore, const unsigned int & nodeFromID, const unsigned int & nodeToID)
{
	//if(currentEvalNumber >= startNo && currentEvalNumber <= endNo) //only calc some for parallel process
	if(currentEvalNumber % processSize == processRank)
	{
		double scoreDiff = evaluateEdgeChangeDiff(networkScore, nodeFromID, nodeToID);
	
		if(scoreDiff > 0)
		{
			evaluations.insert(make_pair(scoreDiff, make_pair(nodeFromID, nodeToID)));			
			//evaluationsScoreDiff[recordedEvals] = scoreDiff;			
		}
		/*else
		{
			evaluationsScoreDiff[recordedEvals] = 0;
		};				

		evaluationsFromNode[recordedEvals] = nodeFromID;
		evaluationsToNode[recordedEvals] = nodeToID;	
		recordedEvals++;*/
	};

	currentEvalNumber++;
};

//! Adds new evaluations based one one node changing a parent in parallel using Open MPI.
void NetworkBNLearn::addEdgeChangeDiffOneNodeParallel(const double & networkScore, const unsigned int & nodeID)
{
	double scoreDiff = 0;

	for(map<unsigned int, Node *>::iterator nd1 = allNetworkNodes.begin(); nd1 != allNetworkNodes.end(); ++nd1)
	{
		if(nd1->first != nodeID)
		{
			//if(currentEvalNumber >= startNo && currentEvalNumber <= endNo) //only calc some for parallel process
			if(currentEvalNumber % processSize == processRank)
			{
				scoreDiff = evaluateEdgeChangeDiff(networkScore, nodeID, nd1->first);
		
				if(scoreDiff > 0)
				{
					evaluations.insert(make_pair(scoreDiff, make_pair(nd1->first, nodeID)));				
					//evaluationsScoreDiff[recordedEvals] = scoreDiff;					
				}
				/*else
				{
					evaluationsScoreDiff[recordedEvals] = 0;
				};				

				evaluationsFromNode[recordedEvals] = nd1->first;
				evaluationsToNode[recordedEvals] = nodeID;	
				recordedEvals++;*/
			};

			currentEvalNumber++;
		};
	};

	map<unsigned int, Node *> parents = getNetworkNode(nodeID)->getParents();

	for(map<unsigned int, Node *>::const_iterator pa = parents.begin(); pa != parents.end(); ++pa)
	{
		//if(currentEvalNumber >= startNo && currentEvalNumber <= endNo) //only calc some for parallel process
		if(currentEvalNumber % processSize == processRank)
		{
			scoreDiff = evaluateEdgeChangeDiff(networkScore, pa->first, nodeID);
		
			if(scoreDiff > 0)
			{
				evaluations.insert(make_pair(scoreDiff, make_pair(nodeID, pa->first)));			
				//evaluationsScoreDiff[recordedEvals] = scoreDiff;				
			}
			/*else
			{
				evaluationsScoreDiff[recordedEvals] = 0;
			};				

			evaluationsFromNode[recordedEvals] = nodeID;
			evaluationsToNode[recordedEvals] = pa->first;	
			recordedEvals++;*/
		};

		currentEvalNumber++;
	};
};

//! Updates evaluations for some node pairings are keeps ones that have not changed in parallel using Open MPI.
void NetworkBNLearn::evaluateAllEdgeChangesDiffsUpdateParallel(const unsigned int & parentID, const unsigned int & childID, unsigned int & updateType)
{
	//updateType, 1 = delete, 2 = reverse, 3 = add
	double origScore = getScoreFromCache();

	//remove evaluations that need to be updated firstly
	unsigned int nodeUpdate1 = childID;
	unsigned int nodeUpdate2 = 0;

	if(updateType==2) nodeUpdate2 = parentID;

	for(multimap<double, pair<unsigned int, unsigned int> >::iterator eval = evaluations.begin(); eval != evaluations.end(); )
	{
		if(eval->second.second == nodeUpdate1 || 
		   eval->second.second == nodeUpdate2 ||
		   (getNetworkNode(childID)->parentExists(eval->second.second) && (eval->second.first == nodeUpdate1 || eval->second.first == nodeUpdate2)) ||
		   (updateType == 2 && getNetworkNode(parentID)->parentExists(eval->second.second) && (eval->second.first == nodeUpdate1 || eval->second.first == nodeUpdate2)) ||
		   (updateType == 1 && eval->second.first == childID && eval->second.second == parentID)		   
		   )
		{
			evaluations.erase(eval++);
		}		
		else
		{
			++eval;
		};

	};

	
	currentEvalNumber = 1;

	/////////////////////////////Start parallel evaluations
	//now add in updated evaluations
	addEdgeChangeDiffOneNodeParallel(origScore, nodeUpdate1);
	if(nodeUpdate2 != 0) addEdgeChangeDiffOneNodeParallel(origScore, nodeUpdate2); 
	
	//If an edge was deleted then need to reevaluate the edge also in the other direction (as would have been taken previously from a reverse edge and childNode may have changed parents since then)
	if(updateType == 1)
	{		
		addEdgeChangeDiffOneEdgeParallel(origScore, childID, parentID); 
	};
	/////////////////////////////End parallel evaluations

};
#endif

//! Calculates score differences in changing the edges, in a more efficient way for bnlearn network.
void NetworkBNLearn::evaluateAllEdgeChanges()
{
	//now done elsewhere on the fly

	////cache all network node evals
	//double origScore = getScoreFromCache();

	////clear list of previous evaluations of edge changes
	//evaluations.clear();
	//double score;
	//
	////loop thro' pairs of nodes 
	//for(map<unsigned int, Node *>::iterator nd1 = allNetworkNodes.begin(); nd1 != allNetworkNodes.end(); ++nd1)
	//{
	//	for(map<unsigned int, Node *>::iterator nd2 = allNetworkNodes.begin(); nd2 != allNetworkNodes.end(); ++nd2)
	//	{
	//		if(nd1->first != nd2->first)
	//		{ 
	//
	//		if(nd2->second->parentExists(nd1->first)) // nd1 --> nd2 exists
	//		{
	//			if((!isWhiteEdge(nd1->first, nd2->first) || (isWhiteEdge(nd1->first, nd2->first) && isWhiteEdge(nd2->first, nd1->first))) //can switch white edge if white in both directions					
	//				&& edgeAllowed(nd2->first, nd1->first)
	//				)
	//			{
	//				//try REVERSING edge, nd1 <-- nd2	
	//				removeEdge(nd1->first, nd2->second);
	//				addEdge(nd2->second, nd1->second);

	//				if(!hasLoop())
	//				{						
	//					score = calcScore2DifferentNodes(origScore, nd1->first, nd2->first);
	//					evaluations.insert(make_pair(score, make_pair(nd2->first, nd1->first)));						
	//				};

	//				//change back: now put back the edge, nd1 --> nd2 
	//				removeEdge(nd2->first, nd1->second);
	//				addEdge(nd1->second, nd2->second);
	//			};
	//		}
	//		else if(nd1->second->parentExists(nd2->first))  // nd1 <-- nd2 exists
	//		{
	//			if(!isWhiteEdge(nd2->first, nd1->first))
	//			{
	//				//try DELETING edge, nd2 --> nd1
	//				removeEdge(nd2->first, nd1->second);
	//				score = calcScore2DifferentNodes(origScore, nd1->first, 0);					
	//				evaluations.insert(make_pair(score, make_pair(nd2->first, nd1->first)));
	//				
	//				//change back: now put back the edge, nd2 --> nd1 	
	//				addEdge(nd2->second, nd1->second);		
	//			};
	//		}
	//		else
	//		{			
	//			//try ADDING edge nd1 <-- nd2
	//			if(edgeAllowed(nd2->first, nd1->first))
	//			{							
	//				addEdge(nd2->second, nd1->second);
	//
	//				if(!hasLoop())
	//				{					
	//					score = calcScore2DifferentNodes(origScore, nd1->first, 0);					
	//					evaluations.insert(make_pair(score, make_pair(nd2->first, nd1->first)));				
	//				};

	//				//change back: remove edge				
	//				removeEdge(nd2->first, nd1->second);
	//			};
	//	
	//		};
	//		
	//	
	//		};//end of nodes not the same
	//	};		
	//			
	//};

};


//! Update network according to the two given nodes, returns if reversed or not.
void Network::updateNetworkEdge(unsigned int & parentNo, unsigned int & childNo)
{
	//setup nodes	
	Node * parentNode;	
	Node * childNode;
	
	parentNode = getNetworkNode(parentNo);
	childNode = getNetworkNode(childNo);

	//make change, first node is parent and the second the child
	if(edgeExists(parentNo, childNo)) //parent --> child
	{
		//if the edges exists then delete it
		removeEdge(parentNo, childNode);
	}
	else if(edgeExists(childNo, parentNo)) // parent <-- child
	{
		//if edge the other way exists then reverse the edge
		reverseEdge(childNode, parentNode);
		//parentNode->setDefaultInitPrior();
		setScoreBitInCache(parentNo, parentNode->calcUpdatedScoreBit(this));
	}
	else
	{
		//no edge exists so add a new one
		addEdge(parentNode, childNode);
	};

	//childNode->setDefaultInitPrior();
	//calculatePriorsAndPosteriorsDeal();

	setScoreBitInCache(childNo, childNode->calcUpdatedScoreBit(this));
};

#include <time.h>


//! Update network according to the two given nodes, returns if has loop or not.
bool NetworkBNLearn::updateNetworkEdgeLoopCheck(const unsigned int & parentNo, const unsigned int & childNo, unsigned int & updateType)
{
	
	//setup nodes	
	Node * parentNode;	
	Node * childNode;
	bool updateParent = false;

	parentNode = getNetworkNode(parentNo);
	childNode = getNetworkNode(childNo);

	//make change, first node is parent and the second the child
	if(edgeExists(parentNo, childNo)) //parent --> child
	{
		//if the edges exists then delete it
		removeEdge(parentNo, childNode);
		updateType = 1;	
	}
	else if(edgeExists(childNo, parentNo)) // parent <-- child
	{
		//if edge the other way exists then reverse the edge
		reverseEdge(childNode, parentNode);
		updateParent = true;
		updateType = 2;	
	}
	else
	{
		//no edge exists so add a new one
		addEdge(parentNode, childNode);
		updateType = 3;		
	};

	if(hasLoop())
	{
		//undo change
		if(updateType==1) addEdge(parentNode, childNode);
		else if(updateType==2) reverseEdge(parentNode, childNode);
		else removeEdge(parentNo, childNode);

		return false;
	};

	//childNode->setDefaultInitPrior();
	//calculatePriorsAndPosteriorsDeal();

	setScoreBitInCache(childNo, childNode->calcUpdatedScoreBit(this));//calcScoreBit(this));
	if(updateParent)
	{
		//parentNode->setDefaultInitPrior();
		setScoreBitInCache(parentNo, parentNode->calcUpdatedScoreBit(this));//calcScoreBit(this));
	};

	return true; //edge has updated
};

//! Update the network according to evaluations of edges.
bool Network::updateNetwork(double & networkScore)
{
	if(evaluations.size() == 0) return false; //exitErr("Failed to evaluate any network changes!\n");

	multimap<double, pair<unsigned int, unsigned int> >::const_reverse_iterator eval;
	eval = evaluations.rbegin();

	//skip any NaNs
	while(eval != evaluations.rend() && !(eval->first*0 == 0))
	{
		eval++;
	};

	double initEval = 0;
	
	//pick between equivalent edge changes randomly or will create a bias in the search
	set<pair<unsigned int, unsigned int> > possibleChanges; //parentID, childID

	if(eval != evaluations.rend() && eval->first > networkScore + 1e-8) //do not want to flip to an equivalent network possible already searched so add numerical accuarcy constant
	{
		initEval = eval->first;

		possibleChanges.insert(make_pair(eval->second.first, eval->second.second));

		eval++;

		while(eval != evaluations.rend() && fabs(eval->first - initEval) < 1e-8)
		{
			possibleChanges.insert(make_pair(eval->second.first, eval->second.second));
			eval++;
		};

		unsigned int possChangesSize = possibleChanges.size();
		unsigned int bestParentID, bestChildID;
		set<pair<unsigned int, unsigned int> >::const_iterator pc = possibleChanges.begin();

		if(possChangesSize == 1) 
		{
			bestParentID = pc->first;
			bestChildID = pc->second;
		}
		else
		{
			unsigned int ranNo = rand() % possChangesSize;
			for(unsigned int i = 1; i <= ranNo; ++i)
			{
				pc++;
			};
			bestParentID = pc->first;
			bestChildID = pc->second;
		};

		networkScore = initEval;
		updateNetworkEdge(bestParentID, bestChildID);
		return true;
	};

	return false;
};

//! Update the network according to evaluations of edges.
bool NetworkBNLearn::updateNetwork(double & networkScore)
{
	if(evaluations.size() == 0) return false; //exitErr("Failed to evaluate any network changes!\n");

	multimap<double, pair<unsigned int, unsigned int> >::const_reverse_iterator eval;
	eval = evaluations.rbegin();

	//skip any NaNs
	while(eval != evaluations.rend() && !(eval->first*0 == 0))	eval++;
	
	double initEval = 0;
	unsigned int updateType;

	//pick between equivalent edge changes randomly or will create a bias in the search
	set<pair<unsigned int, unsigned int> > possibleChanges; //parentID, childID
	unsigned int possChangesSize = possibleChanges.size();
	unsigned int bestParentID, bestChildID;
	set<pair<unsigned int, unsigned int> >::const_iterator pc = possibleChanges.begin();

	//if(eval != evaluations.rend() && eval->first > networkScore + 1e-8) //do not want to flip to an equivalent network possible already searched so add numerical accuarcy constant
	//{
	while(eval != evaluations.rend())
	{
		//set initial possible edge change
		if(possChangesSize == 0)
		{
			if(eval->first > 1e-8)
			{
				initEval = eval->first;
				possibleChanges.insert(make_pair(eval->second.first, eval->second.second));
			}
			else
			{
				return false; //all subsequent edge changes are worse
			};
			
			eval++;

			while(eval != evaluations.rend() && fabs(eval->first - initEval) < 1e-8)
			{
				possibleChanges.insert(make_pair(eval->second.first, eval->second.second));			
				eval++;
			};
		};

		possChangesSize = possibleChanges.size();
		pc = possibleChanges.begin();
		

		if(possChangesSize == 1) 
		{
			bestParentID = pc->first;
			bestChildID = pc->second;
		}
		else
		{
			unsigned int ranNo = rand() % possChangesSize;
			for(unsigned int i = 1; i <= ranNo; ++i)
			{
				pc++;
			};
			bestParentID = pc->first;
			bestChildID = pc->second;
		};

		if(updateNetworkEdgeLoopCheck(bestParentID, bestChildID, updateType))
		{
			networkScore += initEval; //update score
			evaluateAllEdgeChangesDiffsUpdate(bestParentID, bestChildID, updateType);
			
			return true;
		};
		
		if(possChangesSize == 1)
		{
			possibleChanges.clear();			
		}
		else
		{			
			possibleChanges.erase(pc);
		};

		possChangesSize = possibleChanges.size();
	};
	
	return false;
};

//! Calculates priors and posteriors for dela network.
void NetworkDeal::calculatePriorsAndPosteriors()
{ 
	clearPriorsAndPosteriors();
	setupLocalPriProbs("");
	calcDiscreteJointPriorDist(); //includes calculating prior local prob. dists
	calcDiscreteMasterPriors(); //set up master prior first, loop thro' nodes
	calcDiscreteLocalParaPris();  
	calcDiscreteLocalParaPosts();

	calcCtsJointPriorDist(); //includes calculating prior local prob. dists
	calcCtsMasterPriors();   //set up master prior first, loop thro' nodes
	calcCtsLocalParaPris();
	calcCtsLocalParaPosts();
};

//! Sets up master priors for each discrete node.
void NetworkDeal::calcDiscreteMasterPriors()
{
	//loop thro' nodes and calculate the master priors
	for(set<DiscreteDealNode *>::iterator nd = discreteDealNodes.begin(); nd != discreteDealNodes.end(); ++nd)
	{
		(*nd)->calcDiscreteMasterPrior(this);
	};

};

//! Sets up master priors for each discrete node.
void NetworkDeal::calcCtsMasterPriors()
{
	//loop thro' nodes and calculate the master priors
	for(set<CtsDealNode *>::iterator nd = ctsDealNodes.begin(); nd != ctsDealNodes.end(); ++nd)
	{
		(*nd)->calcCtsMasterPrior(this);
	};
};

//! Sets up the discrete local parameter priors.
void NetworkDeal::calcDiscreteLocalParaPris()
{
	//loop thro' nodes and calculate the Local Parameter priors
	for(set<DiscreteDealNode *>::iterator nd = discreteDealNodes.begin(); nd != discreteDealNodes.end(); ++nd)
	{
		(*nd)->calcDiscreteLocalParameterPrior();
	};
};

//! Sets up the discrete local parameter priors.
void NetworkDeal::calcDiscreteLocalParaPosts()
{
	//loop thro' nodes and calculate the Local Parameter posteriors
	for(set<DiscreteDealNode *>::iterator nd = discreteDealNodes.begin(); nd != discreteDealNodes.end(); ++nd)
	{
		(*nd)->calcDiscreteLocalParameterPost();
	};

};

//! Outputs the network to file.
void Network::outputNetwork(string & filename)
{
	ofstream outNet(filename.c_str());

	//output nodes
	for(map<unsigned int, Node *>::const_iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		if(!nd->second->getIsFactorChildNode())
		{
			outNet << nd->second->getDisplayName() << "\n";
		};
	};

	//output edges
	for(map<unsigned int, Node *>::const_iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{		
		nd->second->outputEdgeNames(outNet);		
	};

	outNet.close();
};

//! Returns PDAG string for network.
void Network::getPDAGData(const string & PDAGNetwork)
{
	if(freeClearMemory) list<pair<string, string> >().swap(nodesSwappableInit); else nodesSwappableInit.clear();
	string networkPDAGString = PDAGNetwork;
	bool foundNode;
	
	string nodeParentsStr;
	string nodeStr;
	string parentsStr, aParent;
	set<string> someParents;
	
	size_t pos1, pos2, posBar, posDir, posV, posDiv, posColon;

	do{

		foundNode = false;
		pos1 = networkPDAGString.find_first_of('(');
		pos2 = networkPDAGString.find_first_of(')');

		if(pos1 < networkPDAGString.length() && pos1 != string::npos && pos2 < networkPDAGString.length() && pos2 != string::npos && (pos1 + 1) < pos2)
		{
			foundNode = true;
			nodeParentsStr = networkPDAGString.substr(pos1+1, pos2-pos1-1);
			//extract node name and parents, so A|B:C gives A with B and C as parents
			posBar = nodeParentsStr.find_first_of('|');
			posDir = nodeParentsStr.find_first_of('^');
			posV = nodeParentsStr.find_first_of('<');

			//get node name
			if(posDir < nodeParentsStr.length() && posDir != string::npos) posDiv = posDir;
			else if(posV < nodeParentsStr.length() && posV != string::npos) posDiv = posV;
			else posDiv = posBar;

			if(posDiv < nodeParentsStr.length() && posDiv != string::npos) nodeStr = nodeParentsStr.substr(0, posDiv);
			else nodeStr = nodeParentsStr;

			someParents.clear();

			//check if parents exist
			if(posDiv < nodeParentsStr.length() && posDiv != string::npos)
			{
				parentsStr = nodeParentsStr.substr(posDiv+1);

				do{
					//foundParent = false;					
					
					posColon = parentsStr.find_first_of(':');

					if(posColon < parentsStr.length() && posColon != string::npos) aParent = parentsStr.substr(0, posColon);
					else aParent = parentsStr;

					//add edge 
					someParents.insert(aParent);

					if(posColon < parentsStr.length() && posColon != string::npos)
					{
						//chop off previous found parent node
						parentsStr = parentsStr.substr(posColon+1);
					}
					else 
						parentsStr = "";

				}while(parentsStr != "");

			};

			if(someParents.size() > 0 && posBar < nodeParentsStr.length() && posBar != string::npos)
			{
				for(set<string>::const_iterator sp = someParents.begin(); sp != someParents.end(); ++sp)
				{
					nodesSwappableInit.push_back(make_pair(nodeStr, *sp));				
				};	
			};

			//chop off previous found node
			networkPDAGString = networkPDAGString.substr(pos2+1);
		};

	}while(foundNode);

};

//! Outputs the list of equivalent networks to file.
unsigned int Network::outputEquivalentNetworks(const bool & returnList, list<string> & eqNetworks, const string & nodeNamePre)
{
	
	unsigned int noEquivNets = 1;
	ofstream outEquivNets;
	
	string PDAGNetwork = getPDAGNetworkString();
	
	string PDAGNetworkEquiv;

	string nodeNamePre2;

	if(!returnList)
	{
		outEquivNets.open(nodeNamePre.c_str());
		outEquivNets << "Partially_directed_graph_(PDAG):\n";
		outEquivNets << PDAGNetwork << "\n";
		outEquivNets << "Network(s):\n";
		outEquivNets << getNetworkString(0) << "\n";
		nodeNamePre2 = "";
	}
	else
	{
		eqNetworks.push_back(getNetworkString(0));
		nodeNamePre2 = nodeNamePre;
	};
	
	getPDAGData(PDAGNetwork);

	//list<pair<unsigned int, unsigned int> > nodesSwappable;
	list<bool> swapped;

	list<bool>::const_iterator sw;
	list<unsigned int>::const_iterator ns2; 

	unsigned int nodeNumber1; 
	unsigned int nodeNumber2;
	bool isDiscrete1;
	bool isDiscrete2;
	set<string> someNodes;
	string nodeName1, nodeName2;


	list<unsigned int> nodesSwappable1 = list<unsigned int>();
	list<unsigned int> nodesSwappable2 = list<unsigned int>(); 

	//add if dis --> someNode
	for(list<pair<string, string> >::iterator nsi = nodesSwappableInit.begin(); nsi != nodesSwappableInit.end(); ++nsi)
	{
		someNodes.clear();
		nodeName1 = nodeNamePre2 + nsi->first;//avoid temporary parameter literal problem in below function by putting name in variable
		nodeName2 = nodeNamePre2 + nsi->second;
		nodeNumber1 = allNodeData->getNodeDataNumber(nodeName1);
		nodeNumber2 = allNodeData->getNodeDataNumber(nodeName2);
		isDiscrete1 = getNetworkNode(nodeNumber1)->getIsDiscreteNode();
		isDiscrete2 = getNetworkNode(nodeNumber2)->getIsDiscreteNode();

		if(((isDiscrete1 && isDiscrete2) || (!isDiscrete1 && !isDiscrete2))
			&& 
				!isWhiteEdge(nodeNumber1, nodeNumber2) && !isWhiteEdge(nodeNumber2, nodeNumber1)
				&& edgeAllowed(nodeNumber1, nodeNumber2) && edgeAllowed(nodeNumber2, nodeNumber1)
				
			)
		{
			//add in order they are currently in the network			
			if(edgeExists(nodeNumber1, nodeNumber2))
			{
				nodesSwappable1.push_back(nodeNumber1);
				nodesSwappable2.push_back(nodeNumber2);
			}
			else
			{
				nodesSwappable1.push_back(nodeNumber2);
				nodesSwappable2.push_back(nodeNumber1);
			};
			
			swapped.push_back(false);			
		};
	};

	

	string node1, node2;

	list<bool>::iterator swp = swapped.begin();

	bool equivNetworkIsOK;

	//output equiv networks network
	while(swp != swapped.end())
	{
		
		//update swapped record of swappable nodes
		if(!*swp) *swp = true;
		else
		{
			while(swp != swapped.end() && *swp)
			{
				*swp = false;
				++swp;				
			};

			if(swp != swapped.end() && !*swp)
			{
				*swp = true;
				swp = swapped.begin();
			};
		};
	
		if(swp == swapped.end()) break;

		sw = swapped.begin();

		ns2 = nodesSwappable2.begin();

		for(list<unsigned int>::const_iterator ns1 = nodesSwappable1.begin(); ns1 != nodesSwappable1.end(); ++ns1, ++ns2, ++sw)
		{
			if(!*sw)
			{
				nodeNumber1 = *ns1; 
				nodeNumber2 = *ns2; 
			}
			else
			{
				nodeNumber2 = *ns1; 
				nodeNumber1 = *ns2;
			};			

			//update actual network object
			if(edgeExists(nodeNumber2, nodeNumber1)) removeEdge(nodeNumber2, nodeNumber1);
			if(!edgeExists(nodeNumber1, nodeNumber2)) addEdge(nodeNumber1, nodeNumber2);
		};
		
		//check if there are any loops
		equivNetworkIsOK = !hasLoop();

		//check if the v-strs are the same
		if(equivNetworkIsOK)
		{
			PDAGNetworkEquiv = getPDAGNetworkString();
			equivNetworkIsOK = (PDAGNetworkEquiv == PDAGNetwork); //should result in the same PDAG network
		};

		if(equivNetworkIsOK)
		{
			if(returnList) eqNetworks.push_back(getNetworkString(0));
			else outEquivNets << getNetworkString(0) << "\n";
			noEquivNets++;
		};
		
	};

	ns2 = nodesSwappable2.begin();

	//put back network as it was found
	for(list<unsigned int>::const_iterator ns1 = nodesSwappable1.begin(); ns1 != nodesSwappable1.end(); ++ns1, ++ns2)
	{
		nodeNumber1 = *ns1; 
		nodeNumber2 = *ns2; 

		//update actual network object
		if(edgeExists(nodeNumber2, nodeNumber1)) removeEdge(nodeNumber2, nodeNumber1);
		if(!edgeExists(nodeNumber1, nodeNumber2)) addEdge(nodeNumber1, nodeNumber2);
	};

	if(freeClearMemory) list<pair<string, string> >().swap(nodesSwappableInit); else nodesSwappableInit.clear();

	if(!returnList) outEquivNets.close();

	return noEquivNets;
};

//! Calculates the recall and precision of the otherNetwork compared to this one as original.
void Network::calcRecallPrecision(Network * otherNetwork, double & recallOrPrecision, const string & preNodeNameOther, const bool & isRecall)
{
	//get list of edges which may be reversed
	list<pair<string, string> > nodesSwappableList;
	string PDAGStr;
	
	if(!isRecall)
	{
		PDAGStr = otherNetwork->getPDAGNetworkString();
		otherNetwork->getPDAGData(PDAGStr);
		nodesSwappableList = otherNetwork->getNodesSwappable();
	}
	else {
		PDAGStr = getPDAGNetworkString();
		getPDAGData(PDAGStr);
		nodesSwappableList = getNodesSwappable();
	};

	map<unsigned int, Node *> someParents;
	string aNodeName, paNodeName;
	double edgesFound = 0;
	double totalEdges = 0;

	set<string> nodesSwappable;
	set<string>::const_iterator ns;

	for(list<pair<string, string> >::iterator nsi = nodesSwappableList.begin(); nsi != nodesSwappableList.end(); ++nsi)
	{
		nodesSwappable.insert(nsi->first + " " + nsi->second);		
	};

	//calc recall
	for(map<unsigned int, Node *>::const_iterator ann = allNetworkNodes.begin(); ann != allNetworkNodes.end(); ++ann)
	{
		aNodeName = preNodeNameOther + ann->second->getDisplayName();
		someParents = ann->second->getParents();
		for(map<unsigned int, Node *>::const_iterator pa = someParents.begin(); pa != someParents.end(); ++pa)
		{
			paNodeName = preNodeNameOther + pa->second->getDisplayName();			
			
			if(otherNetwork->edgeExistsInit(paNodeName, aNodeName)) edgesFound++;
			else
			{
				//check if may be other direction
				ns = nodesSwappable.find(ann->second->getDisplayName() + " " + pa->second->getDisplayName());
				if(ns == nodesSwappable.end()) ns = nodesSwappable.find(pa->second->getDisplayName() + " " + ann->second->getDisplayName()); 
				if(ns != nodesSwappable.end() && otherNetwork->edgeExistsInit(aNodeName, paNodeName)) edgesFound++; 
			};

			totalEdges++;
		};
	};
	
	//add to prev recall as may be taking ave later
	if(totalEdges != 0) recallOrPrecision += edgesFound / totalEdges;
	else recallOrPrecision += 1;
};

//! Copies missingness pattern to the other network data from the corresponding network nodes of this network. 
void Network::copyMissingness(Network * otherNetwork, const string & preNodeNameOther, const bool & reverse)
{
	Node * otherNode;
	string aNodeName;
	for(map<unsigned int, Node *>::const_iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		if(reverse)
		{
			aNodeName = nd->second->getDisplayName(); //.substr(preNodeNameOther.length());
		}
		else
		{
			aNodeName = preNodeNameOther + nd->second->getDisplayName();
		};
		otherNode = otherNetwork->getNetworkNode(allNodeData->getNodeDataNumber(aNodeName));
		nd->second->copyMissingness(otherNode);
	};
};

//! Outputs the network edges for igraph.
void Network::outputNetworkIGraphEdges(string & filename)
{
	//calc the signif of each edge
	calcEdgeSignifs();

	ofstream edgeFile(filename.c_str());

	//header of file
	edgeFile << "from to chisq\n";

	for(map<unsigned int, Node *>::const_iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{
		nd->second->outputEdges(edgeFile, nd->first);
	};

	edgeFile.close();
};

//! Outputs the network nodes for igraph.
void Network::outputNetworkIGraphNodes(string & filename)
{
	ofstream nodesFile(filename.c_str());

	//header of file
	nodesFile << "id name type fileno\n";
	string type;

	for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{	
		if(!nd->second->getIsFactorChildNode())
		{
			if(nd->second->getIsFactorNode()) type = "f"; else type = nd->second->getNodeType();
			nodesFile << nd->first << " " << nd->second->getDisplayName() << " " << type << " " << nd->second->getFileNo() << "\n";
		};
	};

	nodesFile.close();
};

//! Outputs network node data.
void Network::outputNodeData(const string & nodeDataFilePrefix, const bool & outputBedFile, bool & wasDisData, bool & wasCtsData, bool & wasSNPDisData, bool & wasSNPCtsData, bool & wasSNPData, const unsigned int & startIndivNo, const unsigned int & endIndivNo, const unsigned int & jobNo)
{
	if(outputBedFile)
	{
		/*if(startIndivNo != 0)
		{
			outErr("Outputting SNP data as a .bed file for a subset of individuals is not available!");
			exit(1);
		};*/

		outputSNPNodeData(nodeDataFilePrefix, outputBedFile, wasSNPData, startIndivNo, endIndivNo); //may set IDs for other files
	}
	else
	{
		outputDiscreteNodeData(nodeDataFilePrefix, outputBedFile, true, wasDisData, wasSNPDisData, startIndivNo, endIndivNo, jobNo);

		outputCtsNodeData(nodeDataFilePrefix, outputBedFile, true, wasCtsData, wasSNPCtsData, startIndivNo, endIndivNo, jobNo);	
	};

	outputDiscreteNodeData(nodeDataFilePrefix, outputBedFile, false, wasDisData, wasSNPDisData, startIndivNo, endIndivNo, jobNo);

	outputCtsNodeData(nodeDataFilePrefix, outputBedFile, false, wasCtsData, wasSNPCtsData, startIndivNo, endIndivNo, jobNo);	

};

//! Outputs discrete network node data.
void Network::outputDiscreteNodeData(const string & nodeDataFilePrefix, const bool & outputBedFile, const bool & outSNPData, bool & wasDisData, bool & wasSNPDisData, const unsigned int & startIndivNo, const unsigned int & endIndivNo, const unsigned int & jobNo)
{
	list<string> header;
	list<DiscreteNode *> outDiscreteNodes;
	list< list<unsigned int>::const_iterator > dataLevels;
	list< list<bool>::const_iterator > dataMissing;

	//set up data iterators
	for(map<unsigned int, DiscreteNode *>::const_iterator dn = discreteNodes.begin(); dn != discreteNodes.end(); ++dn)
	{
		//if(!outputBedFile || !dn->second->getIsSNPNode())
		if((!outSNPData && !dn->second->getIsSNPNode()) || (outSNPData && dn->second->getIsSNPNode()))
		{
			header.push_back(dn->second->getDisplayName());
			outDiscreteNodes.push_back(dn->second);
			dataLevels.push_back(dn->second->getDiscreteData()->values.begin());
			dataMissing.push_back(dn->second->getDiscreteData()->missingValues.begin());
		};
	};

	list<CtsNode *> outCtsFactorNodes;
	list< list<double>::const_iterator > dataFactorValues;
	map<unsigned int, list<double>::const_iterator > dataChildFactorValues; //node ID, iterator for data of child factor data
	map<unsigned int, list<double>::const_iterator >::iterator dcfv;
	list<unsigned int> childFactorNodes;

	//set up data iterators for factor nodes
	for(map<unsigned int, CtsNode *>::const_iterator dn = ctsNodes.begin(); dn != ctsNodes.end(); ++dn)
	{

		if(((!outSNPData && !dn->second->getIsSNPNode()) || (outSNPData && dn->second->getIsSNPNode())) && (dn->second->getIsFactorNode()))
		{
			header.push_back(dn->second->getDisplayName());		
			dataMissing.push_back(dn->second->getCtsData()->missingValues.begin());	
			outCtsFactorNodes.push_back(dn->second);
			dataFactorValues.push_back(dn->second->getCtsData()->values.begin());

			childFactorNodes = allNodeData->getFactorDataGroup(dn->first);
			for(list<unsigned int>::const_iterator cfn = childFactorNodes.begin(); cfn != childFactorNodes.end(); ++cfn)
			{
				dataChildFactorValues[*cfn] = getNetworkNode(*cfn)->getCtsData()->values.begin();
			};
		};
	};

	if(header.empty()) return; //if no data then do not write to file

	if(outSNPData) wasSNPDisData = true;
	else wasDisData = true;

	unsigned int length;
	if(!discreteNodes.empty()) length = discreteNodes.begin()->second->getDiscreteData()->values.size();
	else length = ctsNodes.begin()->second->getCtsData()->values.size();

	//open up file and write data
	string discreteFileName;
	if(outSNPData) discreteFileName = nodeDataFilePrefix+"-SNPs-discrete.dat";
	else discreteFileName = nodeDataFilePrefix+"-discrete.dat";

	ofstream outDisFile(discreteFileName.c_str());

	list<string> nameOfIDs = allNodeData->getNameOfIDs();
	unsigned int noIDS = nameOfIDs.size();

	if(jobNo <= 1) //do not write header if split into jobs unless at start
	{ 
		//write name of IDs		
		for(list<string>::const_iterator nid = nameOfIDs.begin(); nid != nameOfIDs.end(); ++nid)
		{
			outDisFile << *nid << " ";		
		};

		//write header
		for(list<string>::const_iterator hd = header.begin(); hd != header.end(); )
		{
			outDisFile << *hd;
			++hd;
			if(hd != header.end()) outDisFile << " "; 
		};

		outDisFile << "\n";
	};

	list<string> idData = allNodeData->getDataIDs();
	list<string>::const_iterator id = idData.begin();
	list<DiscreteNode *>::const_iterator dnd; 
	list<CtsNode *>::const_iterator cnd; 
	CtsNode * headFactorNode;
	list< list<bool>::const_iterator >::iterator md;
	bool writeRow;
	unsigned int factorLevel;
	unsigned int factorCount;

	unsigned int count = 1;
	do{

		writeRow = startIndivNo == 0 || (count >= startIndivNo && count <= endIndivNo);
	
		//write IDs before data
		for(unsigned int i = 1; i <= noIDS; i++, ++id)
		{
			if(writeRow) outDisFile << *id << " ";
		};
		
		dnd = outDiscreteNodes.begin();
		md = dataMissing.begin();

		//output one row of data
		for(list< list<unsigned int>::const_iterator >::iterator dl = dataLevels.begin(); dl != dataLevels.end(); )
		{		
			if(writeRow)
			{
				if(!*(*md)) outDisFile << (*dnd)->getLevelName(*(*dl));
				else outDisFile << "NA";
			};

			//advance iterator to the next row
			++(*dl); ++(*md);
		
			//next node/data and column
			++dl; ++md;
			if(writeRow && (dl != dataLevels.end() || !outCtsFactorNodes.empty())) outDisFile << " "; 			
			++dnd;
		};
		
		cnd = outCtsFactorNodes.begin();

		//output the factor node data (if any)
		for(list< list<double>::const_iterator >::iterator dfv = dataFactorValues.begin(); dfv != dataFactorValues.end(); )
		{
			
			headFactorNode = *cnd;
			
			//head IDs of child factor nodes
			childFactorNodes = allNodeData->getFactorDataGroup(headFactorNode->getNodeID());

			//get factor level
			factorLevel = 0;
			if(*(*dfv) == 1) factorLevel = 1;
			factorCount = 1;

			for(list<unsigned int>::const_iterator cfn = childFactorNodes.begin(); cfn != childFactorNodes.end(); ++cfn)
			{
				factorCount++;
				//get data for this child factor node
				dcfv = dataChildFactorValues.find(*cfn);
				if(dcfv != dataChildFactorValues.end())
				{
					if(*(dcfv->second) == 1) factorLevel = factorCount;
					++(dcfv->second); //advance data for next row
				};

			};

			if(writeRow)
			{	
				if(!*(*md)) outDisFile << headFactorNode->getFactorLevelName(factorLevel);
				else outDisFile << "NA";	
			};

			//advance iterator to the next row
			++(*dfv); ++(*md);
		
			//next node/data and column
			++dfv; ++md;
			if(writeRow && dfv != dataFactorValues.end()) outDisFile << " "; 			

			++cnd; //next head factor node
		};
		if(writeRow) outDisFile << "\n";

		count++;
	}while(count <= length);

	outDisFile.close();
};

//! Outputs cts network node data.
void Network::outputCtsNodeData(const string & nodeDataFilePrefix, const bool & outputBedFile, const bool & outSNPData, bool & wasCtsData, bool & wasSNPCtsData, const unsigned int & startIndivNo, const unsigned int & endIndivNo, const unsigned int & jobNo)
{
	list<string> header;
	list<CtsNode *> outCtsNodes;
	list< list<double>::const_iterator > dataValues;
	list< list<bool>::const_iterator > dataMissing;

	//set up data iterators
	for(map<unsigned int, CtsNode *>::const_iterator dn = ctsNodes.begin(); dn != ctsNodes.end(); ++dn)
	{
		//if(!outputBedFile || !dn->second->getIsSNPNode())
		if(((!outSNPData && !dn->second->getIsSNPNode()) || (outSNPData && dn->second->getIsSNPNode())) && !(dn->second->getIsFactorNode() || dn->second->getIsFactorChildNode()))
		{
			header.push_back(dn->second->getDisplayName());
			outCtsNodes.push_back(dn->second);
			dataValues.push_back(dn->second->getCtsData()->values.begin());
			dataMissing.push_back(dn->second->getCtsData()->missingValues.begin());
		};
	};

	if(header.empty()) return; //if no data then do not write to file

	if(outSNPData) wasSNPCtsData = true;
	else wasCtsData = true;

	unsigned int length = ctsNodes.begin()->second->getCtsData()->values.size(); //number of indivs

	string ctsFileName;
	if(outSNPData) ctsFileName = nodeDataFilePrefix+"-SNPs-cts.dat";
	else ctsFileName = nodeDataFilePrefix+"-cts.dat";

	ofstream outCtsFile(ctsFileName.c_str());

	list<string> nameOfIDs = allNodeData->getNameOfIDs();
	unsigned int noIDS = nameOfIDs.size();


	if(jobNo <= 1) //if split into jobs just write header for first job
	{
		//write name of IDs	
		for(list<string>::const_iterator nid = nameOfIDs.begin(); nid != nameOfIDs.end(); ++nid)
		{
			outCtsFile << *nid << " ";	
		};

		//write header
		for(list<string>::const_iterator hd = header.begin(); hd != header.end(); )
		{
			outCtsFile << *hd;	
			++hd;
			if(hd != header.end()) outCtsFile << " "; 
		};

		outCtsFile << "\n";
	};

	list<string> idData = allNodeData->getDataIDs();
	list<string>::const_iterator id = idData.begin();
	list< list<bool>::const_iterator >::iterator md;
	bool writeRow;

	unsigned int count = 1;
	do{

		writeRow = startIndivNo == 0 || (count >= startIndivNo && count <= endIndivNo);
	
		//write IDs before data
		for(unsigned int i = 1; i <= noIDS; i++, ++id)
		{
			if(writeRow) outCtsFile << *id << " ";
		};
		

		md = dataMissing.begin();

		//output one row of data
		for(list< list<double>::const_iterator >::iterator dl = dataValues.begin(); dl != dataValues.end(); )
		{
			if(writeRow)
			{
				if(!*(*md)) outCtsFile << *(*dl);
				else outCtsFile << "NA";
			};

			//advance iterator to the next row
			++(*dl); ++(*md);
		
			//next node/data and column
			++dl; ++md;
			if(writeRow && dl != dataValues.end()) outCtsFile << " "; 			
			
		};
		if(writeRow) outCtsFile << "\n";

		count++;
	}while(count <= length);

	outCtsFile.close();
};

//! Write binary data.
void writeBedGenotype(ofstream & psBedFile, unsigned int & psBitCount, int & aBit, const bool & allele1, const bool & allele2)
{

	aBit = aBit >> 1; //shift bits to the right
	if(allele1) aBit = aBit | 128; 
	aBit = aBit >> 1;  //shift bits to the right
	if(allele2) aBit = aBit | 128;
	
	psBitCount += 2;

	//write to file if byte is finished
	if(psBitCount == 8)
	{
		//write to file
		char buffer[1];
		buffer[0] = aBit;
		psBedFile.write(buffer, 1);
		psBitCount = 0;
		aBit = 0;
	};

};

//! Write last byte and resets the Bit and bit counter.
void writeLastByteBeforeNextSNP(ofstream & psBedFile, unsigned int & psBitCount, int & aBit)
{
	if(psBitCount == 0) return; //last byte may already be written

	//shift right bits over
	while(psBitCount < 8)
	{
		aBit = aBit >> 1; //shift bits to the right
		psBitCount++;
	};
	
	//write to file
	char buffer[1];
	buffer[0] = aBit;
	psBedFile.write(buffer, 1);
	psBitCount = 0;
	aBit = 0;
};

//! Outputs discrete and cts network node data to .bed/.bim/.fam files.
void Network::outputSNPNodeData(const string & nodeDataFilePrefix, const bool & outputBedFile, bool & wasSNPData, const unsigned int & startIndivNo, const unsigned int & endIndivNo)
{
	wasSNPData = false;

	//check there is some SNP data
	for(map<unsigned int, Node *>::const_iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
	{		
		wasSNPData = wasSNPData || nd->second->getIsSNPNode(); 
	};

	if(!wasSNPData) return; //no SNP data
	
	//Output .bim/.fam/.bed files
	string outputBimFilename = nodeDataFilePrefix + ".bim";
	string outputFamFilename = nodeDataFilePrefix + ".fam";
	string outputBedFilename = nodeDataFilePrefix + ".bed";
	
	ofstream outBimFile(outputBimFilename.c_str());
	ofstream outBedFile(outputBedFilename.c_str(), ios::binary);

	//write out initial binary pedigree file bytes, first 2 byte are magic numbers the third to indicate SNP major (subjects x SNPs)
	unsigned int psBitCount = 0;
	int aBit = 0;
	char buffer[3];	
	buffer[0] = 108;
	buffer[1] = 27;
	buffer[2] = 1;
	outBedFile.write(buffer, 3);

	CtsData * ctsData;
	list<bool>::const_iterator mi;
	unsigned int basePos = 10;
	unsigned int indivNo;

	//Output cts SNPs to .bed file
	for(map<unsigned int, CtsNode *>::const_iterator cn = ctsNodes.begin(); cn != ctsNodes.end(); ++cn)
	{
		if(cn->second->getIsSNPNode())
		{
			//do bim file first
			outBimFile << 1 << " " << cn->second->getDisplayName() << " " <<  0  << " " << basePos << " 1 2\n";

			ctsData = cn->second->getCtsData();
			mi = ctsData->missingValues.begin();
			indivNo = 1;
			for(list<double>::const_iterator cd = ctsData->values.begin(); cd != ctsData->values.end(); ++cd, ++mi, ++indivNo)
			{
				if(((startIndivNo == 0) || indivNo >= startIndivNo) && ((endIndivNo == 0) || indivNo <= endIndivNo))
				{
					if(*mi)	writeBedGenotype(outBedFile, psBitCount, aBit, 1, 0); //missing
					else
					{
						if(*cd == 0) writeBedGenotype(outBedFile, psBitCount, aBit, 0, 0);
						else if(*cd == 1) writeBedGenotype(outBedFile, psBitCount, aBit, 0, 1);
						else writeBedGenotype(outBedFile, psBitCount, aBit, 1, 1);
					};
				};
			};

			writeLastByteBeforeNextSNP(outBedFile, psBitCount, aBit);
		};
	};

	//set up which is what
	unsigned int geno00Level;
	unsigned int geno01Level;
	unsigned int geno11Level;
	
	DiscreteData * disData;
	//Output dis SNPs to .bed file
	for(map<unsigned int, DiscreteNode *>::const_iterator dn = discreteNodes.begin(); dn != discreteNodes.end(); ++dn)
	{
		if(dn->second->getIsSNPNode())
		{
			//do bim file first
			outBimFile << 1 << " " << dn->second->getDisplayName() << " " <<  0  << " " << basePos << " 1 2\n";			
			basePos += 10;

			geno00Level = dn->second->getLevelNo("0");
			geno01Level = dn->second->getLevelNo("1");
			geno11Level = dn->second->getLevelNo("2");

			disData = dn->second->getDiscreteData();
			mi = disData->missingValues.begin();
			indivNo = 1;
			for(list<unsigned int>::const_iterator dd = disData->values.begin(); dd != disData->values.end(); ++dd, ++mi, ++indivNo)
			{
				if(((startIndivNo == 0) || indivNo >= startIndivNo) && ((endIndivNo == 0) || indivNo <= endIndivNo))
				{
					if(*mi)	writeBedGenotype(outBedFile, psBitCount, aBit, 1, 0); //missing
					else
					{
						if(*dd == geno00Level) writeBedGenotype(outBedFile, psBitCount, aBit, 0, 0);
						else if(*dd == geno01Level) writeBedGenotype(outBedFile, psBitCount, aBit, 0, 1);
						else writeBedGenotype(outBedFile, psBitCount, aBit, 1, 1);
					};
				};
			};

			writeLastByteBeforeNextSNP(outBedFile, psBitCount, aBit);
		};
	};

	outBedFile.close();
	outBimFile.close();

	//output basic .fam file
	unsigned int noIndivs = allNodeData->getAmountOfData();

	//output IDs if given
	list<string> nameOfIDs = allNodeData->getNameOfIDs();
	unsigned int noIDs = nameOfIDs.size();

	list<string> idData = allNodeData->getDataIDs();
	list<string>::const_iterator id = idData.begin();
	bool outputIndiv;

	ofstream outFamFile(outputFamFilename.c_str());
	for(unsigned int i = 1; i <= noIndivs; ++i)
	{
		outputIndiv = ((startIndivNo == 0) || i >= startIndivNo) && ((endIndivNo == 0) || i <= endIndivNo);
		
		if(noIDs==2)
		{
			//write IDs before data
			for(unsigned int i2 = 1; i2 <= noIDs; i2++, ++id)
			{
				if(outputIndiv) outFamFile << *id << " ";
			};

			if(outputIndiv) outFamFile << "0 0 1 0\n";
		}
		else
		{
			if(outputIndiv) outFamFile << i << " " << i << " 0 0 1 0\n";
		};
		
	};
	outFamFile.close();

	//set IDs for discrete or cts data to match the fam file if not defined
	if(noIDs==0)
	{
		list<string> dataIDs; 
		list<string> nameOfIDs;
		nameOfIDs.push_back("FAMID");
		nameOfIDs.push_back("ID");

		for(unsigned int i = 1; i <= noIndivs; ++i)
		{
			dataIDs.push_back(toString(i));
			dataIDs.push_back(toString(i));
		};

		allNodeData->setNameOfIDs(nameOfIDs);
		allNodeData->setIDs(dataIDs);
	};
};

//! Outputs priors for deal network for testing.
//void NetworkDeal::outputPriorsTEST()
//{	
//	//output the joint priors first
//	cout << "Prior Local Probability Distributions:\n";
//	cout << "======================================\n\n";
//	
//
//	//output the joint priors first
//	cout << "Joint Priors:\n";
//	cout << "=============\n\n";
//	cout << "Discrete Joint Priors:\n";
//
//	unsigned int covarSize;
//
//	cout << "\nContinuous Joint Priors:\n";
//	for(map<unsigned int, CtsMultiDistJoint *>::const_iterator cjp = ctsJointProbDist.begin(); cjp != ctsJointProbDist.end(); ++cjp)
//	{
//		cout << " Group no. " << cjp->first << ":\n";
//		
//		cout << "  Means:\n   ";
//		for(map<unsigned int, double>::const_iterator m = cjp->second->means.begin(); m != cjp->second->means.end(); ++m)
//		{
//			cout << m->first << ": " << m->second << " ";
//		};
//
//		cout << "\n  Covariance matrix:\n";
//		covarSize = cjp->second->variances.size();
//
//		for(unsigned int nd1 = 1; nd1 <= covarSize; ++nd1)
//		{
//			cout << "   ";
//			for(unsigned int nd2 = 1; nd2 <= covarSize; ++nd2)
//			{
//				cout << cjp->second->getCovariance(nd1, nd2) << " ";
//			};
//
//			cout << "\n";
//		};
//		
//	};
//
//};

//! Outputs the priors for the deal network.
void NetworkDeal::outputPriors(string & fileName)
{
	setPriorPostNames();

	ofstream priorsFile(fileName.c_str());

	//output the joint priors first
	priorsFile << "Prior Local Probability Distributions:\n";
	priorsFile << "======================================\n\n";
	for(set<DiscreteDealNode *>::const_iterator dnd = discreteDealNodes.begin(); dnd != discreteDealNodes.end(); ++dnd)
	{
		(*dnd)->outputPriorLocalProbDists(priorsFile);
	};

	for(set<CtsDealNode *>::const_iterator cnd = ctsDealNodes.begin(); cnd != ctsDealNodes.end(); ++cnd)
	{
		(*cnd)->outputPriorLocalProbDists(priorsFile);
	};

	//output the joint priors first
	priorsFile << "Joint Priors:\n";
	priorsFile << "=============\n\n";
	priorsFile << "Discrete Joint Priors:\n";
	priorsFile << " DISCRETE PARENTS: ";
	for(map<unsigned int, DiscreteNode *>::const_iterator dnd0 = discreteNodes.begin(); dnd0 != discreteNodes.end();)
	{
		priorsFile << dnd0->second->getDisplayName();
		++dnd0;
		if(dnd0 != discreteNodes.end()) priorsFile << ":";
	};
	priorsFile << "\n";
	if(discreteJointNames.size() != discreteJointProbDist.size()) exitErr("Problem setting up discrete prior names!");
	map<unsigned int, string>::const_iterator djn = discreteJointNames.begin();
	for(map<unsigned int, double>::const_iterator djp = discreteJointProbDist.begin(); djp != discreteJointProbDist.end(); ++djp, ++djn)
	{
		priorsFile << djn->second << ": " << djp->second << "\n"; 
	};

	unsigned int covarSize;

	priorsFile << "\nContinuous Joint Priors:\n";
	priorsFile << " DISCRETE PARENTS: ";
	for(map<unsigned int, DiscreteNode *>::const_iterator dnd1 = discreteNodes.begin(); dnd1 != discreteNodes.end();)
	{
		priorsFile << dnd1->second->getDisplayName();
		++dnd1;
		if(dnd1 != discreteNodes.end()) priorsFile << ":";
	};
	priorsFile << "\n";
	map<unsigned int, string>::const_iterator mnm = meanNames.begin();
	for(map<unsigned int, CtsMultiDistJoint *>::const_iterator cjp = ctsJointProbDist.begin(); cjp != ctsJointProbDist.end(); ++cjp)
	{
		priorsFile << cjp->second->name << ":\n";
		
		priorsFile << "  Means:\n   ";

		if(cjp->second->means.size() != meanNames.size()) exitErr("Problem setting up continuous joint prior means!");
		mnm = meanNames.begin();
		for(map<unsigned int, double>::const_iterator m = cjp->second->means.begin(); m != cjp->second->means.end(); ++m, ++mnm)
		{
			priorsFile << mnm->second << ": " << m->second << " ";
		};

		priorsFile << "\n  Covariance matrix:\n";
		covarSize = cjp->second->variances.size();

		for(unsigned int nd1 = 1; nd1 <= covarSize; ++nd1)
		{
			priorsFile << "   ";
			for(unsigned int nd2 = 1; nd2 <= covarSize; ++nd2)
			{
				priorsFile << cjp->second->getCovariance(nd1, nd2) << " ";
			};

			priorsFile << "\n";
		};
		
	};

	//output master priors
	priorsFile << "\nMaster Priors:\n";
	priorsFile << "==============\n\n";

	for(set<DiscreteDealNode *>::const_iterator dnd = discreteDealNodes.begin(); dnd != discreteDealNodes.end(); ++dnd)
	{
		(*dnd)->outputMasterPriors(priorsFile);
	};

	for(set<CtsDealNode *>::const_iterator cnd = ctsDealNodes.begin(); cnd != ctsDealNodes.end(); ++cnd)
	{
		(*cnd)->outputMasterPriors(priorsFile);
	};

	priorsFile << "Local Parameter Priors:\n";
	priorsFile << "==============\n\n";

	for(set<DiscreteDealNode *>::const_iterator dnd = discreteDealNodes.begin(); dnd != discreteDealNodes.end(); ++dnd)
	{
		(*dnd)->outputLocalParameterPriors(priorsFile);
	};

	for(set<CtsDealNode *>::const_iterator cnd = ctsDealNodes.begin(); cnd != ctsDealNodes.end(); ++cnd)
	{
		(*cnd)->outputLocalParameterPriors(priorsFile);
	};

	priorsFile.close();
};

//! Outputs the posteriors for deal network.
void NetworkDeal::outputPosteriors(string & fileName)
{
	setPriorPostNames();

	ofstream posteriorsFile(fileName.c_str());

	posteriorsFile << "Local Parameter Posteriors:\n";
	posteriorsFile << "===========================\n\n";

	for(set<DiscreteDealNode *>::const_iterator dnd = discreteDealNodes.begin(); dnd != discreteDealNodes.end(); ++dnd)
	{
		(*dnd)->outputLocalParameterPosteriors(posteriorsFile);
	};

	for(set<CtsDealNode *>::const_iterator cnd = ctsDealNodes.begin(); cnd != ctsDealNodes.end(); ++cnd)
	{
		(*cnd)->outputLocalParameterPosteriors(posteriorsFile);
	};

	posteriorsFile.close();
};

//! Calculates priors and posteriors for BN learns.
void NetworkBNLearn::calculatePriorsAndPosteriors()
{
	setupDefaultInitPriors();

	calculateDiscretePosteriors();
	calculateCtsPosteriors();
};

//! Calculates discrete posteriors for bnlearn network.
void NetworkBNLearn::calculateDiscretePosteriors()
{
	//done elsewhere 
};

//! Calculates cts posteriors for bnlearn network.
void NetworkBNLearn::calculateCtsPosteriors()
{
	//done elsewhere 
};

//! Simulates data for the network.
void NetworkBNLearn::simulateData(unsigned int & noSims, bool & snpDataIntegers)
{
	
	if(freeClearMemory) list<Node*>().swap(orderedNodes); else orderedNodes.clear();

	clearVisitedNodes();

	unsigned int noVisitedThisLoop = 0;
	unsigned int totalVisited = 0;
	unsigned int totalNodes = allNetworkNodes.size();

	bool updated = false;	

	//try to visit every node if the parents are given (if any), if no loops then can visit every loop
	//output even if a string is detected, but try to do in order first
	do{
		noVisitedThisLoop = 0;

		for(map<unsigned int, Node *>::iterator nd = allNetworkNodes.begin(); nd != allNetworkNodes.end(); ++nd)
		{
			if(!nd->second->getVisited() && nd->second->allParentsVisited())
			{
				orderedNodes.push_back(nd->second);
				nd->second->setVisited(true);
				++noVisitedThisLoop;
			};

		};//end of looping thro' nodes

		if(noVisitedThisLoop == 0)
		{
			exitErr("Loop detected in network while trying to simulate network data!\n");
		};

		totalVisited += noVisitedThisLoop;
	}while(totalVisited < totalNodes);

	//do the simulating of data
	for(unsigned int simNo = 1; simNo <= noSims; ++simNo)
	{
		simulateDataOne(snpDataIntegers);
	};

	updateNetworkMissingData();

	//simulate some more if there was missing data
	pair<unsigned int, unsigned int> missingNotMissing = getNoMissingNotMissing();

	//if there was a problem simulating data, say rare discrete combinations appear where linear reg failed so sim failed, so try again
	if(missingNotMissing.first > 0)
	{
		unsigned int count = 1; //just in case, not possible to simulate any non missing data avoid inf loop
		do{
			simulateDataOne(snpDataIntegers);
			count++;
			updateNetworkMissingData();
		}while(getNoMissingNotMissing().second < noSims && count < 100000);

	};

};

//! Simulates data for the network.
void NetworkBNLearn::simulateDataOne(bool & snpDataIntegers)
{

	//loop thro' nodes in order of dependency
	for(list<Node *>::const_iterator on = orderedNodes.begin(); on != orderedNodes.end(); ++on)
	{
		(*on)->simData(snpDataIntegers);
	};
	
};

//! Sets some data to missing used when sim'ing data.
void NetworkBNLearn::setSomeMissingData(const map<string, unsigned int> & missingFirst, map<string, unsigned int> & missingLast, map<string, double> & missingProb)
{
	Data * theData = 0;
	list<bool> missingValues;
	//bool missing;
	unsigned int n;
	unsigned int firstMissingNo;
	double ranNum;

	for(map<string, unsigned int>::const_iterator mf = missingFirst.begin(); mf != missingFirst.end(); ++mf)
	{
		theData = allNodeData->getNodeData(mf->first);
		n = theData->getAmountOfData();

		for(unsigned int i = 1; i <= n; ++i) missingValues.push_back(i <= mf->second);

		theData->setMissingData(missingValues);
	};

	for(map<string, unsigned int>::const_iterator ml = missingLast.begin(); ml != missingLast.end(); ++ml)
	{
		theData = allNodeData->getNodeData(ml->first);
		n = theData->getAmountOfData();
		firstMissingNo = n - ml->second + 1;

		for(unsigned int i = 1; i <= n; ++i) missingValues.push_back(i >= firstMissingNo);

		theData->setMissingData(missingValues);
	};

	for(map<string, double>::const_iterator mp = missingProb.begin(); mp != missingProb.end(); ++mp)
	{
		theData = allNodeData->getNodeData(mp->first);
		n = theData->getAmountOfData();
		ranNum = (double)rand()/(double)RAND_MAX;

		for(unsigned int i = 1; i <= n; ++i) missingValues.push_back(ranNum <= mp->second);

		theData->setMissingData(missingValues);
	};

	//if(theData != 0) theData->setMissingData(missingValues);

}

//! Adds nodes to a BN learn network.
void NetworkBNLearn::addNode(const unsigned int & nodeNumber)
{
	Data * theData = allNodeData->getNodeData(nodeNumber);

	if(theData->dataType == 1) //discrete
	{
		DiscreteData * disData = static_cast<DiscreteData*>(theData);
		DiscreteBNLearnNode * aNode = new DiscreteBNLearnNode(theData->name, disData);
		allNetworkNodes[nodeNumber] = aNode;
		discreteNodes[nodeNumber] = aNode;
		aNode->setNodeID(nodeNumber);
		discreteBNLearnNodes.insert(aNode);
	}
	else if(theData->dataType == 2) //cts
	{
		CtsData * ctsData = static_cast<CtsData*>(theData);		
		CtsBNLearnNode * aNode = new CtsBNLearnNode(theData->name, ctsData);
		allNetworkNodes[nodeNumber] = aNode;
		ctsNodes[nodeNumber] = aNode;
		aNode->setNodeID(nodeNumber);
		ctsBNLearnNodes.insert(aNode);
	}
	else if(theData->dataType == 3)	addNoDataNode(theData->name, static_cast<NoData*>(theData));
	
};

//! Outputs priors for a BN learn network.
void NetworkBNLearn::outputPriors(string & fileName)
{
	setPriorPostNames();

	ofstream priorsFile(fileName.c_str());

	//output the priors first
	priorsFile << "Priors for nodes are not used in BayesNetty when using a bnlearn type network\n";
	priorsFile << "=============================================================================\n\n";

	priorsFile.close();
}

//! Outputs posteriors for a BN learn network.
void NetworkBNLearn::outputPosteriors(string & fileName)
{
	setPriorPostNames();

	ofstream posteriorsFile(fileName.c_str());

	posteriorsFile << "Posteriors:\n";
	posteriorsFile << "===========\n\n";

	for(set<DiscreteBNLearnNode *>::const_iterator dnd = discreteBNLearnNodes.begin(); dnd != discreteBNLearnNodes.end(); ++dnd)
	{
		(*dnd)->outputPosteriors(posteriorsFile);
	};

	for(set<CtsBNLearnNode *>::const_iterator cnd = ctsBNLearnNodes.begin(); cnd != ctsBNLearnNodes.end(); ++cnd)
	{
		(*cnd)->outputPriorLocalProbDists(posteriorsFile); //these are the posteriors for bnlearn (but called prior local prob dist in deal)
	};

	posteriorsFile.close();
};
