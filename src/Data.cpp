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


/*! \file Data.cpp
    \brief This file contains the source for manipulating SNP data.
    
*/

//using namespace std; // initiates the "std" or "standard" namespace

#include "main.h"
#include "Data.h"
#include <math.h>

//! Sets up the was imputed list.
void Data::setupWasImputed()
{
	if(wasImputed.empty()) 
	{
		for(list<bool>::const_iterator m = missingValues.begin(); m != missingValues.end(); ++m)
		{
			wasImputed.push_back(false);
		};
	}
	else
	{
		for(list<bool>::iterator wi = wasImputed.begin(); wi != wasImputed.end(); ++wi)
		{
			*wi = false;
		};
	}
};

//! Updates missing data as non missing of it was imputed.
void Data::updateImputedDataAsNonMissing()
{
	list<bool>::iterator m = missingValues.begin();
	for(list<bool>::const_iterator wi = wasImputed.begin(); wi != wasImputed.end(); ++wi, ++m)
	{
		if(*wi) *m = false;
	};
};

//! Returns amount of missing data.
unsigned int NetworkMissingData::getNoMissing()
{
	unsigned int noMissing = 0;

	for(list<bool>::const_iterator m = missing.begin(); m != missing.end(); ++m) if(*m) noMissing++;

	return noMissing;
};

//! Adds cts data.
void CtsData::addData(const double & num, const bool & missing)
{
	values.push_back(num);
	missingValues.push_back(missing);
};

//! Returns mean of cts data.
double CtsData::getMean()
{
	double total = 0;
	double count = 0;

	list<bool>::const_iterator mv = missingValues.begin(); 
	for(list<double>::const_iterator v = values.begin(); v != values.end(); ++v, ++mv)
	{
		if(!*mv)
		{
			total += *v;
			count++;
		};

	};
	
	return total/count;
};

//! Returns standard dev. of cts data.
double CtsData::getStDev()
{
	double mean = getMean();
	double total = 0;
	double count = 0;
	double diff;

	list<bool>::const_iterator mv = missingValues.begin(); 
	for(list<double>::const_iterator v = values.begin(); v != values.end(); ++v, ++mv)
	{
		if(!*mv)
		{
			diff = (*v - mean);
			total += diff*diff;
			count++;
		};

	};
	
	return sqrt(total/(count-1));
};

//! Return the level name when the cts data is a factor data (for head factor node).
string CtsData::getFactorLevelName(const unsigned int & lv)
{
	map<unsigned int, string>::const_iterator i = levels.find(lv);
	if(i != levels.end())
	{
		return i->second;
	};

	return "?";
};

//!adds missing data, if already missing then remains missing.
void CtsData::setMissingData(const list<bool> & md)
{
	if(md.size() !=  missingValues.size())
	{
		outErr("Problem setting missing data!");
	};

	list<bool>::iterator mv = missingValues.begin();
	for(list<bool>::const_iterator newmv = md.begin(); newmv != md.end(); ++newmv, ++mv)
	{
		*mv = *mv || *newmv;
	};

};

//! Returns a random value from the complete data.
double CtsData::getRandomCompleteValue()
{

	unsigned int noCompleteValues = completeValues.size();

	if(completeValues.size()==0)
	{
		unsigned int count = 1;		
		list<bool>::const_iterator mv = missingValues.begin(); 
		for(list<double>::const_iterator v = values.begin(); v != values.end(); ++v, ++mv)
		{
			if(!*mv)
			{						
				completeValues[count] = *v;
				count++;
			};
		};

		noCompleteValues = completeValues.size();
	};

	if(noCompleteValues==0) return nan("");

	//get the value randomly
	unsigned int ranNum = rand() % noCompleteValues + 1;
	
	double randomValue = completeValues[ranNum];

	return randomValue;
};

//! Adds discrete data.
void DiscreteData::addData(const string & value, bool & missing)
{
	map<string, unsigned int>::const_iterator l = levels.find(value);

	if(l != levels.end()) values.push_back(l->second);
	else
	{
		if(!missing)
		{
			//get new level
			unsigned int newLevel = levels.size() + 1;
			levels[value] = newLevel;
			values.push_back(newLevel);
		}
		else
		{
			values.push_back(0); //not a valid level because missing data
		};

	};	
	
	missingValues.push_back(missing);
};

//! Adds missing data, if already missing then remains missing.
void DiscreteData::setMissingData(const list<bool> & md)
{
	if(md.size() !=  missingValues.size())
	{
		outErr("Problem setting missing data!");
	};

	list<bool>::iterator mv = missingValues.begin();
	for(list<bool>::const_iterator newmv = md.begin(); newmv != md.end(); ++newmv, ++mv)
	{
		*mv = *mv || *newmv;
	};

};

//! Adds new level to discrete data without actually addding any data.
unsigned int DiscreteData::addLevel(const string & value)
{
	map<string, unsigned int>::const_iterator l = levels.find(value);

	if(l == levels.end())	
	{		
		//get new level
		unsigned int newLevel = levels.size() + 1;
		levels[value] = newLevel;
		return newLevel;
	};	

	return l->second;
};

//! Adds sim'd discrete data.
void DiscreteData::addSimData(const unsigned int & level, const bool & missing)
{
	values.push_back(level);
	missingValues.push_back(missing);
};

//! Returns level name, for outputting priors/posts.
string DiscreteData::getLevelName(const unsigned int & level)
{
	for(map<string, unsigned int>::const_iterator l = levels.begin(); l != levels.end(); ++l)
	{
		if(l->second == level) return l->first;
	};

	outErr("Cannot find level "); out(level); out(" for a discrete node.");
	exitErr("Unable to find level name for a discrete name!");

	return "";
};

//! Returns a random level
unsigned int DiscreteData::getRandomLevel()
{

	unsigned int noLevels = levels.size();

	if(levelProbsCum.size()==0)
	{
		double cumProb = 0;		
		double total = 0;
		map<unsigned int, unsigned int> levCounts;
		list<bool>::const_iterator mv = missingValues.begin(); 
		for(list<unsigned int>::const_iterator v = values.begin(); v != values.end(); ++v, ++mv)
		{
			if(!*mv)
			{						
				levCounts[*v]++;
				total++;
			};
		};

		unsigned int lev;
		for(map<unsigned int, unsigned int>::const_iterator lv = levCounts.begin(); lv != levCounts.end(); )
		{
			cumProb += (double)lv->second/total;
			lev = lv->first;
			++lv;
			if(lv != levCounts.end()) levelProbsCum[cumProb] = lev; else levelProbsCum[1] = lev;
		};

	};

	//get the level randomly
	double ranNum = (double)rand()/(double)RAND_MAX;
	unsigned int randomLevel = 1;

	for(map<double, unsigned int>::const_iterator lvp = levelProbsCum.begin(); lvp != levelProbsCum.end(); 	++lvp)
	{
		if(ranNum <= lvp->first)
		{
			randomLevel = lvp->second;			
			break;
		};
	
	};

	return randomLevel;
};

// Check if node data is added.
bool AllNodeData::isNodeDataAdded(string & dataName)
{
	return (nameNumber.find(dataName) != nameNumber.end());
};

//! Adds a node.
void AllNodeData::addNodeData(Data * data)
{
	map<string, unsigned int>::const_iterator nn = nameNumber.find(data->name);
	
	if(nn != nameNumber.end())
	{
		string mess = "Attempt to add variable \"" + data->name + "\" twice to the data set!";
		exitErr(mess);
	};

	unsigned int nodeNumber = nameNumber.size() + 1;
	allNodeData[nodeNumber] = data;
	nameNumber[data->name] = nodeNumber;
};

//Deletes node data for given id.
void AllNodeData::deleteNodeData(const unsigned int & nodeNumber)
{
	map<unsigned int, Data *>::iterator oldNode = allNodeData.find(nodeNumber);

	if(oldNode == allNodeData.end())
	{
		string mess = "Attempt to delete data but the node could not be found!";
		exitErr(mess);
	};

	map<string, unsigned int>::iterator nn = nameNumber.find(oldNode->second->name);
	if(nn != nameNumber.end()) nameNumber.erase(nn);

	map<unsigned int, list<unsigned int> >::iterator fg = factorGroups.find(oldNode->first);
	if(fg != factorGroups.end()) factorGroups.erase(fg);

	delete oldNode->second;
	allNodeData.erase(oldNode);
};

//! Replaces existing node data, used for updating NoData data.
void AllNodeData::replaceNodeData(Data * data)
{
	map<string, unsigned int>::const_iterator nn = nameNumber.find(data->name);
	
	if(nn == nameNumber.end())
	{
		string mess = "Attempt to replace data \"" + data->name + "\" but it does not exist!";
		exitErr(mess);
	};

	map<unsigned int, Data *>::iterator oldNode = allNodeData.find(nn->second);

	if(oldNode == allNodeData.end())
	{
		string mess = "Attempt to replace data \"" + data->name + "\" but the node could not be found!";
		exitErr(mess);
	};

	delete oldNode->second;

	allNodeData[nn->second] = data;
};

//! Gets new node no without adding any data.
unsigned int AllNodeData::getNewNodeNo(const string & nodeName)
{
	unsigned int nodeNumber = nameNumber.size() + 1;
	
	nameNumber[nodeName] = nodeNumber;
	return nodeNumber;
};

//! Returns node data for a node.
Data * AllNodeData::getNodeData(const unsigned int & nodeNo) const
{
	map<unsigned int, Data *>::const_iterator nd = allNodeData.find(nodeNo);

	if(nd == allNodeData.end())
	{
		string mess = "Attempt to get node data but no such node data exists!";
		exitErr(mess);
	};

	return nd->second;
};

//! Returns node data number, but returns 0 if does not exist, allows nodes with no data for sim'ing data
unsigned int AllNodeData::getNodeDataNumberInit(const string & nodeName) const
{
	map<string, unsigned int>::const_iterator nn = nameNumber.find(nodeName);
	
	if(nn == nameNumber.end())
	{
		return 0; //does not exist
	};

	return nn->second;
};

//! Returns node data number (ID) given the name of a node.
unsigned int AllNodeData::getNodeDataNumber(const string & nodeName) const
{
	map<string, unsigned int>::const_iterator nn = nameNumber.find(nodeName);
	
	if(nn == nameNumber.end())
	{		

		string mess = "Attempt to find node data number of node data \"" + nodeName + "\", but no such data exists!\n"
			+"Have you set too many ID columns or forgot to insert an ID column?\n"
			+"Check that you have set the correct input format for the network.\n";

		exitErr(mess);
	};

	return nn->second;
};

//! Returns the node name given the node data number (ID).
string AllNodeData::getNodeDataName(const unsigned int & nodeDataNo) const
{
	for(map<string, unsigned int>::const_iterator nn = nameNumber.begin(); nn != nameNumber.end(); ++nn)
	{
		if(nn->second == nodeDataNo) return nn->first;

	};

	out("Attempt to find node data name of node data "); out(nodeDataNo); out(", but no such data exists!\n");
	exitErr("");

	return "";
};

//! Returns the node name given the node data number (ID).
string AllNodeData::getNodeDataNameError(const unsigned int & nodeDataNo) const
{
	for(map<string, unsigned int>::const_iterator nn = nameNumber.begin(); nn != nameNumber.end(); ++nn)
	{
		if(nn->second == nodeDataNo) return nn->first;
	};

	return "";
};

//! Returns the node data given the node name.
Data * AllNodeData::getNodeData(const string & nodeName) const
{
	map<string, unsigned int>::const_iterator nn = nameNumber.find(nodeName);
	
	if(nn == nameNumber.end())
	{
		string mess = "Attempt to get node data \"" + nodeName + "\", but no such data exists!";
		exitErr(mess);
	};
	
	map<unsigned int, Data *>::const_iterator nd = allNodeData.find(nn->second);

	if(nd == allNodeData.end())
	{
		string mess = "Attempt to get node data \"" + nodeName + "\", but no such data exists in node list!";
		exitErr(mess);
	};

	return nd->second;
};

//! Check data before doing stuff, may not catch all errors.
void AllNodeData::checkData()
{

	if(allNodeData.size() >= 2 && !dataChecked)
	{	
		map<unsigned int, Data *>::const_iterator nd = allNodeData.begin();
		unsigned int prevNum = nd->second->getAmountOfData();
		unsigned int num;
		string prevNodeName = nd->second->name;
		++nd;

		for( ; nd != allNodeData.end(); ++nd)
		{
			num = nd->second->getAmountOfData();

			if(num != prevNum)
			{
				outErr("Variable "); outErr(prevNodeName); outErr(" has "); outErr(prevNum); outErr(" data entries, but ");
				outErr("variable "); outErr(nd->second->name); outErr(" has "); outErr(num); outErr("!");
				exitErr("");
			};

			prevNum = num;
			prevNodeName = nd->second->name;
		};
	};

	dataChecked = true;
};

//! Update data names to remove ":" and other things maybe later.
void AllNodeData::updateDataNames()
{
	bool updated = false;
	bool updated1 = false;
	size_t posColon;
	string name = "";
	map<string, unsigned int>::iterator nn; // nameNumber

	for(map<unsigned int, Data *>::iterator nd = allNodeData.begin(); nd != allNodeData.end(); ++nd)
	{	
		updated1 = false;
		do{
			name = nd->second->name;
			posColon = nd->second->name.find_first_of(':');
			if(posColon < nd->second->name.length() && posColon != string::npos)
			{
				updated = true; updated1 = true;
				nd->second->name = nd->second->name.substr(0, posColon) + "_";
				if((posColon + 1) < name.length()) nd->second->name = nd->second->name + name.substr(posColon+1);	
			};
		}while(posColon < nd->second->name.length() && posColon != string::npos);

		if(updated1)
		{
			nn = nameNumber.find(nd->second->name);
			if(nn != nameNumber.end())
			{
				out("\nWARNING: variables with colons, \":\", are not allowed and have been replaced with underscores, \"_\".\n\n");
				string mess = "Attempt to add variable \"" + nd->second->name + "\" twice to the data set!";
				exitErr(mess);
			};
			nameNumber[nd->second->name] = nd->first; //add reference
		};
	};
	
	if(updated)
	{
		out("\nWARNING: variables with colons, \":\", are not allowed and have been replaced with underscores, \"_\".\n\n");
	};	
};

//! Returns number of indivs.
unsigned int AllNodeData::getAmountOfData() const
{
	if(allNodeData.size() > 0) return allNodeData.begin()->second->getAmountOfData();
	else return 0;
};

//! Adds a group of data set that belong to the same vairable stored as factors.
void AllNodeData::addFactorDataGroup(string headFactor, list<string> otherFactors)
{
	unsigned int headFactorNumber = getNodeDataNumber(headFactor);
	list<unsigned int> otherFactorsNumbers;

	for(list<string>::const_iterator of = otherFactors.begin(); of != otherFactors.end(); ++of)
	{	
		otherFactorsNumbers.push_back(getNodeDataNumber(*of));
	};

	factorGroups[headFactorNumber] = otherFactorsNumbers;
};

//! Returns the cts node IDs of the factors for a head factor node.
list<unsigned int> AllNodeData::getFactorDataGroup(unsigned int headNodeID)
{
	map<unsigned int, list<unsigned int> >::const_iterator fg = factorGroups.find(headNodeID);

	if(fg != factorGroups.end()) return fg->second;
	
	return list<unsigned int>();
};

//! Checks IDs in different files are the same.
void AllNodeData::checkIDs(string & filename, list<string> & dataIDsToCheck)
{
	//setup initial list of data IDs to check sebsequent data ids against
	if(dataIDs.empty())
	{
		dataIDs = dataIDsToCheck;
		dataFilename1 = filename;
		return;
	};

	if(dataIDs.size() != dataIDsToCheck.size()) exitErr("The amount of ID data in file "+dataFilename1+" is different to that in file "+filename+"!");

	list<string>::const_iterator d2 = dataIDsToCheck.begin();
	for(list<string>::const_iterator d = dataIDs.begin(); d != dataIDs.end(); ++d, ++d2)
	{		
		if(*d2 != *d) exitErr("The data IDs in file "+dataFilename1+" are different to those in file "+filename+"!");
	};

};
