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


/*! \file Data.h
    \brief This file contains classes for manipulating SNP data. 
*/

#ifndef __DATA
#define __DATA

#include <list>
#include <map>
#include <string>
#include <iostream>
#include "main.h"

using namespace std;

//! Stores SNP data for all subjects for a given traits.
struct Data
{
	string name; //column name of data
	unsigned int fileNo; //where data came from, 0 = not from file as it is from sim'd data
	unsigned int dataType; // 0 = IDs, 1 = discrete, 2 = cts, 3 = no data
	bool isSNPData;	
	list<bool> missingValues; //true if missing
	list<bool> wasImputed; //true if recently imputed

	virtual unsigned int getAmountOfData() {return 0;};
	virtual void clearData() {};
	virtual bool getIsFactorNode() const {return false;};
	virtual bool getIsFactorChildNode() const {return false;};
	virtual void setMissingData(const list<bool> & md) {};
	void setupWasImputed();
	void updateImputedDataAsNonMissing();

	Data(const string & nm, const unsigned int & fno, const unsigned int & type) : name(nm), fileNo(fno), dataType(type), isSNPData(false),  wasImputed(), missingValues() {};

	virtual ~Data()
	{

	};
	
};

//! Keeps track of whether there is any missing for any node in network.
struct NetworkMissingData
{
	list<bool> missing; //true if missing data for any node in network

	NetworkMissingData() : missing() {};

	~NetworkMissingData() {};

	unsigned int getNoMissing();
};

//! Stores SNP data for all subjects for a given SNP.
struct CtsData : public Data
{
	
	list<double> values;      //subject no, quantitive trait
	
	bool isFactorNode;
	bool isFactorChildNode;
	string nodeName;
	string factorName;
	map<unsigned int, string> levels; //level no., name of factor, 
	map<unsigned int, double> completeValues; // list of complete values used for sampling from when imputing

	CtsData(const string & nm, unsigned int & fno) : Data(nm, fno, 2), values(), isFactorNode(false), isFactorChildNode(false), nodeName(""), factorName("")  {};

	~CtsData()
	{
		 clearData();
	};
	
	unsigned int getAmountOfData() {return values.size();};
	void addData(const double & num, const bool & missing);
	void clearData()
	{
		if(freeClearMemory)
		{
			list<double>().swap(values); list<bool>().swap(missingValues); list<bool>().swap(wasImputed);
		}
		else
		{
			values.clear(); missingValues.clear(); wasImputed.clear();
		};
	};
	double getMean();
	double getStDev();
	void setMissingData(const list<bool> & md);
	string getFactorLevelName(const unsigned int & level);
	void setLevelName(const unsigned int & lv, const string & lvName) {levels[lv] = lvName;};
	double getRandomCompleteValue();
};

//! Stores SNP data for all subjects for a given SNP.
struct NoData : public Data
{
	
	NoData(const string & nm) : Data(nm, 0, 3) {};

	~NoData()
	{


	};
	
};

//! Stores Discrete data for all subjects, i.e. data falls into categories.
struct DiscreteData : public Data
{
	
	list<unsigned int> values; //ordered values, missing value given by 0
	
	map<string, unsigned int> levels; //name of factor, level no.	

	map<double, unsigned int> levelProbsCum; //cumlative probs of getting levels 1, 2, 3... used for random selection when imputing

	DiscreteData(const string & nm, unsigned int & fno) : Data(nm, fno, 1), values() {};

	~DiscreteData()
	{
		 clearData();		
	};
	
	unsigned int getAmountOfData() {return values.size();};
	void clearData()
	{
		if (freeClearMemory)
		{
			list<unsigned int>().swap(values); list<bool>().swap(missingValues); list<bool>().swap(wasImputed);
		}
		else
		{
			values.clear(); missingValues.clear(); wasImputed.clear();
		};
	};
	void addData(const string & value, bool & missing);
	void addSimData(const unsigned int & level, const bool & missing);
	unsigned int addLevel(const string & value);
	string getLevelName(const unsigned int & level);
	void setMissingData(const list<bool> & md);
	unsigned int getRandomLevel();
};

//! Covariate data for all subjects, used for node parent data.
struct CovariateData
{
	CtsData * nodeValues;
	list<CtsData *> covariateDataAllSubjects; //list of parent data for cts nodes (if any)
	list<bool> notInGroupOrMissingData;       //false if not in discrete parent group or missing data anywhere in parent data

	CovariateData() : nodeValues() {};
	
	void clear()
	{ 
		if(freeClearMemory)
		{
			list<CtsData*>().swap(covariateDataAllSubjects);
			list<bool>().swap(notInGroupOrMissingData);
		}
		else
		{
			covariateDataAllSubjects.clear(); notInGroupOrMissingData.clear();
		};
	};

	~CovariateData()
	{

	};

};

//! Class to act as central storage and linking place for all data.
class AllNodeData
{
private:

	map<unsigned int, Data *> allNodeData; //potential node data for any network
	map<string, unsigned int> nameNumber;

	list<string> dataIDs; //list of data IDs (multiple IDs one after another), used for check when readin data from file
	list<string> nameOfIDs; //name of the ID columns
	string dataFilename1; //name of first data file used to set up the dataIDs
	map<unsigned int, list<unsigned int> > factorGroups; //data for a factor stored in multiple CtsData variables, first one is the head factor 

	bool dataChecked;

public:
#


	AllNodeData() : dataChecked(false)
	{ 
		
	};

	//! Delete Task things
	~AllNodeData()
	{
		for(map<unsigned int, Data *>::iterator nd = allNodeData.begin(); nd != allNodeData.end(); ++nd) delete nd->second;
		allNodeData.clear();
	};


	void addNodeData(Data * data);
	void replaceNodeData(Data * data);
	void deleteNodeData(const unsigned int & nodeNumber);
	bool isNodeDataAdded(string & dataName);
	unsigned int getNewNodeNo(const string & nodeName);
	Data * getNodeData(const string & nodeDataName) const;
	Data * getNodeData(const unsigned int & nodeDataNo) const;
	map<string, unsigned int> getNameNumber() const {return nameNumber;};
	unsigned int getNodeDataNumber(const string & nodeDataName) const;
	unsigned int getNodeDataNumberInit(const string & nodeName) const;
	string getNodeDataName(const unsigned int & nodeDataNo) const;
	string getNodeDataNameError(const unsigned int & nodeDataNo) const;
	void checkData();
	unsigned int getAmountOfData() const;
	unsigned int getNoNodeData() const {return nameNumber.size();};
	void checkIDs(string & filename, list<string> & dataIDsToCheck);
	void setNameOfIDs(list<string> & nids) {nameOfIDs = nids;};
	void setIDs(list<string> & ids) {dataIDs = ids;};
	list<string> getNameOfIDs() {return nameOfIDs;};
	list<string> getDataIDs() {return dataIDs;};
	void addFactorDataGroup(string headFactor, list<string> otherFactors);
	list<unsigned int> getFactorDataGroup(unsigned int headNodeID);
};

#endif
