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


/*! \file Search.h
    \brief This file searches the possible networks.
*/


#ifndef __SEARCH
#define __SEARCH

#include <list>
#include <string>
#include <map>
#include <set>
#include "main.h"
#include "Network.h"

using namespace std; // initiates the "std" or "standard" namespace

//! Search networks for best fitting model.
class SearchNetworks
{
private:

protected:

	bool hasLoop;
	ofstream outSearch;
	bool outputSearch;
	
public:

	SearchNetworks() : hasLoop(false), outputSearch(false) {};

	~SearchNetworks()
	{
	
	};

	//general search methods
	virtual double doSearch(Network * network) {return 0;};
	
	//search network models methods
	virtual void setSearchOutput(string & filename)
	 {		
		outSearch.open(filename.c_str());
		outSearch << "logLike Model\n";
		outputSearch = true;
	};
	void closeFileSearch() {outSearch.close();};
	void outputSearchText(const string & txt) {outSearch << txt;};
	void outputSearchNum(const double & num) {outSearch << num;};
	void outputSearchNum(const unsigned int& num) { outSearch << num; };
};

//! Search networks for best fitting model.
class GreedySearch : public SearchNetworks
{
private:

	
public:

	GreedySearch() : SearchNetworks() {};

	~GreedySearch() {};

	//general methods
	double doSearch(Network * network);
	double doSearchSerial(Network * network);
#ifdef USING_OPEN_MPI
	double doSearchParallel(Network * network);
#endif
};


#endif
