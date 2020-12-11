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

/*! \file Search.cpp
    \brief This file contains the methods for searching the networks.
    
*/

#include <iostream>
#include <sstream>

using namespace std; // initiates the "std" or "standard" namespace

#include "main.h"
#include "Search.h"

#ifdef USING_OPEN_MPI
#include "mpi.h"
#endif 

//! Output the search for a greedy search.
double GreedySearch::doSearch(Network * network)
{
#ifdef USING_OPEN_MPI
	return doSearchParallel(network);
#else
	return doSearchSerial(network);
#endif
};

//! Output the search for a greedy search in serial (not parallel).
double GreedySearch::doSearchSerial(Network * network)
{
	double networkScore;

	//setup initial priors for initial network
	hasLoop = network->hasLoop();

	if(!hasLoop)
	{
		network->calculatePriorsAndPosteriorsBNL();//only for bnlearn

		networkScore = network->cacheAllNodes(); 
	}
	else
	{
		return 0;
	};

	bool updated = true;

	network->evaluateAllEdgeChangesDiffs();

	while(updated)
	{			
		if(outputSearch && networkScore*0 == 0)
		{
			outSearch << networkScore << " " << network->getNetworkString(0) << "\n";
		};

		//consider all changes (only for deal networks)
		network->evaluateAllEdgeChanges();
		
		updated = network->updateNetwork(networkScore);
	
	};

	return networkScore;
};

#ifdef USING_OPEN_MPI
//! Output the search for a greedy search in parallel using Open MPI.
double GreedySearch::doSearchParallel(Network * network)
{
	
	double networkScore;

	//setup initial priors for initial network
	hasLoop = network->hasLoop();

	if(!hasLoop)
	{
		network->calculatePriorsAndPosteriorsBNL();//only for bnlearn

		networkScore = network->cacheAllNodes(); 
	}
	else
	{
		return 0;
	};

	bool updated = true;

	network->evaluateAllEdgeChangesDiffsParallel();

	while(updated)
	{	
			
		if(outputSearch && networkScore*0 == 0)
		{
			outSearch << networkScore << " " << network->getNetworkString(0) << "\n";
		};

		//consider all changes (only for deal networks)
		network->evaluateAllEdgeChanges();
		
		updated = network->updateNetworkParallel(networkScore);
	
	};


	return networkScore;

};
#endif
