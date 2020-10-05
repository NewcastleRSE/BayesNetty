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


/*! \file main.h
    \brief This file defines a global variable to suppress output to screen. 
    
*/

// Comment out if not using Open MPI for parallel processing
//#ifndef USING_OPEN_MPI
//#define USING_OPEN_MPI
//#endif //OPEN_MPI

#ifndef __MAIN
#define __MAIN

#include <iostream>
#include <ostream>
#include <iomanip>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <string>

using namespace std;

extern bool outputToScreen;
extern bool freeClearMemory; //free up memory after a clear, nec for large networks when imputing or taking average networks etc, when network fitted many times

#ifdef USING_OPEN_MPI
extern int processRank; //for parallel MPI programming, rank 0 is master process
extern int processSize; //number of processes 
#endif

extern string logFilename;
extern ofstream logFile; 
extern stringstream stringLogFile; //avoid logfile going "missing", bus errors for long jobs with virtual memory pages disappearing in linux

template<typename T>

//! Outputs message to screen and log file.
void out(const T & text)
{	
	if(outputToScreen) cout << text;
	stringLogFile << text;
};

template<typename T>

//! Outputs error message to screen and log file.
void outErr(const T & text)
{
	cerr << text;
	stringLogFile << text;
};

//! Exit program due to error.
void exitErr(const string & message);
void header();


#endif
