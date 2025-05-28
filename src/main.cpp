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


/*! \file main.cpp
    \brief This file reads in the initial input files and options.
    
    This file also outputs usage instructions and program details.
*/

#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <string>

using namespace std; // initiates the "std" or "standard" namespace
 
#include "main.h"
#include "Analysis.h"

#ifdef USING_OPEN_MPI
#include "mpi.h"
int processRank = 0;
int processSize = 1;
#endif 

bool outputToScreen = true; 
bool freeClearMemory = false;

ofstream logFile;
stringstream stringLogFile;
string logFilename;

//! Outputs error messages on exit.
void exitErr(const string & message)
{
	outErr("\n"); outErr(message); outErr("\n\n");
	logFile.open(logFilename.c_str());
	logFile << stringLogFile.str();
	logFile.close();
    exit(1); //exit with 1 signal to say not normal termination
};

//! Output program title to screen.
void header()
{
	out("\nBayesNetty: Bayesian Network software, v1.1.1\n");
	out("--------------------------------------------------\n");
	out("Copyright 2015-present Richard Howey, GNU General Public License, v3\n");
	out("Institute of Genetic Medicine, Newcastle University\n\n");
};

//! Output program usage to screen.
void usage()
{
		header();
	 	
		out("Usage:\n\t ./bayesnetty [options] parameterfile.txt\n\n");
		
		out("General Options:\n");
		out("  -log file.log -- output log file (default: bayesnetty.log)\n");	
		out("  -so           -- suppress output to screen\n");
		out("  -seed         -- random number generator seed\n");
		out("\n");
		out("Input Data Options:\n");
		out("  -input-data                          -- do a task to input data\n");	
		out("  -input-data-name name                -- label the task with a name (default: Task-n)\n");
		out("  -input-data-file data.dat            -- the file containing the data for each network node\n");	
		out("  -input-data-include-file nodes.dat   -- a list of nodes/variables from the data file to be included\n");	
		out("  -input-data-exclude-file nodes.dat   -- a list of nodes/variables to be excluded from the network\n");	
		out("  -input-data-cts                      -- set the data file as containing continuous data\n");	
		out("  -input-data-cts-snp                  -- set the .bed data file as containing SNP data to be treated as continuous data\n");	
		out("  -input-data-cts-snp2                 -- set the data file as containing continuous SNP data \n");	
		out("  -input-data-cts-missing-value x      -- set the value of missing data for continuous data to x\n");		
		out("  -input-data-discrete                 -- set the data file as containing discrete data\n");	
		out("  -input-data-discrete-snp             -- set the .bed data file as containing SNP data to be treated as discrete data\n");	
		out("  -input-data-discrete-snp2            -- set the data file as containing discrete SNP data \n");
		out("  -input-data-discrete-missing-value x -- set the value of missing data for discrete data to x (default: NA)\n");
		out("  -input-data-factor                   -- set the data file as containing discrete data encoded using factor variables\n");	
		out("  -input-data-factor-snp               -- set the data file as containing SNP data to be treated as discrete factor data\n");	
		out("  -input-data-factor-snp2              -- set the data file as containing discrete factor SNP data \n");
		out("  -input-data-factor-missing-value x   -- set the value of missing data for discrete data to x (default: NA)\n");
		out("  -input-data-ids n                    -- the number of ID columns in each data file (default: 2)\n");
		out("  -input-data-csv                      -- set the data file as a comma separated file, .csv\n");
		out("\n");
		out("Input Network Options:\n");
		out("  -input-network                                         -- do a task to input a network\n");	
		out("  -input-network-name name                               -- label the task and network with a name (default: Task-n)\n");
		out("  -input-network-type t                                  -- the type of network, choose between bnlearn or deal (default: bnlearn)\n");
		out("  -input-network-file network.dat                        -- input the network in a format where the nodes and then the edges are listed\n");	
		out("  -input-network-file2 network2.dat                      -- input the network in this style of format: [a][b|a][c|a:b]\n");	
		out("  -input-network-igraph-file-prefix mygraph              -- input the network from igraph format files consisting of mygraph-nodes.dat and mygraph-edges.dat\n");	
		out("  -input-network-empty                                   -- set the network to one with no edges and one node for every data variable.\n");	
		out("  -input-network-whitelist-file whitelist.dat            -- a list of edges that must be included in any network\n");	
		out("  -input-network-blacklist-file blacklist.dat            -- a list of edges that must not be included in any network\n");	
		out("  -input-network-blacklist-edge-type dataName1 dataName2 -- edge types that may not be included in any network. The collection of nodes are given by the data input name\n");	
		out("  -input-network-imaginary-sample-size i                 -- for deal networks this sets the imaginary sample size (default: 10)\n");
		out("  -input-network-score score                             -- for a bnlearn network choose between loglike, AIC or BIC (default: BIC)\n");
		//out("  -input-network-score-fix                               -- when too few points for linear regression, choose between none, average or skip (default: none)\n\n");
		out("  -input-network-no-parents-node nodeX                   -- nodeX must not have any parents (except for white edges)\n");
		out("  -input-network-no-children-node nodeY                  -- nodeY must not have any children (except for white edges)\n");
		out("  -input-network-prob-edge node1 node2 prob              -- set the probability of edge direction of node1 to node2 as prob\n");
		out("  -input-network-prob-edge-type nodeType1 nodeType2 prob -- set the probabilities of edge direction of nodeType1 to nodeType2 as prob\n");
		out("\n");
		out("Calculate Posterior Options\n");
		out("  -calc-posterior                      -- do a task to calculate the posterior\n");	
		out("  -calc-posterior-name name            -- label the task with a name	Task-n\n");
		out("  -calc-posterior-network-name network -- the name of the network to calculate the posterior of (default: previous network)\n");
		out("\n");
		out("Calculate Network Score Options\n");
		out("  -calc-network-score                               -- do a task to calculate the score\n");	
		out("  -calc-network-score-name name                     -- label the task with a name (default: Task-n)\n");	
		out("  -calc-network-score-network-name network          -- the name of the network to calculate the score of (default: previous network)\n");
		out("  -calc-network-score-file file                     -- write the score to this file\n");
		out("  -calc-network-score-all-scores network-scores.dat -- calculate the scores of every possible network and record the results in network-scores.dat\n");
		out("\n");
		out("Search Models Options\n");
		out("  -search-models                      -- do a task to search network models\n");	
		out("  -search-models-name name            -- label the task with a name (default: Task-n)\n");
		out("  -search-models-network-name network -- the name of the network to start the search from (default: previous network)\n");
		out("  -search-models-file search.dat      -- record the network models and scores in the search path to file search.dat\n");	
		out("  -search-models-random-restarts n    -- do another n searches starting from a random network (default: 0)\n");
		out("  -search-models-jitter-restarts m    -- after the initial search and every random restart search do another m searches jittered from the recently found network (default: 0)\n");
		out("\n");
		out("Average Networks Options\n");
		out("  -average-networks                            -- do a task to calcualte an average network\n");	
		out("  -average-networks-name name                  -- label the task with a name (default: Task-n)\n");
		out("  -average-networks-network-name network       -- the name of the network to calculate the average network from (default: previous network)\n");
		out("  -average-networks-file average-network.dat   -- the name of the file to record the average network in\n");	
		out("  -average-networks-igraph-file-prefix mygraph -- output igraph format files consisting of mygraph-nodes.dat, mygraph-edges.dat and R code mygraph-plot.R\n");	
		out("  -average-networks-threshold thres            -- the strength threshold used to include an edge in the drawn average network (default: estimated)\n");
		out("  -average-networks-bootstraps k               -- the number of bootstraps used to calculate the average network (default: 100)\n");	
		out("  -average-networks-random-restarts n          -- for each network fit do another n searches starting from a random network (default: 0)\n");
		out("  -average-networks-jitter-restarts m          -- for each network fit after the initial search and every random restart search do another m searches jittered from the recently found network (default: 0)\n");
		out("  -average-networks-use-weight-method          -- use edge chi sq values to weight edge strengths\n");
		out("  -average-networks-use-score-method           -- use network score method instead of bootstrapping (slow)\n");
		out("  -average-networks-likelihood-file            -- output likelihoods of separate bootstrap networks to file (slow)\n");
		out("\n");
		out("Measurement Error Robustness Options\n");
		out("  -measurement-error-robustness                            -- do a task to analyse robustness of the fitted network when accounting for measurement error\n");
		out("  -measurement-error-robustness-name name                  -- label the task with a name (default: Task-n)\n");
		out("  -measurement-error-robustness-network-name network       -- the name of the network to use (default: previous network)\n");
		out("  -measurement-error-robustness-ave-file ave-net.dat       -- the name of the file to record the average network in\n");
		out("  -measurement-error-robustness-stdev e                    -- the standard deviation of the measurement error for all vairaibles (default: 0)\n");
		out("  -measurement-error-robustness-stdev-file stdevs.dat      -- the standard deviation of measurement errors for each variable (default: above)\n");		
		out("  -measurement-error-robustness-node-stdev nd e            -- the standard deviation of the measurement error for node/variable nd\n");
		out("  -measurement-error-robustness-file average-network.dat   -- the name of the file to record the average network in\n");
		out("  -measurement-error-robustness-igraph-file-prefix mygraph -- output igraph format files consisting of mygraph-nodes.dat, mygraph-edges.dat and R code mygraph-plot.R\n");
		out("  -measurement-error-robustness-threshold thres            -- the strength threshold used to include an edge in the drawn average network (default: estimated)\n");
		out("  -measurement-error-robustness-iterations k               -- the number of bootstraps used to calculate the average network (default: 100)\n");
		out("  -measurement-error-robustness-random-restarts n          -- for each network fit do another n searches starting from a random network (default: 0)\n");
		out("  -measurement-error-robustness-jitter-restarts m          -- for each network fit after the initial search and every random restart search do another m searches jittered from the recently found network (default: 0)\n");
		out("\n");
		out("Impute Network Data Options\n");
		out("  -impute-network-data                         -- do a task to impute network data\n");	
		out("  -impute-network-data-name name               -- label the task with a name (default: Task-n)\n");
		out("  -impute-network-data-network-name network    -- the name of the network to impute data for (default: previous network)\n");
		//out("  -impute-network-data-max-missing") getLatestImputeNetworkDataTask()->setMaxMissing(getUIntOptionValue());
		out("  -impute-network-data-min-non-missing-edges x -- the percentage (0 to 100) of non-missing edges required to impute data for an individual (default: 50)\n");
		//out("  -impute-network-data-use-prev-network") getLatestImputeNetworkDataTask()->setUsePrevNetwork();
		out("  -impute-network-subsample-percent n          -- percentage of data to use when taking a subsample for training data (default 90). Set to 0 to use all data and fit only once. Advised for v. large datasets.\n");
		out("  -impute-network-data-complete-training       -- use complete data for training data (default if complete data is at least 40)\n");
	    out("  -impute-network-data-random-training         -- use randomly sampled values for missing data for training data\n"); 
		out("  -impute-network-data-random-restarts n       -- for each bootstrap network fit do another n searches starting from a random network (default: 0)\n");
		out("  -impute-network-data-jitter-restarts m       -- for each bootstrap network fit after the initial search and every random restart search do another m searches jittered from the recently found network (default: 0)\n");
		out("  -impute-network-data-start-indiv a           -- start imputing data from individual a (default: individual 1)\n");
		out("  -impute-network-data-end-indiv b             -- end imputing data at individual b (default: last individual)\n");
		out("  -impute-network-data-job i t                 -- only impute individuals for subset i from a total of t subsets\n");
		out("\n");
		out("Estimate Recall and Precision Before and After Data Imputation Options\n");
		out("  -impute-estimate-recall-precision                         -- do a task to estimate recall and precision before and after data imputation\n");
		out("  -impute-estimate-recall-precision name                    -- label the task with a name (default: Task-n)\n");
		out("  -impute-estimate-recall-precision-network-name network    -- set the name of the initial network when estimating recall and precision (default: previous network)\n");
		out("  -impute-estimate-recall-precision-random-restarts n       -- for each network fit do another n searches starting from a random network (default: 0)\n");
		out("  -impute-estimate-recall-precision-jitter-restarts m       -- for each network fit after the initial search and every random restart search do another m searches jittered from the recently found network (default: 0)\n");
		out("  -impute-estimate-recall-precision-skip-imputation         -- do not estimate recall and precision for imputed data\n");
		out("  -impute-estimate-recall-precision-iterations i            -- estimate the recall and precision i times and take the average (default: 1)\n");
		out("  -impute-estimate-recall-precision                         -- do a task to estimate recall and precision before and after data imputation\n");
		out("  -impute-estimate-recall-precision name                    -- label the task with a name (default: Task-n)\n");
		out("\n");
		out("Calculate Recall and Precision\n");
		out("  -calculate-recall-precision                              -- do a task to calculate the recall and precision of a network\n"); 
		out("  -calculate-recall-precision-name name                    -- label the task with a name (default: Task-n)\n");
		out("  -calculate-recall-precision-network-name network1        -- the name of the network to calculate the recall and precision for (default: previous network)\n");
		out("  -calculate-recall-precision-true-network-name network2   -- the name of the true network to calculate the recall and precision against\n");
		out("  -calculate-recall-precision-file results.dat             -- file to write recall and precision results to\n");
		out("Simulate Network Data Options\n");
		out("  -simulate-network-data                               -- do a task to simulate network data for a given network\n");	
		out("  -simulate-network-data-name name                     -- label the task with a name (default: Task-n)\n");
		out("  -simulate-network-data-network-name network          -- simulate data for this network (default: previous network)\n");
		out("  -simulate-network-data-no-sims n                     -- simulate n replicates of network data (default: 100)\n");
		out("  -simulate-network-data-parameter-file parameters.txt -- network with parameters in bnlearn posteriors file format\n");	
		out("  -simulate-network-data-whitelist-file whitelist.dat  -- a list of edges that must be included in any network\n");	
		out("  -simulate-network-data-blacklist-file blacklist.dat  -- a list of edges that must not be included in any network\n");	
		//out("  -simulate-network-data-imaginary-sample-size i       -- for deal networks this sets the imaginary sample size (default: 10)\n");
		out("  -simulate-network-data-score score                   -- for a bnlearn network choose between loglike, AIC or BIC (default: BIC)\n");
		out("  -simulate-network-data-score-fix f                   -- when too few points for logistic regression, choose between none, average or skip (default: none)\n");
		out("\n");
		out("Markov Blanket Options\n");
		out("  -markov-blanket\n");
		out("  -markov-blanket-name name            -- label the task with a name (default: Task-n)\n");
		out("  -markov-blanket-network-name network -- the name of the network for calculating the Markov blanket (default: previous network)\n");
		out("  -markov-blanket-node-name node       -- calculate the Markov blanket for the node with this name\n");
		out("\n");
		out("Output Network Options\n");
		out("  -output-network                                             -- do a task to output a network to file\n");	
		out("  -output-network-name name                                   -- label the task with a name (default: Task-n)\n");
		out("  -output-network-network-name network                        -- output this network (default: previous network)\n");
		out("  -output-network-file network.dat                            -- output the network in a format where the nodes and then the edges are listed\n");	
		out("  -output-network-file2 network2.dat                          -- output the network in this style of format: [a][b|a][c|a:b]\n");	
		out("  -output-network-equivalent-networks-file equiv-networks.dat -- output a list of equivalent networks to file equiv-networks.dat\n");	
		out("  -output-network-igraph-file-prefix mygraph                  -- output igraph format files consisting of mygraph-nodes.dat, mygraph-edges.dat and R code mygraph-plot.R\n");	
		out("  -output-network-node-data-file-prefix mydata                -- output network discrete data to mydata-discrete.dat and continuous data to mydata-cts.dat\n");	
		out("  -output-network-node-data-bed-file                          -- output SNP data to files mydata.bed/.bim/.fam and not to files mydata-discrete.dat and mydata-cts.dat\n");	
		out("  -output-network-node-data-start-indiv a                     -- start output of data from individual a (default: individual 1)\n");
		out("  -output-network-node-data-end-indiv b                       -- stop output of data at individual b (default: last individual)\n");
		out("  -output-network-node-data-job i t                           -- only output individuals for subset i from a total of t subsets\n");
		out("\n");
		out("Output Priors Options\n");
		out("  -output-priors                      -- do a task to output the priors of a network to file\n");	
		out("  -output-priors-name name            -- label the task with a name (default: Task-n)\n");
		out("  -output-priors-network-name network -- output priors for this network (default: previous network)\n");
		out("  -output-priors-file priors.dat      -- output the priors to file priors.dat (default: priors.dat)\n");
		out("\n");
		out("Output Posteriors Options\n");
		out("  -output-posteriors                      -- do a task to output the posteriors of a network to file\n");	
		out("  -output-posteriors-name name            -- label the task with a name (default: Task-n)\n");
		out("  -output-posteriors-network-name network -- output posteriors for this network (default: previous network)\n");
		out("  -output-posteriors-file posts.dat       -- output the posteriors to file posts.dat (default: posteriors.dat)\n");

};

//! The start of the program.
int main(int argc, char * argv[])
{
	outputToScreen = true;

#ifdef USING_OPEN_MPI
	//set up different processes
	MPI_Init(NULL, NULL);

	//get rank (processer number 0 to size-1) and total number of processes
	MPI_Comm_rank(MPI_COMM_WORLD, &processRank);	
	MPI_Comm_size(MPI_COMM_WORLD, &processSize);

	if(processRank != 0) outputToScreen = false;
#endif

	time_t start,end;
	double dif;
	time(&start);
	
	int argcount = 1;
	string option = "";
	logFilename = "bayesnetty.log";
	
	//create analysis option to add Tasks to from parameter file
	Analysis anAnalysis;
	anAnalysis.setCommandArgs(argcount, argc, argv);
	

	//set given options
	while(argcount < argc)
	{
		option = argv[argcount];
	
		if(option.substr(0,1) == "-")
		{
			if(!anAnalysis.processOption(option)) //process option
			{
				logFile.open(logFilename.c_str());
			
				header();
    			string message = "Unrecognised command line option: " + option + "!";

				exitErr(message);					
			}; 
		}
		else
		{
			anAnalysis.processParameterFile(option); //assume option is a parameter file 
		};	
		
		anAnalysis.updateArgcount(argcount); //add one and keep in snyc with copy in anAnalysis
	};

	if(argc==1)
	{
		usage();
		exit(0);
	};
	
	//output options to screen
	header();

	unsigned int seed = anAnalysis.getRandomSeed();
	if(seed == 0) seed = (unsigned int)time(NULL);

	
#ifdef USING_OPEN_MPI
	 //broadcast seed, ensure all processes have same random seed                                                                                                                                                            
	MPI_Bcast(&seed, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

	out("Number of processes: "); out(processSize); out("\n");
#endif

	anAnalysis.setRandomSeed(seed);
	srand(seed);

	out("Random seed: "); out(seed); out("\n");
	
	//set up tasks
	//do tasks set by the user
	anAnalysis.doTasks();
	
	//output analysis options and summary of results	
	//output run time
	time(&end);
	dif = difftime(end, start);
	out("\nRun time: "); out(getTime(dif)); out("\n\n");

	//write log to file from the string in memory
#ifdef USING_OPEN_MPI
	if(processRank==0)
	{
#endif
		logFile.open(logFilename.c_str());
		logFile << stringLogFile.str();			
		logFile.close();
#ifdef USING_OPEN_MPI
	};

	MPI_Finalize();
#endif

};

