.. _impute-data:

Impute Data
===========

When there is missing data, the standard approach is to remove every individual with missing data before performing any Bayesian network analysis, and this is the default behaviour.

This can be wasteful and undesirable when there are many individuals with missing data, perhaps with only one variable missing, making imputation a natural choice.

BayesNetty includes a new imputation method designed to increase the power to detect causal relationships whilst accounting for model uncertainty.
This method uses a version of nearest neighbour imputation, whereby missing data from one individual is replaced with data from another individual, the nearest neighbour.

An important feature of this approach is that it can be used with both discrete and continuous data.
 
For each individual with missing data, subsets of variables that can be used to find the nearest neighbour are chosen by bootstrapping the complete data to estimate a Bayesian network. 


.. _impute-data-options:

Options
-------

The options are as follows:


.. list-table:: 
    :header-rows: 1

    * - Option
      - Description
      - Default

    * - -impute-network-data
      - do a task to impute network data
      -

    * - -impute-network-data-name name
      - label the task with a name
      - Task-n

    * - -impute-network-data-network-name network
      - the name of the network to impute data for
      - previous network (or the default model given by a node for each data variable and no edges if there is no previous network) 

    * - -impute-network-data-min-non-missing-edges x
      - the percentage (0 to 100) of non-missing edges required to impute data for an individual
      - 0

    * - -impute-network-subsample-percent n
      - percentage of data to use when taking a subsample for training data. Set to 0 to use all data and fit only once - advised for v. large datasets.
      - 90

    * - -impute-network-data-complete-training
      - use complete data for training data (default if complete data is at least 40 percent)
      - 

    * - -impute-network-data-random-training
      - use randomly sampled values for missing data for training data (default if complete data is less than 40 percent)
      - 

    * - -impute-network-data-random-restarts n
      - for each subsample network fit do another n searches starting from a random network
      - 0

    * - -impute-network-data-jitter-restarts m
      - for each subsample network fit after the initial search and every random restart search do another m searches jittered from the recently found network
      - 0

    * - -impute-network-data-start-indiv a
      - start imputing data from individual a
      - 1

    * - -impute-network-data-end-indiv b
      - end imputing data at individual b
      - last individual

    * - -impute-network-data-job i t
      - only impute individuals for subset i from a total of t subsets
      -



The only option that is necessary to impute data is the `-impute-network-data` option. 

Once data has been imputed its missing data is filled in with imputed values and any subsequent analyses in BayesNetty will use this imputed data. 

If an individual has too much missing data then it may not be beneficial to impute the data for this individual,
as the imputed data would be too poor to add value to any analysis. The `-impute-network-min-non-missing-edges` allows the user to change the required amount of edges between non-missing variables to impute data for an individual.
Around **50** percent has been shown to be a suitable value to impute most individuals whilst effectively discarding individuals with too much missing data,
although it may depend on the structure of any fitted networks. If this value is set to 0 then all individuals will have their data imputed, and a value of 100 will result in no data being imputed.
If a dataset has a block of data with non-missing data for only a few variables then it is best to simply remove these individuals before using BayesNetty.

The `-impute-network-data-random-restarts` and `-impute-network-data-jitter-restarts` options can be increased to improve the network search at each step of the algorithm and may potentially increase the quality of the imputed data at the expense of a longer running time.

The `-impute-network-data-start-indiv`, `-impute-network-data-end-indiv` and `-impute-network-data-job` options may be used to only impute a range of individuals.
These options may be useful for large networks to impute data in parallel and then combine later if the data is output to file
(see :ref:`output-network` to output data). See :ref:`impute-parallel-example` for an example. 

The option `-impute-network-data-job` can also be used to only impute data for some individuals and makes it easier to split the imputation into a number of jobs.

The imputation process can be slow and computationally intensive, as it requires identifying a best fit network for each individual. If the dataset that you wish to impute is very large then this can be prohibitive - however,
the imputation process can be sped up by using the option `-impute-network-subsample-percent 0`. This fits only one best fit network for all individuals using all data. This has demonstrated performance approaching that of the standard approach.


.. _impute-example:

Example
-------

The following is an example parameter file to impute network data and search for the best network both before and after imputation.


.. code-block:: none

    #input continuous data
    -input-data
    -input-data-file impute-example-cts.dat
    -input-data-cts

    #input SNP data as discrete data
    -input-data
    -input-data-file impute-example.bed
    -input-data-discrete-snp

    #search network models with the original data
    -search-models

    #impute the missing data
    -impute-network-data

    #search network models with the imputed data
    -search-models


This parameter file, `paras-impute.txt`, and example data for imputation can be found in `impute-example.zip <https://github.com/NewcastleRSE/BayesNetty/raw/refs/heads/main/docs/resources/impute-example.zip>`_ and can be used as follows:


.. code-block:: none

    ./bayesnetty paras-impute.txt


Which should produce output that looks like something as follows:

.. code-block:: none

    BayesNetty: Bayesian Network software, v1.00
    --------------------------------------------------
    Copyright 2015-present Richard Howey, GNU General Public License, v3
    Institute of Genetic Medicine, Newcastle University

    Random seed: 1545221384
    --------------------------------------------------
    Task name: Task-1
    Loading data
    Continuous data file: impute-example-cts.dat
    Number of ID columns: 2
    Including (all) 5 variables in analysis
    Each variable has 1000 data entries
    Missing value: not set
    --------------------------------------------------
    --------------------------------------------------
    Task name: Task-2
    Loading data
    SNP binary data file: impute-example.bed
    SNP data treated as discrete data
    Total number of SNPs: 2
    Total number of subjects: 1000
    Number of ID columns: 2
    Including (all) 2 variables in analysis
    Each variable has 1000 data entries
    --------------------------------------------------
    --------------------------------------------------
    Task name: Task-3
    Searching network models
    --------------------------------------------------
    Loading defaultNetwork network
    Network type: bnlearn
    Network score type: BIC
    Total number of nodes: 7 (Discrete: 2 | Factor: 0 | Continuous: 5)
    Total number of edges: 0
    Network Structure: [bio1][bio2][bio3][trait1][trait2][rs1][rs2]
    Total data at each node: 213
    Missing data at each node: 787
    --------------------------------------------------
    Network: defaultNetwork
    Search: Greedy
    Random restarts: 0
    Random jitter restarts: 0
    Network Structure: [rs1][rs2][trait2|rs2][bio2|trait2][trait1|bio2][bio1|trait1][bio3|bio1:bio2]
    Network score type: BIC
    Network score = -1970.2
    --------------------------------------------------
    --------------------------------------------------
    Task name: Task-4
    Imputing network data
    Network: defaultNetwork
    Network Structure: [rs1][rs2][trait2|rs2][bio2|trait2][trait1|bio2][bio1|trait1][bio3|bio1:bio2]
    Number of individuals with missing data: 787
    Number of individuals imputed: 787
    Percentage of data imputed (when attempted): 98.4466
    Minimum percentage of non-missing edges (or singleton nodes) required to impute individual: 50
    Random restarts: 0
    Random jitter restarts: 0
    --------------------------------------------------
    --------------------------------------------------
    Task name: Task-5
    Searching network models
    Network: defaultNetwork
    Search: Greedy
    Random restarts: 0
    Random jitter restarts: 0
    Network Structure: [bio1][bio2][bio3|bio1:bio2][rs1][rs2][trait1|bio1:rs1][trait2|bio2:rs2]
    Network score type: BIC
    Network score = -9240.19
    --------------------------------------------------

    Run time: 34 seconds



The data is loaded, a search is performed and then the network data is imputed and another search is performed. The run time for performing imputation is longer than most other operations in BayesNetty.
This is because, every individual with missing data, we take a 90% sample (without replacement) of the individuals with complete data at the variables of interest. This sampled data set is used to find a best fit network.
This best fit network determines the variables that are used to choose the nearest neighbour for the individual with missing data,
and then the missing data is filled in from the nearest neighbour.

There are a lot of individuals with missing data in this example data resulting in the incorrect network being estimated initially but after the data is imputed the correct network is found.
That is, the network that the data was simulated from.


It may be possible that some individuals are not imputed as they have too much missing data, or sometimes only partially imputed if the data is not suitable for the imputation algorithm.

.. _impute-parallel-example:

Parallel Example
----------------


As imputing network data is a computationally intensive task, it makes sense to do it in parallel.
This can be done by running the parallel version of BayesNetty as described in :ref:`parallel`,
but a much quicker way is given here by running the non-parallel version of BayesNetty in parallel where each process imputes a subset of the individuals.
The data of the imputed individuals can then be output for each process (see :ref:`output-network`) and then combined into the final imputed data set.

A handy Unix script has been written to do this and is ran as follows:


.. code-block:: none

    ./runImputeParallel paras-impute-parallel.txt imputed-data 20


The first argument is a Bayesnetty parameter file to impute the data (example shown below).
The second argument is a file name (without extension) for the imputed data set to be outputted to. The last argument is the number of processes to run.  


.. code-block:: none

    #input continuous data
    -input-data
    -input-data-file impute-example-cts.dat
    -input-data-cts

    #input SNP data as discrete data
    -input-data
    -input-data-file impute-example.bed
    -input-data-discrete-snp

    #impute the missing data
    -impute-network-data

    #output the network data, set file names on command line
    -output-network


The Unix script `runImputeParallel`, as shown below, runs a number of BayesNetty processes in parallel and outputs separate data files for different subsets of individuals.
As the random number seed is set by default by the execution time, and the processes are set off at the same time, it is necessary to set the seed to different values.
The output files are then combined and the data files from separate processes deleted.


.. code-block:: none
      
    #!/bin/bash                                                                                                                                                                                       
    # $1 = parameter file to impute data in parallel
    # $2 = imputed data file name
    # $3 = no. of processes to run in parallel                                                                                                                                                       
    RANDOM=$$
    #run bayesnetty $3 times for X bootstraps each; processes run simultaneously in the background                                                                
    for i in $(seq 1 $3);
    do

    ./bayesnetty $1 -so -seed $i0$RANDOM -output-network-node-data-file-prefix $2$i-i -output-network-node-data-bed-file -output-network-node-data-job $i $3 -impute-network-data-job $i $3&

    done

    #wait for all processes to finish
    wait

    ##collate files                                                                                                                                                                                   
    if [ -f "$21-i-cts.dat" ]
    then
    > $2-cts.dat
    fi

    if [ -f "$21-i-discrete.dat" ]
    then
    > $2-discrete.dat
    fi

    for j in $(seq 1 $3);
    do

    #collate cts data
    if [ -f "$2$j-i-cts.dat" ]
    then
    cat $2$j-i-cts.dat >> $2-cts.dat
    rm $2$j-i-cts.dat
    fi

    #collate discrete data
    if [ -f "$2$j-i-discrete.dat" ]
    then
    cat $2$j-i-discrete.dat >> $2-discrete.dat
    rm $2$j-i-discrete.dat
    fi


    #collate SNP plink style data
    if [ -f "$2$j-i.fam" ]
    then

    if [ $j == 1 ]
    then
      cp $2$j-i.fam $2.fam
      cp $2$j-i.bim $2.bim
      cp $2$j-i.bed $2.bed
    else
      plink --noweb --silent --bfile $2 --bmerge $2$j-i.bed $2$j-i.bim $2$j-i.fam --make-bed --out $2-merge
      mv $2-merge.bed $2.bed
      mv $2-merge.bim $2.bim
      mv $2-merge.fam $2.fam
      rm $2-merge.log
    fi

    rm $2$j-i.fam
    rm $2$j-i.bim
    rm $2$j-i.bed
    fi

    done


The final imputed data can then be used in any BayesNetty analysis. For example, to search for the best fit network:


.. code-block:: none

    ./bayesnetty paras-search-imputed-data.txt


Where the parameter file is as follows:

.. code-block:: none

    #input imputed continuous data
    -input-data
    -input-data-file imputed-data-cts.dat
    -input-data-cts

    #input imputed SNP data as discrete data
    -input-data
    -input-data-file imputed-data.bed
    -input-data-discrete-snp

    #search network models with the imputed data
    -search-models


The files `paras-impute-parallel.txt`, `runImputeParallel` and `paras-search-imputed-data.txt` can be found in the `impute-example.zip <https://github.com/NewcastleRSE/BayesNetty/raw/refs/heads/main/docs/resources/impute-example.zip>`_ file.
