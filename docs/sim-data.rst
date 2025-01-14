.. _sim-data:

Simulate network data
=====================

It is possible to simulate data for a given network using BayesNetty. The network must be a bnlearn network, see :ref:`bnlearn`,
and can be set using a network file which has the same format as a posterior file and sets the network structure and the network parameters.
Alternatively, rather than reading in a posterior file, the simulating network can be set by first (in the same parameter file) carrying out any analysis task in BayesNetty that results in a bnlearn network with calculated posteriors.

.. _sim-data-options: 

Options
-------

The options are as follows:

.. list-table:: 
    :header-rows: 1

    * - Option
      - Description
      - Default

    * - -calculate-recall-precision
      - do a task to calculate the recall and precision
      -

    * - -simulate-network-data
      - do a task to simulate network data for a given network
      -

    * - -simulate-network-data-name name
      - label the task with a name
      - Task-n

    * - -simulate-network-data-network-name network
      - simulate data for this network & previous network (or the default model given by a node for each data variable and no edges if there is no previous network)
      -

    * - -simulate-network-data-no-sims n
      - simulate n replicates of network data
      - 100

    * - -simulate-network-data-parameter-file parameters.txt
      - network with parameters in bnlearn posteriors file format
      -

    * - -simulate-network-data-whitelist-file whitelist.dat
      - a list of edges that must be included in any network
      -

    * - -simulate-network-data-blacklist-file blacklist.dat
      - a list of edges that must *not* be included in any network
      -

    * - -simulate-network-data-score score
      - for a bnlearn network choose between loglike, AIC or BIC & BIC
      -


If posteriors are not given then default network node parameters are used to simulate data.
 The default probabilities for a discrete SNP node are 0.25, 0.5 and 0.25 for levels 0, 1 and 2 respectively. A minor allele frequency of 0.5 would give these probabilites.
 For a discrete node the levels are given by equal probabilities. For a continuous node the intercept is set to 10 and the coefficients and variance are set to 1.
 **Note:** for edges from discrete variables to continuous variables there will be no effect by default, as the continuous node will be given the same intercept and coefficients for 
 every different level configuration of the discrete variables.

The following sections show some examples of simulating data and can all be found in `example.zip <https://github.com/NewcastleRSE/BayesNetty/raw/refs/heads/main/docs/resources/example.zip>`_.


.. _sim-data-example1: 

Example 1: no data
------------------

This example simulates node data using a given network structure with network parameters.


The network structure with network parameters must be given to simulate network data.
These can be set by using a network parameters file which takes the same format as a posterior file and may be quite complex.
To create this network parameter file it is therefore recommended to first output a posterior file for the network that you wish to simulate data for and then to edit this posterior file as required.


To create a network parameter file, `example-network-parameters.txt`, run the parameter file, `paras-make-network-paras.txt`:


.. code-block:: none

    #input the example network in format 1
    -input-network
    -input-network-file example-network-sim.dat

    #simulate data using default parameters
    -simulate-network-data

    #calculate the posterior of the network using default parameters
    -calc-posterior

    #output the posteriors to be used as a network parameters file
    -output-posteriors
    -output-posteriors-file example-network-parameters.txt


The example network file, `example-network-sim.dat`, is in network format 1, see :ref:`input-network-formats`, except that extra (optional) node information is given as there is no node data to determine the node type.
After each node either `|dis.snp`, `|cts.snp`, `|dis|2` or `|cts` may be written to give a discrete SNP node, a continuous SNP node, a discrete node with the number of levels or a continuous node respectively.
If the type of node is not set then the node is set to a continuous node. If a discrete node is specified without the number of levels then the number of levels is set to 0. The example network file is as follows.


.. code-block:: none

    rs1|dis.snp
    rs2|dis.snp
    mood|dis|2
    express|cts
    pheno|cts
    rs1 pheno
    rs2 pheno
    mood express
    pheno express


The network parameter file, `example-network-parameters.txt`, can then be created by running the parameter file

.. code-block:: none

    ./bayesnetty paras-make-network-paras.txt


The output will look something as follows


.. code-block:: none

    BayesNetty: Bayesian Network software, v1.00
    --------------------------------------------------
    Copyright 2015-present Richard Howey, GNU General Public License, v3
    Institute of Genetic Medicine, Newcastle University

    Random seed: 1551712417
    --------------------------------------------------
    Task name: Task-1
    Loading network
    Network file: example-network-sim.dat
    Network type: bnlearn
    Network score type: BIC
    Total number of nodes: 5 (Discrete: 0 | Factor: 0 | Continuous: 0 | No data: 5)
    Total number of edges: 4
    Network Structure: [rs1][rs2][mood][pheno|rs1:rs2][express|mood:pheno]
    The network has nodes with no data
    --------------------------------------------------
    --------------------------------------------------
    Task name: Task-2
    Simulating network data
    Data simulation network given by network: Task-1
    Number of simulations: 100
    Network: Task-2
    Network type: bnlearn
    Network score type: BIC
    Total number of nodes: 5 (Discrete: 3 | Factor: 0 | Continuous: 2)
    Total number of edges: 4
    Network Structure: [rs1][rs2][mood][pheno|rs1:rs2][express|mood:pheno]
    Total data at each node: 100
    Missing data at each node: 0
    --------------------------------------------------
    --------------------------------------------------
    Task name: Task-3
    Calculating network score
    Network: Task-2
    Network structure: [rs1][rs2][mood][pheno|rs1:rs2][express|mood:pheno]
    Network score type: BIC
    Network score = -600.611
    --------------------------------------------------
    --------------------------------------------------
    Task name: Task-4
    Outputting posteriors
    Network: Task-2
    Network Structure: [rs1][rs2][mood][pheno|rs1:rs2][express|mood:pheno]
    Output posteriors to file: example-network-parameters.txt
    --------------------------------------------------

    Run time: less than one second


The network parameter file, `example-network-parameters.txt`, will look something like as follows:

.. code-block:: none

    Posteriors:
    ===========

    DISCRETE SNP NODE: rs1
      0: 0.26
      1: 0.49
      2: 0.25

    DISCRETE SNP NODE: rs2
      0: 0.25
      1: 0.51
      2: 0.24

    DISCRETE NODE: mood
      0: 0.47
      1: 0.53

    CONTINUOUS NODE: express
    DISCRETE PARENTS: mood
    0:
      Intercept: 9.71024
      Coefficients: pheno: 1.04254 
      Mean: 20.2202
      Variance: 0.811076
    1:
      Intercept: 9.89
      Coefficients: pheno: 1.00832 
      Mean: 19.9452
      Variance: 0.6773

    CONTINUOUS NODE: pheno
    DISCRETE PARENTS: rs1:rs2
    0:0:
      Intercept: 10.4619
      Coefficients: 
      Mean: 10.4619
      Variance: 2.07491
    1:0:
      Intercept: 9.97966
      Coefficients: 
      Mean: 9.97966
      Variance: 1.7134
    2:0:
      Intercept: 9.58534
      Coefficients: 
      Mean: 9.58534
      Variance: 0.532354
    0:1:
      Intercept: 10.0514
      Coefficients: 
      Mean: 10.0514
      Variance: 1.44765
    1:1:
      Intercept: 9.96623
      Coefficients: 
      Mean: 9.96623
      Variance: 0.760351
    2:1:
      Intercept: 9.54098
      Coefficients: 
      Mean: 9.54098
      Variance: 0.559675
    0:2:
      Intercept: 10.7559
      Coefficients: 
      Mean: 10.7559
      Variance: 1.78883
    1:2:
      Intercept: 10.1066
      Coefficients: 
      Mean: 10.1066
      Variance: 0.983235
    2:2:
      Intercept: 10.3842
      Coefficients: 
      Mean: 10.3842
      Variance: 1.50665


The data was simulated using default network node parameters. The simulated node data and subsequent fitted parameters are thus close to these values.


The network parameter file, *code* example-network-parameters.txt */code*, can now be edited using parameters of your choice.
The mean is not required, this simply reports the mean of the node data for continuous nodes. The "Posteriors" title in the file is also not required, but these may be left in the file.
The levels of discrete nodes are labelled 0, 1, 2 etc. These may be renamed to something more meaningful in this file. For example, for node "mood" the levels could be renamed "sad" and "happy".

Finally, network node data may be simulated for a network with chosen network parameters where initially there was no data available for the network.


.. code-block:: none

    #simulate data
    -simulate-network-data
    -simulate-network-data-no-sims 200
    -simulate-network-data-parameter-file example-network-parameters.txt

    #output simulated data
    -output-network
    -output-network-node-data-file-prefix sim-data
    -output-network-node-data-bed-file


The parameter file shown above, `paras-sim-data1.txt`, can then be used to simulate data and output it to several data files.


.. code-block:: none

    ./bayesnetty paras-sim-data1.txt


The output will look something as follows


.. code-block:: none

    BayesNetty: Bayesian Network software, v1.00
    --------------------------------------------------
    Copyright 2015-present Richard Howey, GNU General Public License, v3
    Institute of Genetic Medicine, Newcastle University

    Random seed: 1551712579
    --------------------------------------------------
    Task name: Task-1
    Simulating network data
    Number of simulations: 200
    Parameter file name: example-network-parameters.txt
    Network: Task-1
    Network type: bnlearn
    Network score type: BIC
    Total number of nodes: 5 (Discrete: 3 | Factor: 0 | Continuous: 2)
    Total number of edges: 4
    Network Structure: [rs1][rs2][mood][pheno|rs1:rs2][express|mood:pheno]
    Total data at each node: 200
    Missing data at each node: 0
    --------------------------------------------------
    --------------------------------------------------
    Task name: Task-2
    Outputting network
    Network: Task-1
    Network Structure: [rs1][rs2][mood][pheno|rs1:rs2][express|mood:pheno]
    Network output to file: network.dat
    Node data output to files:
    sim-data-discrete.dat
    sim-data-cts.dat
    sim-data.bed/.bim/.fam
    --------------------------------------------------

    Run time: less than one second


The file `sim-data-discrete.dat` contains the discrete node data, `sim-data-cts.dat` the continuous node data and `sim-data.bed`, `sim-data.bim` and `sim-data.fam` the SNP node data in PLINK binary pedigree format.



.. _sim-data-example2:

Example 2: data and fitted network
----------------------------------


This example inputs some data, sets the network structure, fits network posterior parameters and then simulates network node data using these parameters.
The simulated data is then output to file. The parameter file, `paras-sim-data2.txt`, in the example files does this and is as follows:

.. code-block:: none

    #input continuous data
    -input-data
    -input-data-file example-cts.dat
    -input-data-cts

    #input discrete data
    -input-data
    -input-data-file example-discrete.dat
    -input-data-discrete

    #input SNP data as discrete data
    -input-data
    -input-data-file example.bed
    -input-data-discrete-snp

    #input the example network in format 1
    -input-network
    -input-network-file example-network-format1.dat

    #simulate data
    -simulate-network-data
    -simulate-network-data-no-sims 200

    #output simulated data
    -output-network
    -output-network-node-data-file-prefix sim-data
    -output-network-node-data-bed-file



.. code-block:: none

    ./bayesnetty paras-sim-data2.txt


The output will look something as follows

.. code-block:: none

    BayesNetty: Bayesian Network software, v1.00
    --------------------------------------------------
    Copyright 2015-present Richard Howey, GNU General Public License, v3
    Institute of Genetic Medicine, Newcastle University

    Random seed: 1551716397
    --------------------------------------------------
    Task name: Task-1
    Loading data
    Continuous data file: example-cts.dat
    Number of ID columns: 2
    Including (all) 2 variables in analysis
    Each variable has 1500 data entries
    Missing value: not set
    --------------------------------------------------
    --------------------------------------------------
    Task name: Task-2
    Loading data
    Discrete data file: example-discrete.dat
    Number of ID columns: 2
    Including the 1 and only variable in analysis
    Each variable has 1500 data entries
    Missing value: NA
    --------------------------------------------------
    --------------------------------------------------
    Task name: Task-3
    Loading data
    SNP binary data file: example.bed
    SNP data treated as discrete data
    Total number of SNPs: 2
    Total number of subjects: 1500
    Number of ID columns: 2
    Including (all) 2 variables in analysis
    Each variable has 1500 data entries
    --------------------------------------------------
    --------------------------------------------------
    Task name: Task-4
    Loading network
    Network file: example-network-format1.dat
    Network type: bnlearn
    Network score type: BIC
    Total number of nodes: 5 (Discrete: 3 | Factor: 0 | Continuous: 2)
    Total number of edges: 4
    Network Structure: [mood][rs1][rs2][pheno|rs1:rs2][express|pheno:mood]
    Total data at each node: 1495
    Missing data at each node: 5
    --------------------------------------------------
    --------------------------------------------------
    Task name: Task-5
    Simulating network data
    Data simulation network given by network: Task-4
    Number of simulations: 200
    Network: Task-5
    Network type: bnlearn
    Network score type: BIC
    Total number of nodes: 5 (Discrete: 3 | Factor: 0 | Continuous: 2)
    Total number of edges: 4
    Network Structure: [mood][rs1][rs2][pheno|rs1:rs2][express|pheno:mood]
    Total data at each node: 200
    Missing data at each node: 0
    --------------------------------------------------
    --------------------------------------------------
    Task name: Task-6
    Outputting network
    Network: Task-5
    Network Structure: [mood][rs1][rs2][pheno|rs1:rs2][express|pheno:mood]
    Network output to file: network.dat
    Node data output to files:
    sim-data-discrete.dat
    sim-data-cts.dat
    sim-data.bed/.bim/.fam
    --------------------------------------------------

    Run time: less than one second



As in the previous example the simulated node data is output to a number of files.

.. _sim-data-example3:

Example 3: data and unknown network
-----------------------------------

In this example the network that the data will be simulated for is not known initially.
To do this, some data is input, the best fitting network is chosen using a network search and then node data is simulated using the fitted parameters.
The parameter file, `paras-sim-data3.txt`, in the example files does this and is as follows: 


.. code-block:: none

    #input continuous data
    -input-data
    -input-data-file example-cts.dat
    -input-data-cts

    #input discrete data
    -input-data
    -input-data-file example-discrete.dat
    -input-data-discrete

    #input SNP data as discrete data
    -input-data
    -input-data-file example.bed
    -input-data-discrete-snp

    #search for the best fitting model
    -search-models

    #simulate data
    -simulate-network-data
    -simulate-network-data-no-sims 200

    #output simulated data
    -output-network
    -output-network-node-data-file-prefix sim-data
    -output-network-node-data-bed-file



.. code-block:: none

    ./bayesnetty paras-sim-data3.txt


The output will look something as follows


.. code-block:: none

    BayesNetty: Bayesian Network software, v1.00
    --------------------------------------------------
    Copyright 2015-present Richard Howey, GNU General Public License, v3
    Institute of Genetic Medicine, Newcastle University

    Random seed: 1551716718
    --------------------------------------------------
    Task name: Task-1
    Loading data
    Continuous data file: example-cts.dat
    Number of ID columns: 2
    Including (all) 2 variables in analysis
    Each variable has 1500 data entries
    Missing value: not set
    --------------------------------------------------
    --------------------------------------------------
    Task name: Task-2
    Loading data
    Discrete data file: example-discrete.dat
    Number of ID columns: 2
    Including the 1 and only variable in analysis
    Each variable has 1500 data entries
    Missing value: NA
    --------------------------------------------------
    --------------------------------------------------
    Task name: Task-3
    Loading data
    SNP binary data file: example.bed
    SNP data treated as discrete data
    Total number of SNPs: 2
    Total number of subjects: 1500
    Number of ID columns: 2
    Including (all) 2 variables in analysis
    Each variable has 1500 data entries
    --------------------------------------------------
    --------------------------------------------------
    Task name: Task-4
    Searching network models
    --------------------------------------------------
    Loading defaultNetwork network
    Network type: bnlearn
    Network score type: BIC
    Total number of nodes: 5 (Discrete: 3 | Factor: 0 | Continuous: 2)
    Total number of edges: 0
    Network Structure: [express][pheno][mood][rs1][rs2]
    Total data at each node: 1495
    Missing data at each node: 5
    --------------------------------------------------
    Network: defaultNetwork
    Search: Greedy
    Random restarts: 0
    Random jitter restarts: 0
    Network Structure: [mood][rs1][rs2][express|rs1:rs2][pheno|express:mood]
    Network score type: BIC
    Network score = -8213.45
    --------------------------------------------------
    --------------------------------------------------
    Task name: Task-5
    Simulating network data
    Data simulation network given by network: defaultNetwork
    Number of simulations: 200
    Network: Task-5
    Network type: bnlearn
    Network score type: BIC
    Total number of nodes: 5 (Discrete: 3 | Factor: 0 | Continuous: 2)
    Total number of edges: 4
    Network Structure: [mood][rs1][rs2][express|rs1:rs2][pheno|express:mood]
    Total data at each node: 200
    Missing data at each node: 0
    --------------------------------------------------
    --------------------------------------------------
    Task name: Task-6
    Outputting network
    Network: Task-5
    Network Structure: [mood][rs1][rs2][express|rs1:rs2][pheno|express:mood]
    Network output to file: network.dat
    Node data output to files:
    sim-data-discrete.dat
    sim-data-cts.dat
    sim-data.bed/.bim/.fam
    --------------------------------------------------

    Run time: less than one second


As in the previous example the simulated node data is output to a number of files.

