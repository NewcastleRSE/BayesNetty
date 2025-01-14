.. _calc-score: 

Calculate network score
=======================

The network score is used as a measure of how well the network model describes the data and is used to compare different models when searching through models.
In BayesNetty the network score is based on the log likelihood and higher values imply a better model
(see :ref:`bnlearn-score` for further details). This is calculated assuming that discrete nodes follow a multinomial distribution and continuous nodes a normal distribution.
BayesNetty considers the network score to be a property of the network and its method of calculation is set using the
option `-input-network-score`, see :ref:`input-network`.

.. _calc-score-options:

Options
-------

The options are as follows:

.. list-table:: 
    :header-rows: 1

    * - Option
      - Description
      - Default

    * - -calc-network-score
      - do a task to calculate the score
      -

    * - -calc-network-score-name name
      - label the task with a name
      - Task-n

    * - -calc-network-score-network-name network
      - the name of the network to calculate the score of
      - previous network (or the default model given by a node for each data variable and no edges if there is no previous network) 

    * - -calc-network-score-file file
      - write the score to this file
      - 

    * - -calc-network-score-all-scores network-scores.dat
      - calculate the scores of *every* possible network and record the results in `network-scores.dat`
      - 

  
.. _calc-score-example:

Example
-------

As an example of calculating the score the parameter file `paras-calc-score.txt`, which can be found in `example.zip <https://github.com/NewcastleRSE/BayesNetty/raw/refs/heads/main/docs/resources/example.zip>`_,
calculates the score for the same network but for different score methods.


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
    -input-network-name networkLike
    -input-network-score loglike
    -input-network-file example-network-format1.dat

    #input the example network in format 1
    -input-network
    -input-network-name networkBIC
    -input-network-score BIC
    -input-network-file example-network-format1.dat

    #calculate the network of the network with BIC
    -calc-network-score

    #calculate the network of the network with log likelihood
    -calc-network-score
    -calc-network-score-network-name networkLike

This can be executed as usual

.. code-block:: none

    ./bayesnetty paras-calc-score.txt

and will output something as follows

.. code-block:: none

    BayesNetty: Bayesian Network software, v1.00
    --------------------------------------------------
    Copyright 2015-present Richard Howey, GNU General Public License, v3
    Institute of Genetic Medicine, Newcastle University

    Random seed: 1551700452
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
    Task name: networkLike
    Loading network
    Network file: example-network-format1.dat
    Network type: bnlearn
    Network score type: log likelihood
    Total number of nodes: 5 (Discrete: 3 | Factor: 0 | Continuous: 2)
    Total number of edges: 4
    Network Structure: [mood][rs1][rs2][pheno|rs1:rs2][express|pheno:mood]
    Total data at each node: 1495
    Missing data at each node: 5
    --------------------------------------------------
    --------------------------------------------------
    Task name: networkBIC
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
    Task name: Task-6
    Calculating network score
    Network: networkBIC
    Network structure: [mood][rs1][rs2][pheno|rs1:rs2][express|pheno:mood]
    Network score type: BIC
    Network score = -8519.74
    --------------------------------------------------
    --------------------------------------------------
    Task name: Task-7
    Calculating network score
    Network: networkLike
    Network structure: [mood][rs1][rs2][pheno|rs1:rs2][express|pheno:mood]
    Network score type: log likelihood
    Network score = -8413.75
    --------------------------------------------------

    Run time: less than one second

The above output shows the data input and then two networks input with the same structure but with different scores.
The network with the BIC score is evaluated firstly, as by default the most recent network is used unless otherwise stated.
The network using the log likelihood is then calculated by using the `-calc-network-score-network-name` option to specify which network should be used.
