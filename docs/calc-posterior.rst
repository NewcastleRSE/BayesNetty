.. _calc-posterior:

Calculate posterior
===================

The likelihood is calculated under specific distributional assumptions,
namely that discrete nodes follow a multinomial distribution and continuous nodes a normal distribution,
with distributional parameters determined by the values of the incoming parent nodes. We call the network structure with its corresponding distributional parameters the posterior.

If only the posterior is of interest then the `-calc-posterior` option can be used without the need to perform any other analyses.
One would probably want to also use the `-output-posteriors` option to output the posteriors, see :ref:`output-posteriors`.


.. _calc-posterior-options: 

Options
=======

The options are as follows:

.. list-table:: 
    :header-rows: 1

    * - Option
      - Description
      - Default

    * - -calc-posterior
      - do a task to calculate the posterior
      -

    * - -calc-posterior-name
      - label the task with a name
      - Task-n

    * - -calc-posterior-network-name network
      - the name of the network to calculate the posterior of
      - previous network (or the default model given by a node for each data variable and no edges if there is no previous network) 


.. _calc-post-example: 

Example
-------

The posterior is calculated by simply using the `-calc-posterior` option in the parameter file after the data and network has been set up.
For example:

.. list-table:: 
    :header-rows: 1

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

    #calculate the posterior of the network
    -calc-posterior

Note that the network has not been set for the `-calc-posterior` task as there is only one network, and so by default the most recent network is used.
This parameter file, `paras-calc-post.txt`, can be found in `example.zip <https://github.com/NewcastleRSE/BayesNetty/raw/refs/heads/main/docs/resources/example.zip>`_ and produces the following output:

.. code-block:: none    

    BayesNetty: Bayesian Network software, v1.00
    --------------------------------------------------
    Copyright 2015-present Richard Howey, GNU General Public License, v3
    Institute of Genetic Medicine, Newcastle University

    Random seed: 1551697290
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
    Calculating posterior
    Network: Task-4
    Network Structure: [mood][rs1][rs2][pheno|rs1:rs2][express|pheno:mood]
    --------------------------------------------------

    Run time: less than one second

For an example of calculating and outputting the posterior to file, see :ref:`output-posts-example`.
