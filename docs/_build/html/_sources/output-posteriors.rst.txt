
.. _output-posteriors:

Output posteriors
=================


The posteriors may be output to file for inspection with the `-output-posteriors`.

.. _output-posteriors-options:

Options
-------

The options are as follows:


.. list-table:: 
    :header-rows: 1

    * - Option
      - Description
      - Default

    * -  -output-posteriors
      - do a task to output the posteriors of a network to file
      -

    * - -output-posteriors-name name
      - label the task with a name
      - Task-n

    * - -output-posteriors-network-name network
      - output posteriors for this network
      - previous network (or the default model given by a node for each data variable and no edges if there is no previous network)

    * - -output-posteriors-file posts.dat
      - output the posteriors to file posts.dat
      - posteriors.dat



.. _output-posts-example:

Example
-------

The following is an example parameter file to output the posteriors of a network.

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

    #calculate the posterior of the network
    -calc-posterior

    #output the posteriors to file
    -output-posteriors
    -output-posteriors-file example-posteriors.dat


This parameter file, `paras-output-post.txt`, can be found in `example.zip <https://github.com/NewcastleRSE/BayesNetty/raw/refs/heads/main/docs/resources/example.zip>`_ and can be used as follows:


.. code-block:: none

    ./bayesnetty paras-output-post.txt


Which should produce output that looks like something as follows:

.. code-block:: none

    BayesNetty: Bayesian Network software, v1.00
    --------------------------------------------------
    Copyright 2015-present Richard Howey, GNU General Public License, v3
    Institute of Genetic Medicine, Newcastle University

    Random seed: 1551958097
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
    --------------------------------------------------
    Task name: Task-6
    Outputting posteriors
    Network: Task-4
    Network Structure: [mood][rs1][rs2][pheno|rs1:rs2][express|pheno:mood]
    Output posteriors to file: example-posteriors.dat
    --------------------------------------------------

    Run time: less than one second


The data is loaded, the network input, the posterior is calculated and then output to a file.

