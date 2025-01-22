.. _output-priors:

Output priors (deal only)
=========================

The priors for a deal network may be output to file for inspection with the `-output-priors`.
The deal Bayesian network model has a quite complex default prior which is based on the given network data, structure and imaginary sample size, see :cite:`deal_paper` for details.
The bnlearn Bayesian network, which is the recommended and default Bayesian network model, has no prior to output, see :cite:`bnlearn_paper` for details.

.. _output-priors-options:

Options
-------

The options are as follows:


.. list-table:: 
    :header-rows: 1

    * - Option
      - Description
      - Default

    * - -output-priors
      - do a task to output the priors of a network to file
      - 

    * - -output-priors-name name
      - label the task with a name
      - Task-n

    * - -output-priors-network-name network
      - output priors for this network
      - previous network (or the default model given by a node for each data variable and no edges if there is no previous network)

    * - -output-priors-file priors.dat
      - output the priors to file priors.dat
      - priors.dat


.. _output-priors-example: 

Example
-------

The following is an example parameter file to output the priors of a network.

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
    -input-network-type deal

    #output the priors to file
    -output-priors
    -output-priors-file example-priors.dat


This parameter file, `paras-output-priors.txt`, can be found in `example.zip <https://github.com/NewcastleRSE/BayesNetty/raw/refs/heads/main/docs/resources/example.zip>`_ and can be used as follows:


.. code-block:: none

    ./bayesnetty paras-output-priors.txt


Which should produce output that looks like something as follows:

.. code-block:: none

    BayesNetty: Bayesian Network software, v1.00
    --------------------------------------------------
    Copyright 2015-present Richard Howey, GNU General Public License, v3
    Institute of Genetic Medicine, Newcastle University

    Random seed: 1551957572
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
    Network type: deal
    Total number of nodes: 5 (Discrete: 3 | Factor: 0 | Continuous: 2)
    Total number of edges: 4
    Network Structure: [mood][rs1][rs2][pheno|rs1:rs2][express|pheno:mood]
    Imaginary sample size: 10
    Total data at each node: 1495
    Missing data at each node: 5
    --------------------------------------------------
    --------------------------------------------------
    Task name: Task-5
    Outputting priors
    Network: Task-4
    Network Structure: [mood][rs1][rs2][pheno|rs1:rs2][express|pheno:mood]
    Output priors to file: example-priors.dat
    --------------------------------------------------

    Run time: less than one second



The data is loaded, the network input and then the prior is output to a file.

