.. _search-models:

Search models
=============

Network models can be searched for one that best describes the data as given by the network model assumptions, network constraints, the network score and the data itself.
The search option uses a network to start the search and finishes with it updated to the found best fit network.
(If a network is not set then a default network is used.) Any constraints to the search must be set when the network is setup, see :ref:`input-network`.

.. _search-models-options: 

Options
-------

The options are as follows:

.. list-table:: 
    :header-rows: 1

    * - Option
      - Description
      - Default

    * - -search-models
      - do a task to search network models
      -

    * - -search-models-name name
      - label the task with a name
      - Task-n

    * - -search-models-network-name network
      - the name of the network to start the search from
      - previous network (or the default model given by a node for each data variable and no edges if there is no previous network) 

    * - -search-models-file search.dat
      - record the network models and scores in the search path to file search.dat
      -

    * - -search-models-random-restarts n
      - do another n searches starting from a random network
      - 0

    * - -search-models-jitter-restarts m
      - after the initial search and every random restart search do another m searches jittered from the recently found network
      - 0


.. _search-models-greedy:

Greedy search
-------------


The greedy search algorithm is the default algorithm for searching through network models, and is currently the only search algorithm.

.. _search-models-greedy-restart: 

Number of random restarts for the greedy algorithm
--------------------------------------------------

The greedy algorithm can be ran a further number of times from a random starting network. The number of random restarts is set by using the option `-search-models-random-restarts`.

.. _search-models-greedy-jitter: 

Number of jitter restarts for the greedy algorithm
--------------------------------------------------

Once the greedy algorithm has converged on a final best fit network, the algorithm can be restarted at a network given by slightly modifying the best fit network, also called *jittering*.
This may be useful to avoid the algorithm sticking in a local maximum whilst still retaining more or less the same network.
The number of times times the search should be jittered is set by using the option `-search-models-jitter-restarts`.


Random restarts and jittered restarts can be used together, if there are :math:`n` random restarts and :math:`m` jittered restarts then there will be :math:`(n + 1)` times :math:`m` searches. 


.. _search-models-example: 

Example
-------

As an example of searching through network models the parameter file `paras-search.txt`,
which can be found in `example.zip <https://github.com/NewcastleRSE/BayesNetty/raw/refs/heads/main/docs/resources/example.zip>`_,
searches through network models starting from the default model given by a node for each data variable and no edges.


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

    #search network models
    -search-models


This can be executed as usual

.. code-block:: none

    ./bayesnetty paras-search.txt


and will output something as follows

.. code-block:: none

    BayesNetty: Bayesian Network software, v1.00
    --------------------------------------------------
    Copyright 2015-present Richard Howey, GNU General Public License, v3
    Institute of Genetic Medicine, Newcastle University

    Random seed: 1551700554
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

    Run time: less than one second


The above shows the data input and then the default network input consisting of a node for each data variable given by the data and no edges.
The network with the highest network score is shown in the output.

