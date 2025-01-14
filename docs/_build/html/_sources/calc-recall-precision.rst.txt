.. _calc-recall-precision:

Calculate Recall and Precision
==============================

It is possible to calculate the recall and precision of a network against the "true" network structure.
Typically this option will be used when the true network structure is chosen and used to simulate data.
A best fit network can then be found for this data and the accuaracy assessed by calculating the recall and precision.  


For a network the recall is the percentage of edges found from the original true network.
The precision is the percentage of edges in the network that are also in the original true network. For an edge to be correct it must be in the correct direction.
However, if the true network has equivalent networks such that some edges may be in either direction, then these edges are considered correct if they are in either direction. 


.. _calc-recall-precision-options: 

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

    * - -calculate-recall-precision-name
      - label the task with a name
      - Task-n

    * - -calculate-recall-precision-network-name network1
      - the name of the network to calculate the recall and precision for
      -

    * - -calculate-recall-precision-true-network-name network2
      - the name of the true network to calculate the recall and precision against
      -

    * - -calculate-recall-precision-file file.dat
      - file to write recall and precision results to
      -



.. _calc-recall-precision-example: 

Example
-------


An example of calculating the recall and precision is contained in the parameter file `paras-example-calc-recall-precision.txt`,
which can be found in `example.zip <https://github.com/NewcastleRSE/BayesNetty/raw/refs/heads/main/docs/resources/example.zip>`_.
The network is again taken from the *code* child */code* network from the bnlearn repository *cite* bnlearn */cite*.

For example, the following parameter file:

.. code-block:: none

    #input the network to calculate the recall and precision of
    -input-network
    -input-network-name example-net
    -input-network-file example-net-calc-recall-pre.dat

    #input true network structure for these nodes
    -input-network
    -input-network-name true-net
    -input-network-file example-true-net-child.dat

    #calculate the recall and precision
    -calculate-recall-precision
    -calculate-recall-precision-network-name example-net
    -calculate-recall-precision-true-network-name true-net
    -calculate-recall-precision-file recall-precision.dat

can be ran in the usual way

.. code-block:: none

    ./bayesnetty paras-example-calc-recall-precision.txt

and will output something as follows

.. code-block:: none

    BayesNetty: Bayesian Network software, v1.1
    --------------------------------------------------
    Copyright 2015-present Richard Howey, GNU General Public License, v3
    Institute of Genetic Medicine, Newcastle University

    Random seed: 1605806239
    --------------------------------------------------
    Task name: example-net
    Loading network
    Network file: example-net-calc-recall-pre.dat
    Network type: bnlearn
    Network score type: BIC
    Total number of nodes: 20 (Discrete: 0 | Factor: 0 | Continuous: 0 | No data: 20)
    Total number of edges: 16
    Network Structure: [BirthAsphyxia][LVHreport][RUQO2][Sick][XrayReport][Age|Sick]
    [ChestXray|XrayReport][LVH|LVHreport][LowerBodyO2|LVHreport][CO2Report|Age]
    [HypDistrib|LowerBodyO2][LungFlow|ChestXray][CO2|CO2Report][DuctFlow|LungFlow]
    [Disease|DuctFlow][CardiacMixing|Disease][HypoxiaInO2|CardiacMixing][GruntingReport|Hy...
    The network has nodes with no data
    --------------------------------------------------
    --------------------------------------------------
    Task name: true-net
    Loading network
    Network file: example-true-net-child.dat
    Network type: bnlearn
    Network score type: BIC
    Total number of nodes: 20 (Discrete: 0 | Factor: 0 | Continuous: 0 | No data: 20)
    Total number of edges: 25
    Network Structure: [BirthAsphyxia][Disease|BirthAsphyxia][CardiacMixing|Disease][DuctFlow|Disease]
    [LVH|Disease][LungFlow|Disease][LungParench|Disease][Sick|Disease][Age|Disease:Sick][CO2|LungParench]
    [ChestXray|LungFlow:LungParench][Grunting|LungParench:Sick][HypDistrib|CardiacMixing:DuctFlow][HypoxiaInO2|CardiacMixing...
    The network has nodes with no data
    --------------------------------------------------
    --------------------------------------------------
    Task name: Task-3
    Calculating the recall and precision
    Network: example-net
    Network Structure: [BirthAsphyxia][LVHreport][RUQO2][Sick][XrayReport][Age|Sick][ChestXray|XrayReport]
    [LVH|LVHreport][LowerBodyO2|LVHreport][CO2Report|Age][HypDistrib|LowerBodyO2][LungFlow|ChestXray]
    [CO2|CO2Report][DuctFlow|LungFlow][Disease|DuctFlow][CardiacMixing|Disease][HypoxiaInO2|CardiacMixing][GruntingReport|Hy...
    True Network: true-net
    True Network Structure: [BirthAsphyxia][Disease|BirthAsphyxia][CardiacMixing|Disease][DuctFlow|Disease]
    [LVH|Disease][LungFlow|Disease][LungParench|Disease][Sick|Disease][Age|Disease:Sick][CO2|LungParench]
    [ChestXray|LungFlow:LungParench][Grunting|LungParench:Sick][HypDistrib|CardiacMixing:DuctFlow][HypoxiaInO2|CardiacMixing...
    Recall and precision written to file: recall-precision.dat

    Recall: the percentage of edges found from the original true network
    Precision: the percentage of edges in the network that are also in the original true network

    Recall: 32
    Precision: 50

    Recall and precision written to file: recall-precision.dat
    --------------------------------------------------

    Run time: less than one second


In this example the network for which we wish we calculate the recall and precision is input into BayesNetty.
The true network structure is then also input into BayesNetty. Finally the recall and precision is calculated and the results output to a file.
This file simply contains 2 numbers: the recall followed by the precision.
