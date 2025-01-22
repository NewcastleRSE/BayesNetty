.. _markov-blanket:

Markov blanket
==============


The Markov blanket for a node contains all the variables that shield the node from the rest of the network.
This means that the Markov blanket of a node is the only knowledge needed to predict the behaviour of that node and its children.
This may be useful for large networks where some nodes are of particular interest. See :cite:`bnlearn_paper` for more details.

The `-markov-blanket` option can be used to calculate the Markov blanket for a given node.
A sub-network is created for the given node and its Markov blanket which may be output using the `-output-network` option, see :ref:`output-network`, or used in any other network analysis.

.. _markov-blanket-options: 

Options
-------

The options are as follows:

.. list-table:: 
    :header-rows: 1

    * - Option
      - Description
      - Default

    * - -markov-blanket
      - do a task to calculate the Markov blanket for a given node
      -

    * - -markov-blanket-name name
      - label the task with a name
      - Task-n

    * - -markov-blanket-network-name network
      - the name of the network to calculate the Markov blanket
      - previous network (or the default model given by a node for each data variable and no edges if there is no previous network)

    * - -markov-blanket-node-name node
      - calculate the Markov blanket for the node with this name
      -


.. _calc-blanket-example:

Example
-------

The Markov blanket for a given node is calculated by using the `-markov-blanket` option together with the `-markov-blanket-node-name` option to choose the node. For example:


.. code-block:: none

    #input the example network
    -input-network
    -input-network-file example-network-format1.dat

    #calculate the Markov Blanket
    -markov-blanket
    -markov-blanket-node-name express

    #output the network
    -output-network
    -output-network-file express-markov-blanket.dat


As the Markov blanket does not depend on the data of the nodes it is possible to calculate the Markov blanket without data.
In the example the network structure is input, the Markov blanket calculated for the "express" node and then output to the file `express-markov-blanket.dat`.
This parameter file, `paras-blanket.txt`, can be found in `example.zip <https://github.com/NewcastleRSE/BayesNetty/raw/refs/heads/main/docs/resources/example.zip>`_ and produces output which should look something as follows:


.. code-block:: none

    BayesNetty: Bayesian Network software, v1.00
    --------------------------------------------------
    Copyright 2015-present Richard Howey, GNU General Public License, v3
    Institute of Genetic Medicine, Newcastle University

    Random seed: 1551956198
    --------------------------------------------------
    Task name: Task-1
    Loading network
    Network file: example-network-format1.dat
    Network type: bnlearn
    Network score type: BIC
    Total number of nodes: 5 (Discrete: 0 | Factor: 0 | Continuous: 0 | No data: 5)
    Total number of edges: 4
    Network Structure: [rs1][rs2][mood][pheno|rs1:rs2][express|mood:pheno]
    The network has nodes with no data
    --------------------------------------------------
    --------------------------------------------------
    Task name: Task-2
    Calculating Markov blanket
    Network: Task-1
    Node: express
    Network structure: [rs1][rs2][mood][pheno|rs1:rs2][express|mood:pheno]
    Markov blanket network structure: [mood][pheno][express|mood:pheno]
    --------------------------------------------------
    --------------------------------------------------
    Task name: Task-3
    Outputting network
    Network: Task-2
    Network Structure: [mood][pheno][express|mood:pheno]
    Network output to file: express-markov-blanket.dat
    --------------------------------------------------

    Run time: less than one second
