.. _input-network:

Input network
=============


A network may be specified using the ``-input-network`` task. The network may be used as a starting point for analyses, such as searches, or to perform an analysis on this network.
Only nodes in input files will be used in the network so that a subset of the data may be specified.


Any network **constraints** must be set using the ``-input-network`` option.
These constraints then belong to the network and will be used in any subsquent analysis, including searches, calculating average networks etc. 


If no network is specified then a network with no edges and a node for every data variable (as given by the input data) will be created and named "defaultNetwork".

.. _input-network-options:

Options
-------

The options are as follows:


.. list-table:: 
    :header-rows: 1

    * - Option
      - Description
      - Default

    * - -input-network
      - do a task to input a network
      -

    * - -input-network-name name
      - label the task and network with a name
      - Task-n

    * - -input-network-type t
      - the type of Bayesian network, choose between bnlearn or deal
      - bnlearn 

    * - -input-network-file network.dat
      - input the network in a format where the nodes and then the edges are listed
      -

    * - -input-network-file2 network2.dat
      - input the network in this style of format: ``[a][b|a][c|a:b]``
      -

    * - -input-network-igraph-file-prefix mygraph
      - input the network from igraph format files consisting of mygraph-nodes.dat and mygraph-edges.dat
      -

    * - -input-network-empty
      - set the network to one with no edges and one node for every data variable. An input network file is not required if this option is used
      -

    * - -input-network-whitelist-file whitelist.dat
      - a list of edges that must be included in any network
      -

    * - -input-network-blacklist-file blacklist.dat
      - a list of edges that must *not* be included in any network
      -

    * - -input-network-blacklist-edge-type dataName1 dataName2
      - edge types that may *not* be included in any network. The collection of nodes are given by the data input name, and so the data types must be given in different files
      -

    * - -input-network-no-parents-node nodeX
      - nodeX must not have any parents (except for white edges)
      -

    * - -input-network-no-children-node nodeY
      - nodeY must not have any children (except for white edges)
      -

    * - -input-network-prob-edge node1 node2 prob
      - set the prior probability of edge direction of node1 to node2 as prob
      -

    * - -input-network-prob-edge-type nodeType1 nodeType2 prob
      - set the prior probabilities of edge direction of nodeType1 to nodeType2 as prob
      -

    * - -input-network-imaginary-sample-size i
      - for deal networks this sets the imaginary sample size
      - 10

    * - -input-network-score score
      - for a bnlearn network choose between loglike, AIC or BIC
      - BIC


.. experimental option
.. -input-network-score-fix & fix for likelihood calculation when discrete data (or combinations) has too few points for linear regression, choose between none, average or skip & none

.. _input-network-black:

Black lists
-----------

A black list can be given using the ``-input-network-blacklist-file`` option to define a list of edges that must not be included in any network.
The text file should be formatted as follows:

.. code-block:: none

    node1 node2
    node2 node1
    node1 node3


such that the two nodes of each blacklisted edge are on one line. The nodes are ordered so the first line states that the edge node1 to node2 is not permitted.
The next line states that the edge in the reverse direction is also not permitted.

Any searches will ignore these blacklisted edges and attempting to use a network with a blacklisted edge will result in the edge being removed.

Edges between different types of nodes may also be blacklisted. This can be done using the ``-input-network-blacklist-edge-type`` option. It can be used as follows:

.. code-block:: none
        
    -input-data
    -input-data-name horses
    -input-data-file horses.dat
    -input-data-cts

    -input-data
    -input-data-name whips
    -input-data-file whips.dat
    -input-data-cts

    -input-network
    -input-network-name race
    -input-network-file model.dat
    -input-network-blacklist-edge-type horses whips 



Firstly the different node types must be loaded separately and given names using the ``-input-data-name`` option.
Then, when initially loading a network, the ``-input-network-blacklist-edge-type`` can be used to forbid any edge from one data set to another data set (or the same data if desired).
In the above example the network may not have any edge that goes from a horse to a whip, that is, a whip node may not have a horse node as a parent.
In any search that is performed these edges will not be considered. 



.. _input-network-white:

White lists
-----------

A white list can be given using the ``-input-network-whitelist-file`` option to define a list of edges that must be included in any network. The text file should be formatted as follows:

.. code-block:: none

    node1 node3
    node1 node2
    node2 node1


such that the two nodes of each whitelisted edge are on one line. The nodes are ordered so the first line states that the edge node1 to node3 must be included.
If both directions are included between two nodes then the edge must be included but may be in any direction.


If the whitelist and blacklist contradict one another then an error will be given.

.. _input-network-soft-con:

Soft Constraints
----------------

Soft constraints provide a way that the direction of an edge may be influenced but not with certainty, unlike blacklisted edges or whitelisted edges as described above.
An example parameter file setting a soft constraint, such that the prior probability of variable ``express`` to variable ``pheno`` is believed to be 0.8 is shown below.


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
    -input-network-name myNetwork
    -input-network-file example-network-format1.dat
    -input-network-prob-edge express pheno 0.8

    #search network models with the soft constraint
    -search-models



This parameter file, ``paras-soft-constraints.txt``, can be found in the file `example.zip <https://github.com/NewcastleRSE/BayesNetty/raw/refs/heads/main/docs/resources/example.zip>`_.


Any searches will use this prior probability.


If you wish to blacklist or whitelist an edge you should use those options rather than setting the prior probability to 0 or 1 for the sake of computational efficiency.

.. _input-network-formats: 

Network formats
---------------

The network may be defined using one of 3 different formats.

.. _input-network-formats-format1: 

Network file format 1
^^^^^^^^^^^^^^^^^^^^^

The first format is given by using the ``-input-network-file`` option and the network text file should be formatted as follows:

.. code-block:: none

    node1
    node2
    node3
    node2 node1
    node3 node1


where the nodes are listed first followed by the directed edges. In the above example there are 3 nodes and 2 edges, which are node2 to node1 and node3 to node1.

.. _input-network-formats-format2:

Network file format 2
^^^^^^^^^^^^^^^^^^^^^

The second format is given by using the ``-input-network-file2`` option and the network text file should be formatted as follows:

.. code-block:: none

    [node2][node3][node1|node2:node3]



where the nodes are listed in order of dependency. The independent nodes node2 and node3 are list first followed by node1 which is a child node of both node2 and node3.
This is the format that is typically output for searches and such like.


.. _input-network-formats-format3:

Network file format 3
^^^^^^^^^^^^^^^^^^^^^

The third format is given by using the ``-input-network-igraph-file-prefix`` option using the files that were output to draw the network in ``R``, see :ref:`igraph`.
There will be one file for the nodes and one for the edges, for example ``myNetwork-nodes.dat`` and ``myNetwork-edges.dat`` respectively. The node file will look something as follows:

.. code-block:: none

    id name type fileno
    1 node1 c 1
    2 node2 c 1
    3 node3 c 1


and the edges file will look like something as follows:


.. code-block:: none

    from to chisq
    2 1 6860.83
    3 1 5709.51


.. _input-network-example:

Example
-------

The following is an example parameter file to input a network.

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
    -input-network-name myNetwork
    -input-network-file example-network-format1.dat


This parameter file, ``paras-input-network.txt``, can be found in the file `example.zip <https://github.com/NewcastleRSE/BayesNetty/raw/refs/heads/main/docs/resources/example.zip>`_ and can be used as follows:

.. code-block:: none

    ./bayesnetty paras-input-network.txt


Which should produce output that looks like something as follows:

.. code-block:: none
        
    BayesNetty: Bayesian Network software, v1.00
    --------------------------------------------------
    Copyright 2015-present Richard Howey, GNU General Public License, v3
    Institute of Genetic Medicine, Newcastle University

    Random seed: 1551697141
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
    Task name: myNetwork
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

    Run time: 1 second



The data is loaded and then the network is loaded. The network has been named "myNetwork", and basic information about the network is output. 


Similarly, the network may be input using format 2 and 3 as given in parameter files ``paras-input-network2.txt`` and ``paras-input-network3.txt`` respectively.
