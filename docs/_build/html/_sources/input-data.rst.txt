.. _input-data:

Input data
==========

All data must be input using the ``-input-data`` task.

.. _input-data-options:

Options
-------

The options are as follows:


.. list-table:: 
   :header-rows: 1

   * - Option
     - Description
     - Default

   * - -input-data
     - do a task to input data
     -

   * - -input-data-name name
     - label the task with a name
     - Task-n

   * - -impute-network-data-network-name network
     - the name of the network to impute data for
     - previous network (or the default model if there is no previous network) given by a node for each data variable and no edges

   * - -input-data-include-file nodes.dat
     - a list of nodes/variables from the data file to be included in the network. Only the nodes in this list will be used in any analysis.
     -

   * - -input-data-exclude-file nodes.dat
     - a list of nodes/variables to be excluded from the network
     -

   * - -input-data-cts
     - set the data file as containing continuous data
     -

   * - -input-data-cts-snp
     - set the .bed data file as containing SNP data to be treated as continuous data (0, 1, 2) in the network
     -

   * - -input-data-cts-snp2
     - set the data file as containing continuous SNP data (taking any continuous values)
     -

   * - -input-data-cts-missing-value x
     - set the value of missing data for continuous data to x
     -

   * - -input-data-discrete
     - set the data file as containing discrete data
     -

   * - -input-data-discrete-snp
     - set the .bed data file as containing SNP data to be treated as discrete data in the network
     -

   * - -input-data-discrete-snp2
     - set the data file as containing discrete SNP data (taking any number of discrete values)
     -

   * - -input-data-discrete-missing-value x
     - set the value of missing data for discrete data to x
     - NA

   * - -input-data-factor
     - set the data file as containing discrete data encoded using factor variables
     -

   * - -input-data-factor-snp
     - set the .bed file as containing SNP data to be treated as discrete factor data in the network
     -

   * - -input-data-factor-snp2
     - set the data file as containing discrete factor SNP data (taking any number of discrete values)
     -

   * - -input-data-factor-missing-value x
     - set the value of missing data for discrete factor data to x
     - NA

   * - -input-data-ids n
     - the number of ID columns in each data file
     - 2 

   * - -input-data-csv
     - set the data file as a comma separated file, .csv 
     -


.. _input-data-discrete:

Discrete data
-------------

Discrete data is automatically constrained to have no parent nodes that are continuous.


Discrete data is input by using the ``-input-data`` task, setting the data file and setting the data file to discrete. For example, the following


.. code-block:: none

  -input-data
  -input-data-file example-discrete.dat
  -input-data-discrete


could be used and the output would be something as follows:


.. code-block:: none

  BayesNetty: Bayesian Network software, v1.00
  --------------------------------------------------
  Copyright 2015-present Richard Howey, GNU General Public License, v3
  Institute of Genetic Medicine, Newcastle University

  Random seed: 1551695824
  --------------------------------------------------
  Task name: Task-1
  Loading data
  Discrete data file: example-discrete.dat
  Number of ID columns: 2
  Including the 1 and only variable in analysis
  Each variable has 1500 data entries
  Missing value: NA
  --------------------------------------------------

  Run time: less than one second


This parameter file can be found ``paras-input-discrete.txt`` in the examples, `example.zip <https://github.com/NewcastleRSE/BayesNetty/raw/refs/heads/main/docs/resources/example.zip>`_.


.. _input-data-cts:

Continuous data
---------------


Continuous data is input by using the ``-input-data`` task, setting the data file and setting the data file to continuous. For example, the following


.. code-block:: none

  -input-data
  -input-data-file example-cts.dat
  -input-data-cts


could be used and the output would be something as follows:

.. code-block:: none

  BayesNetty: Bayesian Network software, v1.00
  --------------------------------------------------
  Copyright 2015-present Richard Howey, GNU General Public License, v3
  Institute of Genetic Medicine, Newcastle University

  Random seed: 1551695897
  --------------------------------------------------
  Task name: Task-1
  Loading data
  Continuous data file: example-cts.dat
  Number of ID columns: 2
  Including (all) 2 variables in analysis
  Each variable has 1500 data entries
  Missing value: not set
  --------------------------------------------------

  Run time: less than one second



This parameter file can be found ``paras-input-cts.txt`` in the examples, `example.zip <https://github.com/NewcastleRSE/BayesNetty/raw/refs/heads/main/docs/resources/example.zip>`_.

.. _input-data-factor:

Factor data
-----------

Another way of handling discrete data is with the use of *factors*. Indicator variables are created, one for each different discrete category minus one.
These are treated as continuous explanatory variables in the linear regressions when they are parent nodes.
A restriction of using discrete data with factors is that they cannot be child nodes of other nodes. Input by using the ``-input-data`` task, setting the data file and setting the data file to factor.
For example, the following

.. code-block:: none

  -input-data
  -input-data-file example-discrete.dat
  -input-data-factor


could be used and the output would be something as follows:

.. code-block:: none

  BayesNetty: Bayesian Network software, v1.00
  --------------------------------------------------
  Copyright 2015-present Richard Howey, GNU General Public License, v3
  Institute of Genetic Medicine, Newcastle University

  Random seed: 1551695930
  --------------------------------------------------
  Task name: Task-1
  Loading data
  Discrete data file: example-discrete.dat
  Data treated as factors
  Number of ID columns: 2
  Including the 1 and only variable in analysis
  Each variable has 1500 data entries
  Missing value: NA
  --------------------------------------------------

  Run time: less than one second


This parameter file, ``paras-input-factor.txt``, can be found in the examples, `example.zip <https://github.com/NewcastleRSE/BayesNetty/raw/refs/heads/main/docs/resources/example.zip>`_.

.. _input-data-snp:

SNP data
--------

SNP data is automatically constrained to have no parent nodes.


SNP data may be input as a `binary <https://zzz.bwh.harvard.edu/plink/binary.shtml>`_ PLINK format pedigree file, a ``.bed`` file, see :cite:`purcell:etal:07` .
This requires that the corresponding ``.bim`` and ``.fam``, files are also available. A text PLINK pedigree file, ``.ped``, with corresponding map file, ``.map``, may be used to create a binary file using PLINK as follows:

.. code-block:: none

  plink --noweb --file mydata --make-bed --out myfile


This will create the binary pedigree file, ``myfile.bed``, map file, ``myfile.bim``, and family file, ``myfile.fam`` required.

The SNP data is input by using the ``-input-data`` task, setting the PLINK binary file and setting the data file to a SNP file in discrete mode or continuous mode. For example, in discrete mode, the following


.. code-block:: none

  -input-data
  -input-data-file example.bed
  -input-data-discrete-snp


could be used and the output would be something as follows:

.. code-block:: none

  BayesNetty: Bayesian Network software, v1.00
  --------------------------------------------------
  Copyright 2015-present Richard Howey, GNU General Public License, v3
  Institute of Genetic Medicine, Newcastle University

  Random seed: 1551695984
  --------------------------------------------------
  Task name: Task-1
  Loading data
  SNP binary data file: example.bed
  SNP data treated as discrete data
  Total number of SNPs: 2
  Total number of subjects: 1500
  Number of ID columns: 2
  Including (all) 2 variables in analysis
  Each variable has 1500 data entries
  --------------------------------------------------

  Run time: less than one second


This parameter file can be found ``paras-input-snp.txt`` in the examples, `example.zip <https://github.com/NewcastleRSE/BayesNetty/raw/refs/heads/main/docs/resources/example.zip>`_.


.. _input-data-missing:

Missing data
------------


Missing data is determined by any data matching the given missing value as defined by ``-input-data-discrete-missing-value`` and ``-input-data-cts-missing-value``
when inputting discrete and continuous data respectively (or ``-input-data-factor-missing-value`` when inputting factor data).
When continuous data has an invalid entry this will also be set to missing, for example a value of "NaN" will be set to missing since a numerical value is required.
Missing data for SNP data is given as defined by the PLINK `binary <https://zzz.bwh.harvard.edu/plink/binary.shtml>`_ pedigree format.
When there is missing data for a node for a certain individual then data for this certain individual is considered as missing for *every* node in the network.
Therefore the amount of missing data depends on which nodes are in the network.


Consider a network with 2 continuous nodes, with structure as given by network file ``example-network-missing1.dat`` and input using parameter file ``paras-input-missing1.txt``
as given in the file `example.zip <https://github.com/NewcastleRSE/BayesNetty/raw/refs/heads/main/docs/resources/example.zip>`_, then the output will be will look something as follows:


.. code-block:: none

  BayesNetty: Bayesian Network software, v1.00
  --------------------------------------------------
  Copyright 2015-present Richard Howey, GNU General Public License, v3
  Institute of Genetic Medicine, Newcastle University

  Random seed: 1551694585
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
  Task name: myNetwork
  Loading network
  Network file: example-network-missing1.dat
  Network type: bnlearn
  Network score type: BIC
  Total number of nodes: 2 (Discrete: 0 | Factor: 0 | Continuous: 2)
  Total number of edges: 1
  Network Structure: [express][pheno|express]
  Total data at each node: 1500
  Missing data at each node: 0
  --------------------------------------------------

  Run time: less than one second



As indicated in the network details there is no missing data. However, if the SNP node, ``rs1``, is added (network file ``example-network-missing2.dat``) then the following is given:


.. code-block:: none

  BayesNetty: Bayesian Network software, v1.00
  --------------------------------------------------
  Copyright 2015-present Richard Howey, GNU General Public License, v3
  Institute of Genetic Medicine, Newcastle University

  Random seed: 1551696539
  --------------------------------------------------
  Task name: Task-1
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
  Task name: Task-2
  Loading data
  Continuous data file: example-cts.dat
  Number of ID columns: 2
  Including (all) 2 variables in analysis
  Each variable has 1500 data entries
  Missing value: not set
  --------------------------------------------------
  --------------------------------------------------
  Task name: myNetwork
  Loading network
  Network file: example-network-missing2.dat
  Network type: bnlearn
  Network score type: BIC
  Total number of nodes: 3 (Discrete: 1 | Factor: 0 | Continuous: 2)
  Total number of edges: 2
  Network Structure: [rs1][express|rs1][pheno|express]
  Total data at each node: 1497
  Missing data at each node: 3
  --------------------------------------------------

  Run time: less than one second



This example is given in network file ``example-network-missing2.dat``and parameter file ``paras-input-missing2.txt``.
The amount of missing data for the network is now 3, indicating that 3 individuals have missing SNP data for ``rs1``.
Adding in another SNP node, ``rs2`` (network file ``example-network-missing3.dat``), results in the following:


.. code-block:: none

  BayesNetty: Bayesian Network software, v1.00
  --------------------------------------------------
  Copyright 2015-present Richard Howey, GNU General Public License, v3
  Institute of Genetic Medicine, Newcastle University

  Random seed: 1551696644
  --------------------------------------------------
  Task name: Task-1
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
  Task name: Task-2
  Loading data
  Continuous data file: example-cts.dat
  Number of ID columns: 2
  Including (all) 2 variables in analysis
  Each variable has 1500 data entries
  Missing value: not set
  --------------------------------------------------
  --------------------------------------------------
  Task name: myNetwork
  Loading network
  Network file: example-network-missing3.dat
  Network type: bnlearn
  Network score type: BIC
  Total number of nodes: 4 (Discrete: 2 | Factor: 0 | Continuous: 2)
  Total number of edges: 3
  Network Structure: [rs1][rs2][express|rs1:rs2][pheno|express]
  Total data at each node: 1495
  Missing data at each node: 5
  --------------------------------------------------

  Run time: less than one second



Similarly, this example is given in network file ``example-network-missing3.dat`` and parameter file ``paras-input-missing3.txt``.
Here we see that the amount of missing data in the network has increased due to missing data for SNP node ``rs2``.
This node also has missing data for 3 individuals, with the result that the total amount of missing data for each node is 5.


.. _input-data-ids: 

Data IDs
--------


By default the first two columns of a data file should be IDs and match those in any other data files, and be in the same order (although the ID names in the header do not need to match).
The number of ID columns can be changed using the ``-input-data-ids`` option, and may be set to zero.
If the data contains SNP data in a PLINK binary pedigree file, ``.bed``, then the number of ID columns must be set to 2.
If the data is a binary pedigree file, ``file.bed``, then the family and individual IDs in the file ``file.fam`` must match the IDs in any other data files, and all SNPs may be used as network nodes.
The IDs in different files are checked to be the same including the order, if not BayesNetty will report an error.
If there are zero IDs then the individuals are assumed to be in the same order in each file and are only checked to have the same number of individuals.


.. _input-data-example: 

Example
-------

See the above sections for examples of inputting data.


 