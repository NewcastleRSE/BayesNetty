.. _using-bayesnetty:

Using BayesNetty
================


BayesNetty is executed as follows:

.. code-block:: none

    ./bayesnetty paras.txt

where ``paras.txt`` is a parameter file as described in the following sections.

.. _basic-options:

Basic Options
-------------

The basic options for BayesNetty are as follows (typing ``./bayesnetty`` with no options will output the available options):

================  =============================  =======================
Option            Description                    Default               
================  =============================  =======================
-log file.log     log file of screen output      bayesnetty.log        
-so               suppress output to screen                           
-seed number      random number generator seed   set by execution time 
================  =============================  =======================

Random Seed
^^^^^^^^^^^

The random seed option ``-seed`` may be used to ensure exactly the same output for testing and reproducibility purposes.
If you do this it is important to use the same number of processes if also using the parallel version of BayesNetty. 

Parameter file
--------------


There are many different things that BayesNetty can do, each of these different things is referred to as a "task".
The parameter file defines which tasks BayesNetty will perform and the order in which they are executed. With the exception of the :ref:`basic options<basic-options>`,
all options in the parameter file define tasks. For example, a task to input some continuous data may be given as follows: 

.. code-block:: none

    -input-data
    -input-data-name myGreatData
    -input-data-file example-cts.dat
    -input-data-cts
    -input-data-ids 1


There are a few basic rules for parameter files: 

1. Each option must be written on a separate line.

2. Each line that does not start with a dash, "-", will be ignored, thus allowing comments to be written.

3. The task must first be declared and then followed by any options for the task. 

4. An option for a task is always written by first writing the name of the task. For example, the task option to give the name of an input data file is given by ``-input-data-file``, which begins with ``-input-data``, the name of the task to input data.

5. A task, ``XXX``, may always be given a task name with the option ``XXX-name``. The task may then be referenced by other tasks (if permitted). This may be useful if there is more than one network.

6. Tasks are executed in order, so any tasks that depend on other tasks must be ordered accordingly.

7. Any tasks that require a network will be default use the previously defined network. Therefore, if there is only one network it is not necessary to name or reference it. By default tasks are name "Task-n", where n is the number of the task.


The following is an extract from an example parameter file where a network is referenced by a task to calculate the network score:

.. code-block:: none

    ...

    # This is my comment
    -input-network
    -input-network-name myNetwork
    -input-network-file network-model.dat

    -calc-network-score
    -calc-network-score-network-name myNetwork


The parameter file could be thought of as in an ``R`` programming style, such that the above would look as follows:

.. code-block:: R

    ...

    # This is my comment
    myNetwork<-input.network(file = "network-model.dat")

    calc.network.score(network.name = myNetwork)


However, as BayesNetty is not an ``R`` package (or a programming language), the parameter file uses an unambiguous, longhand, and easy to parse style of syntax.


The options for all the different tasks may be found in the different task sections of the documentation.




*/subsection*

*************

*subsection*

*subsection-name* simple-example */subsection-name*

*subsection-title* Simple Example */subsection-title*

*
Example data and parameter files can be found in the file *html* <a href="example.zip">example.zip</a> */html* *tex* example.zip */tex* . The example parameter file, ``paras-example.txt``, can be used to perform a simple analysis by typing
*

.. code-block:: none
 ./bayesnetty paras-example.txt


*
The following shows the ``paras-example.txt`` file 
*

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
-search-models-file search-example.dat


*
The parameter file instructs BayesNetty to perform 4 tasks: (i) load continuous data from file ``example-cts.dat``; (ii) load discrete data from file ``example-discrete.dat``; (iii) load SNP data to be treated as discrete data from file ``example.bed``; and finally (iv) search the network models. The screen output, which is also saved to a log file, will look something as follows:
*

.. code-block:: none
BayesNetty: Bayesian Network software, v1.00
--------------------------------------------------
Copyright 2015-present Richard Howey, GNU General Public License, v3
Institute of Genetic Medicine, Newcastle University

Random seed: 1551700145
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
Network search output to file: search-example.dat
--------------------------------------------------

Run time: less than one second



*/subsection*

********************************

*subsection*

*subsection-name* command-line */subsection-name*

*subsection-title* Command-line Options */subsection-title*

*
It is also possible to add options on the command line to modify or add to the options in the parameter file. For example
*

.. code-block:: none
./bayesnetty paras-example.txt -seed 1 -log seed-1-results.log


