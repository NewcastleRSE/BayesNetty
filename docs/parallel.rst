.. _parallel: 

Parallel BayesNetty
===================


It is possible to run BayesNetty using Open MPI (*cite* openmpi */cite*), which is an open source Message Passing Interface (MPI) (*cite* mpi */cite*) implementation designed for parallel programming. 


The parallel version of Bayesnetty speeds up the search through network space for the best network by simultaneously evaluating different networks. This is particularly useful for large networks and any type of analysis that depends on network searches, such as network averaging and imputing data.


A much faster way to calculate an average network in parallel is given in the :ref:`Average Network<average-network-parallel>` section.

A much faster way to impute network data in parallel is given in the :ref:`Impute Data<impute-parallel-example>` section. 


After installing Open MPI on your system if it is not already installed, see *cite* openmpi */cite*, a parallel version of BayesNetty needs to be compiled. This can be done by firstly uncommenting a few lines in the ``main.h`` file. So that the following:

.. code-block:: none

    // Comment out if not using Open MPI for parallel processing
    //#ifndef USING_OPEN_MPI
    //#define USING_OPEN_MPI
    //#endif //OPEN_MPI
    
becomes

.. code-block:: none

    // Comment out if not using Open MPI for parallel processing
    #ifndef USING_OPEN_MPI
    #define USING_OPEN_MPI
    #endif //OPEN_MPI
    
then compile the parallel code as follows:


.. code-block:: none

    mpicxx -O3 -o pbayesnetty *.cpp

The code can then be ran using how many processes that you wish, for example to run with 12 processes use

.. code-block:: none

    mpirun -n 12 ./pbayesnetty paras-example.txt

For such a trivial example the code will not run any quicker, and in fact for very small networks one may find that analyses take longer.
There is some overhead in using the MPI libraries, so if trivial networks are used there may be no speed up.

Even for large networks the optimal number of processes to perform the analysis as quick as possible may not be as many processes as you can use. As there is an overhead for processes the best amount to use may be a lot, but not too many...
The best amount will vary depending on the analysis, the data and the computing system that you are using, so some trial and error may be needed.


The output will show the number of processes as well as the random seed. If you wish to reproduce exactly the same results both of these need to be set to the same value.
The seed is set with the ``-seed`` option.


.. code-block:: none

    BayesNetty: Bayesian Network software, v1.00
    --------------------------------------------------
    Copyright 2015-present Richard Howey, GNU General Public License, v3
    Institute of Genetic Medicine, Newcastle University

    Number of processes: 12
    Random seed: 1541430503
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

    ...



.. _compile-parallel-code: 

Compilation Scripts
-------------------


Scripts to compile Bayesnetty as either parallel or non-parallel while automatically uncommenting or commenting the code as appropriate are given below.


Script to compile code in parallel:

.. code-block:: none

    sed -i s://#ifndef\ USING_OPEN_MPI:#ifndef\ USING_OPEN_MPI:g main.h
    sed -i s://#define\ USING_OPEN_MPI:#define\ USING_OPEN_MPI:g main.h
    sed -i s://#endif\ //:#endif\ //:g main.h

    mpicxx -O3 -o pbayesnetty *.cpp


Script to compile code in non-parallel:

.. code-block:: none

    sed -i s://#ifndef\ USING_OPEN_MPI:#ifndef\ USING_OPEN_MPI:g main.h
    sed -i s://#define\ USING_OPEN_MPI:#define\ USING_OPEN_MPI:g main.h
    sed -i s://#endif\ //:#endif\ //:g main.h

    sed -i s:#ifndef\ USING_OPEN_MPI://#ifndef\ USING_OPEN_MPI:g main.h
    sed -i s:#define\ USING_OPEN_MPI://#define\ USING_OPEN_MPI:g main.h
    sed -i s:#endif\ //://#endif\ //:g main.h

    g++ -O3 *.cpp -o bayesnetty

