.. _estimate-impute: 

Estimate Imputation Benefit
===========================


For a given data set with missing data it is natural to wonder how much benefit imputation brings.
We include an option in BayesNetty to attempt to compare different methods of fitting a best fit network to this data set.
We use estimates of the recall and precision to compare the methods.
The recall is the percentage of edges that were recovered from the simulation model and the precision is the percentage of edges in the fitted model that are in the simulation model.
For an edge to be correct it must be in the correct direction. However, if the simulating network has equivalent networks such that some edges may be in either direction, then these edges are considered correct if they are in any direction.

The method follows these steps:


#. An initial network is fitted using imputation.
#. Data is simulated using this network for the same number of individuals in the original data set.
#. A best fit network is found for the full data set.
#. The simulated data set has values set to missing as in the original data set.
#. Best fit networks are found for this data set using: (i) a reduced data set with only complete data; (ii) data imputation; (iii) data imputation with complete training data.
#. Recall and precision are calculated for the 4 different best fit networks against the simulation network.


This method can be repeated a number times as is computationally feasible to take average recall and precision estimates to account for variability in the simulated data.
The recall and precision using the full data set gives an estimate of an upper limit of what may feasibly be achieved using data imputation.
Comparing the recall and precision of the reduced data set with imputation gives an estimate of the increased benefit of using imputation.
Comparing the two imputation methods should show when it is appropriate to use the variant imputation method.


A major drawback of this estimation method is the obvious fact that we do not know the "true" network structure of the data,
we therefore use an estimated network to simulate the data and hope this is sufficiently close for the results to be useful.
In general, we have found that the benefits of imputation are often understated as the simulation network tends to be set without some of the weaker edges that cannot always be detected
(when using data sets where we do actually know the *q* true */q* network).
Even if we cannot be too sure of the exact gain in benefit of imputation this BayesNetty estimation method can give clear confidence of a benefit when there are large differences (and if the variant method using complete training data performs any better).

 
.. _estimate-impute-options:

Options
-------

The options are as follows:

.. list-table:: 
    :header-rows: 1

    * - Option
      - Description
      - Default

    * - -impute-estimate-recall-precision
      - do a task to estimate recall and precision before and after data imputation
      - 

    * - -impute-estimate-recall-precision name
      - label the task with a name
      - Task-n

    * - -impute-estimate-recall-precision-network-name network
      - set the name of the initial network when estimating recall and precision
      - previous network (or the default model given by a node for each data variable and no edges if there is no previous network)

    * - -impute-estimate-recall-precision-random-restarts n
      - for each network fit do another n searches starting from a random network
      - 

    * - -impute-estimate-recall-precision-jitter-restarts m
      - for each network fit after the initial search and every random restart search do another m searches jittered from the recently found network
      - 0

    * - -impute-estimate-recall-precision-skip-imputation
      - do not estimate recall and precision for imputed data
      - 

    * - -impute-estimate-recall-precision-iterations i
      - estimate the recall and precision i times and take the average
      - 1


.. _estimate-impute-example: 

Example
-------

An example of estimating the recall and precision is contained in the parameter file `paras-example-estimate-recall-precision.txt`,
which can be found in `example.zip <https://github.com/NewcastleRSE/BayesNetty/raw/refs/heads/main/docs/resources/example.zip>`_.
For simplicity the example is chosen to be a discrete network and this approach can be used for any kind of data.
The network is the `child` network from the bnlearn repository :cite:`bnlearn`.


.. code-block:: none

    #input example data to estimate recall and precision from
    -input-data
    -input-data-file data-example-est-recall-precision.dat
    -input-data-ids 0
    -input-data-discrete

    #set up network with no edges
    -input-network
    -input-network-empty

    #estimate recall and precision for this data set
    -impute-estimate-recall-precision
    -impute-estimate-recall-precision-iterations 10
    -impute-estimate-recall-precision-random-restarts 2
    -impute-estimate-recall-precision-jitter-restarts 2


This can be executed as usual


.. code-block:: none

    ./bayesnetty paras-example-estimate-recall-precision.txt


and will output something as follows


.. code-block:: none

    BayesNetty: Bayesian Network software, v1.1
    --------------------------------------------------
    Copyright 2015-present Richard Howey, GNU General Public License, v3
    Institute of Genetic Medicine, Newcastle University

    Random seed: 1605624192
    --------------------------------------------------
    Task name: Task-1
    Loading data
    Discrete data file: data-example-est-recall-precision.dat
    Number of ID columns: 0
    Including (all) 20 variables in analysis
    Each variable has 500 data entries
    Missing value: NA
    --------------------------------------------------
    --------------------------------------------------
    Task name: Task-2
    Loading network
    Network set with no edges
    Network type: bnlearn
    Network score type: BIC
    Total number of nodes: 20 (Discrete: 20 | Factor: 0 | Continuous: 0)
    Total number of edges: 0
    Network Structure: [Age][BirthAsphyxia][CO2Report][CO2][CardiacMixing][ChestXray][Disease][DuctFlow]
    [GruntingReport][Grunting][HypDistrib][HypoxiaInO2][LVH][LVHreport][LowerBodyO2][LungFlow][LungParench][RUQO2][Sick][XrayReport]
    Total data at each node: 54
    Missing data at each node: 446
    --------------------------------------------------
    --------------------------------------------------
    Task name: Task-3
    Estimating the recall and precision when imputing network data
    Network: Task-2
    Number of iterations: 10
    Random restarts: 2
    Random jitter restarts: 2
    Minimum percentage of non-missing edges (or singleton nodes) required to impute individual: 0
    Individuals with data: 54
    Individuals with missing data: 446

    Recall: the percentage of edges found from the original true network
    Precision: the percentage of edges in the network that are also in the original true network

                                      Recall     Precision
    No imputation                      40.91      61.66
    Imputation                         78.95      89.03
    Imputation (complete training)     67.75      80.04
    Full data (upper limit)            90.21      95.02

    --------------------------------------------------

    Run time: 1 hour, 39 minutes and 23 seconds


From the example output we can see that with no imputation the recall is estimated to be 40.91 percent and the precision estimated to be is 61.66 percent,
but if the full data were available it would be 90.21 and 95.02 respectively.
Using our imputation method the estimated recall and precision is 78.95 and 89.03 respectively, which is quite a large increase.
Our variant imputation method with complete training data also increases the recall and precision by quite a lot.


Note that the estimation is stochastic due to the stochastic nature of the imputation method, and to a lesser extent the stochastic nature of finding a best fit model,
and so rerunning the analyses may results in slightly different estimates.   
