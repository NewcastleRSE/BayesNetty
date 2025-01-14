.. _output-network:

Output network
==============

A network may be output using the `-output-network` task. The network may be output in one of three different formats.
The `-output-network-igraph-file-prefix` option can also be used to output igraph files, see :ref:`plot-network` for more details.


The output network task can also be used to output node data using the `-output-network-node-data-file-prefix` option and optionally the `-output-network-node-data-bed-file` option. 

.. _output-network-options:

Options
-------

The options are as follows:

.. list-table:: 
    :header-rows: 1

    * - Option
      - Description
      - Default

    * - -output-network-name name
      - label the task with a name
      - Task-n

    * - 
      - 
      -

    * - 
      - 
      -

    * - 
      - 
      -

    * - 
      - 
      -

    * - 
      - 
      -

    * - 
      - 
      -

    * - 
      - 
      -

    * - 
      - 
      -



*tr*
   &  & 
*/tr*

*tr*
  -output-network-network-name network & output this network & previous network (or the default model given by a node for each data variable and no edges if there is no previous network)
*/tr*

*tr*
  -output-network-file network.dat & output the network in a format where the nodes and then the edges are listed &
*/tr*

*tr*
  -output-network-file2 network2.dat & output the network in this style of format: *code* [a][b|a][c|a:b] */code* &
*/tr*

*tr*
-output-network-equivalent-networks-file equiv-networks.dat & output a list of equivalent networks to file equiv-networks.dat &
*/tr*

*tr*
  -output-network-igraph-file-prefix mygraph & output igraph format files consisting of mygraph-nodes.dat, mygraph-edges.dat and R code mygraph-plot.R &
*/tr*

*tr*
-output-network-node-data-file-prefix mydata & output network discrete data to mydata-discrete.dat and continuous data to mydata-cts.dat &
*/tr*

*tr*
-output-network-node-data-bed-file & output SNP data to files mydata.bed, mydata.bim and mydata.fam and not to files mydata-discrete.dat and mydata-cts.dat &
*/tr*

*tr*
-output-network-node-data-start-indiv a & start output of data from individual a & 1
*/tr*

*tr*
-output-network-node-data-end-indiv b & stop output of data at individual b & last individual
*/tr*

*tr*
-output-network-node-data-job i t & only output individuals for subset i from a total of t subsets &
*/tr*

*/table*

*/subsection*

*************

*subsection*

*subsection-name* output-network-example */subsection-name*

*subsection-title* Example */subsection-title*

* 
The following is an example parameter file to output a network.
*

*codeexample*
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

#output the fitted network
-output-network
-output-network-file fittedNetwork.dat
*/codeexample*

*
This parameter file, *code* paras-output-network.txt */code*, can be found in *html* <a href="example.zip">example.zip</a>&nbsp; */html* *tex* example.zip */tex* and can be used as follows:
*

*codeexample*
./bayesnetty paras-output-network.txt
*/codeexample*

*
Which should produce output that looks like something as follows:
*

*codeexample*
BayesNetty: Bayesian Network software, v1.00
--------------------------------------------------
Copyright 2015-present Richard Howey, GNU General Public License, v3
Institute of Genetic Medicine, Newcastle University

Random seed: 1551716789
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
--------------------------------------------------
Task name: Task-5
Outputting network
Network: defaultNetwork
Network Structure: [mood][rs1][rs2][express|rs1:rs2][pheno|express:mood]
Network output to file: fittedNetwork.dat
--------------------------------------------------

Run time: less than one second
*/codeexample*

*
The data is loaded, a search is performed and then the network is output to a file.
*


*/subsection*

*************************

*/section* 
