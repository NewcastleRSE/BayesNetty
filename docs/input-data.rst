.. _input-data:

Input data
==========

All data must be input using the ``-input-data`` task.

.. _input-data-options:

Options
-------

The options are as follows:


=========================================  =========================================  =========================================
Option                                     Description                    Default               
=========================================  =========================================  =========================================
-input-data                                do a task to input data  
-input-data-name name                      label the task with a name                 Task-n
-impute-network-data-network-name network  the name of the network to                 previous network (or the default model if 
                                           impute data for                            there is no previous network) given by a 
                                                                                      node for each data variable and no edges
=========================================  =========================================  =========================================

.. list-table:: 
   :widths: 50 50
   :header-rows: 1

   * - Column 1
     - Column 2
   * - Short
     - A much longer entry
   * - Another
     - Example

*tr*
  -input-data-include-file nodes.dat & a list of nodes/variables from the data file to be included in the network. Only the nodes in this list will be used in any analysis. &
*/tr*

*tr*
  -input-data-exclude-file nodes.dat & a list of nodes/variables to be excluded from the network &
*/tr*

*tr*
  -input-data-cts & set the data file as containing continuous data &
*/tr*

*tr*
  -input-data-cts-snp & set the .bed data file as containing SNP data to be treated as continuous data (0, 1, 2) in the network &
*/tr*

*tr*
  -input-data-cts-snp2 & set the data file as containing continuous SNP data (taking any continuous values) &
*/tr*

*tr*
  -input-data-cts-missing-value x & set the value of missing data for continuous data to x &
*/tr*

*tr*
  -input-data-discrete & set the data file as containing discrete data &
*/tr*

*tr*
  -input-data-discrete-snp & set the .bed data file as containing SNP data to be treated as discrete data in the network &
*/tr*

*tr*
  -input-data-discrete-snp2 & set the data file as containing discrete SNP data (taking any number of discrete values) &
*/tr*

*tr*
  -input-data-discrete-missing-value x & set the value of missing data for discrete data to x & NA
*/tr*

*tr*
-input-data-factor & set the data file as containing discrete data encoded using factor variables &
*/tr*

*tr*
-input-data-factor-snp & set the .bed file as containing SNP data to be treated as discrete factor data in the network &
*/tr*

*tr*
-input-data-factor-snp2 & set the data file as containing discrete factor SNP data (taking any number of discrete values) &
*/tr*

*tr*
-input-data-factor-missing-value x & set the value of missing data for discrete factor data to x & NA
*/tr*

*tr*
  -input-data-ids n & the number of ID columns in each data file & 2 
*/tr*

*tr*
-input-data-csv & set the data file as a comma separated file, .csv &
*/tr*

*/table*

*/subsection*

**********************

*subsection*

*subsection-name* input-data-discrete */subsection-name*

*subsection-title* Discrete data */subsection-title*

*
Discrete data is automatically constrained to have no parent nodes that are continuous.
*

*
Discrete data is input by using the ``-input-data`` task, setting the data file and setting the data file to discrete. For example, the following
*

.. code-block:: none
-input-data
-input-data-file example-discrete.dat
-input-data-discrete


*
could be used and the output would be something as follows:
*

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


*
This parameter file can be found ``paras-input-discrete.txt``in the examples, *html* <a href="example.zip">example.zip</a> */html* *tex* example.zip */tex* .
*

*/subsection*

**********************

*subsection*

*subsection-name* input-data-cts */subsection-name*

*subsection-title* Continuous data */subsection-title*

*
Continuous data is input by using the ``-input-data``task, setting the data file and setting the data file to continuous. For example, the following
*

.. code-block:: none
-input-data
-input-data-file example-cts.dat
-input-data-cts


*
could be used and the output would be something as follows:
*

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


*
This parameter file can be found ``paras-input-cts.txt``in the examples, *html* <a href="example.zip">example.zip</a> */html* *tex* example.zip */tex* .
*

*/subsection*

**********************
**********************

*subsection*

*subsection-name* input-data-factor */subsection-name*

*subsection-title* Factor data */subsection-title*

*
Another way of handling discrete data is with the use of *i* factors */i*. Indicator variables are created, one for each different discrete category minus one. These are treated as continuous explanatory variables in the linear regressions when they are parent nodes. A restriction of using discrete data with factors is that they cannot be child nodes of other nodes. Input by using the ``-input-data``task, setting the data file and setting the data file to factor. For example, the following
*

.. code-block:: none
-input-data
-input-data-file example-discrete.dat
-input-data-factor


*
could be used and the output would be something as follows:
*

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


*
This parameter file can be found ``paras-input-factor.txt``in the examples, *html* <a href="example.zip">example.zip</a> */html* *tex* example.zip */tex* .
*

*/subsection*

**********************

*subsection*

*subsection-name* input-data-snp */subsection-name*

*subsection-title* SNP data */subsection-title*

*
SNP data is automatically constrained to have no parent nodes.
*

*
SNP data may be input as a *html* <a href="http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml">binary</a>&nbsp; */html* *tex* binary */tex* PLINK format pedigree file, a ``.bed``file, see *cite* purcell:etal:07 */cite*. This requires that the corresponding ``.bim``and ``.fam */code*, files are also available. A text PLINK pedigree file, ``.ped */code*, with corresponding map file, ``.map */code*, may be used to create a binary file using PLINK as follows:

.. code-block:: none
plink --noweb --file mydata --make-bed --out myfile


* This will create the binary pedigree file, ``myfile.bed */code*, map file, ``myfile.bim */code*, and family file, ``myfile.fam``required. *

*
The SNP data is input by using the ``-input-data``task, setting the PLINK binary file and setting the data file to a SNP file in discrete mode or continuous mode. For example, in discrete mode, the following
*

.. code-block:: none
-input-data
-input-data-file example.bed
-input-data-discrete-snp


*
could be used and the output would be something as follows:
*

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


*
This parameter file can be found ``paras-input-snp.txt``in the examples, *html* <a href="example.zip">example.zip</a> */html* *tex* example.zip */tex* .
*


*/subsection*

**********************


*subsection*

*subsection-name* input-data-missing */subsection-name*

*subsection-title* Missing data */subsection-title*

*
Missing data is determined by any data matching the given missing value as defined by ``-input-data-discrete-missing-value``and ``-input-data-cts-missing-value``when inputting discrete and continuous data respectively (or ``-input-data-factor-missing-value``when inputting factor data). When continuous data has an invalid entry this will also be set to missing, for example a value of *q* NaN */q* will be set to missing since a numerical value is required. Missing data for SNP data is given as defined by the PLINK *html* <a href="http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml">binary</a>&nbsp; */html* *tex* binary */tex* pedigree format. When there is missing data for a node for a certain individual then data for this certain individual is considered as missing for *i* every */i* node in the network. Therefore the amount of missing data depends on which nodes are in the network.
*

*
Consider a network with 2 continuous nodes, with structure as given by network file ``example-network-missing1.dat``and input using parameter file ``paras-input-missing1.txt``as given in *html* <a href="example.zip">example.zip</a>,&nbsp; */html* *tex* example.zip, */tex* then the output will be will look something as follows:
*

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


*
As indicated in the network details there is no missing data. However, if the SNP node, ``rs1 */code*, is added (network file ``example-network-missing2.dat */code*) then the following is given:
*


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


*
This example is given in network file ``example-network-missing2.dat``and parameter file ``paras-input-missing2.txt */code*. The amount of missing data for the network is now 3, indicating that 3 individuals have missing SNP data for ``rs1 */code*. Adding in another SNP node, ``rs2``(network file ``example-network-missing3.dat */code*), results in the following:
*

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


*
Similarly, this example is given in network file ``example-network-missing3.dat``and parameter file ``paras-input-missing3.txt */code*. Here we see that the amount of missing data in the network has increased due to missing data for SNP node ``rs2 */code*. This node also has missing data for 3 individuals, with the result that the total amount of missing data for each node is 5.
*

*/subsection*
**********************

*subsection*

*subsection-name* input-data-ids */subsection-name*

*subsection-title* Data IDs */subsection-title*

*
By default the first two columns of a data file should be IDs and match those in any other data files, and be in the same order (although the ID names in the header do not need to match). The number of ID columns can be changed using the ``-input-data-ids``option, and may be set to zero. If the data contains SNP data in a PLINK binary pedigree file, ``.bed */code*, then the number of ID columns must be set to 2. If the data is a binary pedigree file, ``file.bed */code*, then the family and individual IDs in the file ``file.fam``must match the IDs in any other data files, and all SNPs may be used as network nodes. The IDs in different files are checked to be the same including the order, if not BayesNetty will report an error. If there are zero IDs then the individuals are assumed to be in the same order in each file and are only checked to have the same number of individuals.
*


*/subsection*
**********************
*subsection*

*subsection-name* input-data-example */subsection-name*

*subsection-title* Example */subsection-title*

* See the above sections for examples of inputting data. *


*/subsection*

*************

*/section* 

********************************

 ********************************

*section*

*section-name* input-network */section-name*

*section-title* Input network */section-title*

*
A network may be specified using the ``-input-network``task. The network may be used as a starting point for analyses, such as searches, or to perform an analysis on this network. Only nodes in input files will be used in the network so that a subset of the data may be specified.
*

*
Any network *b* constraints */b* must be set using the ``-input-network``option. These constraints then belong to the network and will be used in any subsquent analysis, including searches, calculating average networks etc. 
*

*
If no network is specified then a network with no edges and a node for every data variable (as given by the input data) will be created and named *q* defaultNetwork */q*.
*

**********************

*subsection*

*subsection-name* input-network-options */subsection-name*

*subsection-title* Options */subsection-title*

* The options are as follows: *

*tablelopt* *tr* Option & Description & Default */tr*

*tr*
 -input-network  & do a task to input a network &
*/tr*

*tr*
  -input-network-name name & label the task and network with a name & Task-n
*/tr*

*tr*
  -input-network-type t & the type of Bayesian network, choose between bnlearn or deal &  bnlearn 
*/tr*

*tr*
  -input-network-file network.dat & input the network in a format where the nodes and then the edges are listed &
*/tr*

*tr*
  -input-network-file2 network2.dat & input the network in this style of format: ``[a][b|a][c|a:b]``&
*/tr*

*tr*
  -input-network-igraph-file-prefix mygraph & input the network from igraph format files consisting of mygraph-nodes.dat and mygraph-edges.dat &
*/tr*

*tr*
  -input-network-empty & set the network to one with no edges and one node for every data variable. An input network file is not required if this option is used &
*/tr*

*tr*
  -input-network-whitelist-file whitelist.dat & a list of edges that must be included in any network &
*/tr*

*tr*
  -input-network-blacklist-file blacklist.dat & a list of edges that must *i* not */i* be included in any network &
*/tr*

*tr*
  -input-network-blacklist-edge-type dataName1 dataName2 & edge types that may *i* not */i* be included in any network. The collection of nodes are given by the data input name, and so the data types must be given in different files &
*/tr*

*tr*  
  -input-network-no-parents-node nodeX        &       nodeX must not have any parents (except for white edges) &
*/tr*

*tr*
	-input-network-no-children-node nodeY     &       nodeY must not have any children (except for white edges) &
*/tr*

*tr*
	-input-network-prob-edge node1 node2 prob  & set the prior probability of edge direction of node1 to node2 as prob &
*/tr*

*tr*
   -input-network-prob-edge-type nodeType1 nodeType2 prob & set the prior probabilities of edge direction of nodeType1 to nodeType2 as prob &
*/tr*


*tr*
  -input-network-imaginary-sample-size i & for deal networks this sets the imaginary sample size & 10
*/tr*

*tr*
  -input-network-score score & for a bnlearn network choose between loglike, AIC or BIC & BIC
*/tr*

*comment*
*tr*
-input-network-score-fix & fix for likelihood calculation when discrete data (or combinations) has too few points for linear regression, choose between none, average or skip & none
*/tr*
*/comment*

*/table*

*/subsection*

*************

*subsection*

*subsection-name* input-network-black */subsection-name*

*subsection-title* Black lists */subsection-title*

* A black list can be given using the ``-input-network-blacklist-file``option to define a list of edges that must not be included in any network. The text file should be formatted as follows: *

.. code-block:: none
node1 node2
node2 node1
node1 node3


*
such that the two nodes of each blacklisted edge are on one line. The nodes are ordered so the first line states that the edge node1 to node2 is not permitted. The next line states that the edge in the reverse direction is also not permitted.
*

*
Any searches will ignore these blacklisted edges and attempting to use a network with a blacklisted edge will result in the edge being removed.
*

*
Edges between different types of nodes may also be blacklisted. This can be done using the ``-input-network-blacklist-edge-type``option. It can be used as follows:
*

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


*
Firstly the different node types must be loaded separately and given names using the ``-input-data-name``option. Then, when initially loading a network, the ``-input-network-blacklist-edge-type``can be used to forbid any edge from one data set to another data set (or the same data if desired). In the above example the network may not have any edge that goes from a horse to a whip, that is, a whip node may not have a horse node as a parent. In any search that is performed these edges will not be considered. 
*



*/subsection*

*************

*subsection*

*subsection-name* input-network-white */subsection-name*

*subsection-title* White lists */subsection-title*

* A white list can be given using the ``-input-network-whitelist-file``option to define a list of edges that must be included in any network. The text file should be formatted as follows: *

.. code-block:: none
node1 node3
node1 node2
node2 node1


*
such that the two nodes of each whitelisted edge are on one line. The nodes are ordered so the first line states that the edge node1 to node3 must be included. If both directions are included between two nodes then the edge must be included but may be in any direction.
*

*
If the whitelist and blacklist contradict one another then an error will be given.
*


*/subsection*

*************

*subsection*

*subsection-name* input-network-soft-con */subsection-name*

*subsection-title* Soft Constraints */subsection-title*

*
Soft constraints provide a way that the direction of an edge may be influenced but not with certainty, unlike blacklisted edges or whitelisted edges as described above. An example parameter file setting a soft constraint, such that the prior probability of variable ``express``to variable ``pheno``is believed to be 0.8 is shown below.
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

#input the example network in format 1
-input-network
-input-network-name myNetwork
-input-network-file example-network-format1.dat
-input-network-prob-edge express pheno 0.8

#search network models with the soft constraint
-search-models



*
This parameter file, ``paras-soft-constraints.txt */code*, can be found in *html* <a href="example.zip">example.zip</a> */html* *tex* example.zip */tex*.
*


*
Any searches will use this prior probability.
*

*
If you wish to blacklist or whitelist an edge you should use those options rather than setting the prior probability to 0 or 1 for the sake of computational efficiency.
*

*/subsection*

**********************
*subsection*

*subsection-name* input-network-formats */subsection-name*

*subsection-title* Network formats */subsection-title*

* The network may be defined using one of 3 different formats. *


*subsubsection*

*subsubsection-name* input-network-formats-format1 */subsubsection-name*

*subsubsection-title* Network file format 1 */subsubsection-title*

* The first format is given by using the ``-input-network-file``option and the network text file should be formatted as follows: *

.. code-block:: none
node1
node2
node3
node2 node1
node3 node1


*
where the nodes are listed first followed by the directed edges. In the above example there are 3 nodes and 2 edges, which are node2 to node1 and node3 to node1.
*

*/subsubsection*


*subsubsection*

*subsubsection-name* input-network-formats-format2 */subsubsection-name*

*subsubsection-title* Network file format 2 */subsubsection-title*

* The second format is given by using the ``-input-network-file2``option and the network text file should be formatted as follows: *

.. code-block:: none
[node2][node3][node1|node2:node3]


*
where the nodes are listed in order of dependency. The independent nodes node2 and node3 are list first followed by node1 which is a child node of both node2 and node3. This is the format that is typically output for searches and such like.
*


*/subsubsection*


*subsubsection*

*subsubsection-name* input-network-formats-format3 */subsubsection-name*

*subsubsection-title* Network file format 3 */subsubsection-title*

* The third format is given by using the ``-input-network-igraph-file-prefix``option using the files that were output to draw the network in ``R */code*, see *ref* igraph */ref*. There will be one file for the nodes and one for the edges, for example ``myNetwork-nodes.dat``and ``myNetwork-edges.dat``respectively. The node file will look something as follows: *

.. code-block:: none
id name type fileno
1 node1 c 1
2 node2 c 1
3 node3 c 1


*
and the edges file will look like something as follows:
*

.. code-block:: none
from to chisq
2 1 6860.83
3 1 5709.51


*/subsubsection*

*/subsection*

**********************

*************

*subsection*

*subsection-name* input-network-example */subsection-name*

*subsection-title* Example */subsection-title*

* 
The following is an example parameter file to input a network.
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

#input the example network in format 1
-input-network
-input-network-name myNetwork
-input-network-file example-network-format1.dat


*
This parameter file, ``paras-input-network.txt */code*, can be found in *html* <a href="example.zip">example.zip</a>&nbsp; */html* *tex* example.zip */tex* and can be used as follows:
*

.. code-block:: none
./bayesnetty paras-input-network.txt


*
Which should produce output that looks like something as follows:
*

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


*
The data is loaded and then the network is loaded. The network has been named *q* myNetwork */q*, and basic information about the network is output. 
*

*
Similarly, the network may be input using format 2 and 3 as given in parameter files ``paras-input-network2.txt``and ``paras-input-network3.txt``respectively.
*

*/subsection*

*************