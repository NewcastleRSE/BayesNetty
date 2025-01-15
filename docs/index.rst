.. BayesNetty documentation master file, created by
   sphinx-quickstart on Mon Dec  9 10:03:57 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to BayesNetty's webpage!
================================

Please use the menu to the left to navigate through the documentation for BayesNetty or use the PDF of these webpages, `bayesnetty.pdf <https://github.com/NewcastleRSE/BayesNetty/blob/main/docs/_build/latex/bayesnetty.pdf>`_.

Citation
--------

Please use the following papers :cite:`howey:etal:21` and :cite:`howey:etal:20` to reference the BayesNetty software as follows:

.. code-block:: none

   @article{howey:etal:21,
   author={Howey, R AND Clark, A D AND Naamane, N AND Reynard, L N AND Pratt, A G AND Cordell, H J},
   title={{A Bayesian network approach incorporating imputation of missing data enables exploratory analysis of complex causal biological relationships}},
   journal={PLOS Genetics},
   month = {September},
   volume={17},
   number={9},
   doi = {10.1371/journal.pgen.1009811},
   url = {https://doi.org/10.1371/journal.pgen.1009811},
   pages={e1009811},
   year={2021}}

   @article{howey:etal:20,
   author = {Howey, R AND Shin, S-Y AND Relton, C AND Davey Smith, G AND Cordell, H J},
   journal = {PLOS Genetics},
   title = {Bayesian network analysis incorporating genetic anchors complements conventional Mendelian randomization approaches for exploratory analysis of causal relationships in complex data},
   year = {2020},
   month = {March},
   volume = {16},
   url = {https://doi.org/10.1371/journal.pgen.1008198},
   pages = {1-35},
   number = {3},
   doi = {10.1371/journal.pgen.1008198}}


BayesNetty is copyright, 2015-present Richard Howey, GNU General Public License, version 3.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   introduction
   installation
   using
   parallel
   input-data
   input-network
   bnlearn
   deal
   calc-score
   calc-posterior
   search-models
   average-network
   impute-data
   estimate-impute
   calc-recall-precision
   sim-data
   markov-blanket
   output-network
   output-priors
   output-posteriors
   plot-network
   references

.. _contact:

Contact
-------

Please contact `Richard Howey <https://www.staff.ncl.ac.uk/richard.howey/>`_ with any queries about the BayesNetty software.
