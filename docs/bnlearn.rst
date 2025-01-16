.. _bnlearn:

bnlearn network
===============

The default and recommended Bayesian network in BayesNetty is given by the bnlearn algorithm.
All future extensions are intended to be built upon this approach. For a given data set and network structure the likelihood can be calculated under specific distributional assumptions,
namely that discrete nodes follow a multinomial distribution and continuous nodes a normal distribution, with distributional parameters determined by the values of the incoming parent nodes.
The manner in which the likelihood is calculated can vary between Bayesian network algorithms.
See :cite:`bnlearn` and :cite:`bnlearn2` for further details of bnlearn methodology and R package.


.. _bnlearn-score:

Network score
-------------

The network score for a bnlearn network may be set to either the log likelihood, AIC or BIC using the ``-input-network-score`` option, see :ref:`input network options<input-network-options>`.

**NOTE:** The BIC network score is based on the definition used by bnlearn (see :cite:`bnlearn`) such that :math:`\text{BIC} = \log(L) - \frac{d}{2}\log(n)`, where :math:`L` is the likelihood
of the network for the given data set, :math:`d` is the number of parameters and :math:`n` is the number of individuals. This is the original definition used by Swartz in 1978,
see :cite:`swartz:1978`, rather than subsequent definitions of BIC which are multplied by negative two (for example see :cite:`wit:2012`). Therefore in BayesNetty the BIC will
always be negative and higher values of the network score imply a better fit network. (Whichever definition of BIC that one considers, the closer the BIC is to zero the better the model fit.)

The AIC network score in BayesNetty is defined similarly to the BIC such that :math:`\text{AIC} = \log(L) - d`, where :math:`L` and :math:`d` are defined as above.
Therefore higher values of the negatively valued AIC (closer to zero) imply a better network fit to the given data set.


Naturally if only the log likelihood is used then higher values imply a better network fit. 
