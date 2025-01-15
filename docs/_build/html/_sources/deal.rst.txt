 .. _deal:

deal network
============

The deal Bayesian network approach was developed by :cite:`deal_paper` as an approach to model mixed discrete/continuous networks.
It calculates the likelihood differently to bnlearn. However we found several issues with the method, not least that it is no longer actively supported.
Therefore, it is not recommended to use a deal network for network analyses and is included only for comparison purposes.  

.. _imaginary-sample-size:

Imaginary sample size
---------------------

When analysis is performed with a deal network the imaginary sample size (ISS) must be set. The ISS reflects how much confidence
we have in the (in)dependencies expressed in the assumed prior network. This can be set using the `-input-network-imaginary-sample-size` option, see :ref:`input-network-options`.
The results given by deal have been found to be very sensitive to the setting of this parameter and there is no obvious "good" default setting.

The network score in a deal network is based upon the log likelihood and so higher values imply a better network fit to the given data set.
