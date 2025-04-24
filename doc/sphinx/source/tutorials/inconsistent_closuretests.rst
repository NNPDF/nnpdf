How to run an inconsistent closure test
=======================================

Modeling Inconsistencies
------------------------
In a realistic situation it may happen that some source of experimental systematics is overlooked 
or underestimated for a given dataset. In this case, the measured experimental values for this dataset 
may deviate from the true value by an amount that is not reflected by the experimental covariance matrix. 
This will then generate tensions between this dataset and the rest of the data. 
We call such a dataset “inconsistent”. It is interesting to ask how the NNPDF methodology behaves in such case,
and whether the inconsistency can be detected. For more details and a thorough review on the topic,
see the paper `arXiv:2503.17447 <https://arxiv.org/pdf/2503.17447>`_.


In order to study this situation in a closure test, we model the inconsistency as follow. 
We separate off the uncorrelated and correlated parts of the experimental covariance matrix 

.. math::

    (C)_{ij} = \delta_{ij} \sigma_i^{(\rm uncorr)} \sigma_j^{(\rm uncorr)} + \sum_{k=1}^{N_{\rm corr}} \sigma_{i,k}^{(\rm corr)}\sigma^{(\rm corr)}_{j,k}

where :math:`\sigma_i^{(\rm uncorr)}` and :math:`\sigma_{i,k}^{(\rm corr)}` denote respectively the uncorrelated and correlated systematics
for data point :math:`i`. 
We then define a rescaled covariance matrix 

.. math::

    (C^{\lambda})_{ij} = \delta_{ij} \sigma_i^{(\rm uncorr)} \sigma_j^{(\rm uncorr)} + \sum_{k=1}^{N_{\rm corr}} \lambda_{i,k}\sigma_{i,k}^{(\rm corr)} \lambda_{j,k}\sigma^{(\rm corr)}_{j,k},

in which correlated uncertainties have been rescaled.
In the closure test workflow we then generate Level 1 pseudodata using :math:`C`, namely 

.. math::
    L_1 = L_0 + \eta

with :math:`\eta \sim \mathcal{N}(0,C)`. The inconsistency is introduced at Level 2 when fitting the PDF
by using the rescaled covariance matrix :math:`C^{\lambda}` to propagate data uncertainties to the model space:  

.. math::
    L_2 = L_1 + \epsilon^{\lambda}.



Running an inconsistent closure test
------------------------------------

To run a closure test we require a standard closure test runcard. 
As already explained in the previous section the specification which controls 
closure test specific behaviour can be found under ``closuretest``.
An example of a typical level 1 or level 2 ``closuretest`` specification is given by

.. code:: yaml

  closuretest:
    filterseed  : 0   # Random seed to be used in filtering data partitions
    fakedata    : true
    fakepdf     : MMHT2014nnlo68cl
    fakenoise   : true

In order to run an inconsistent closure test, we need to first specify that the closure 
test should be inconsistent by adding the key ``inconsistent_fakedata`` to the ``closuretest``
specification.

.. code:: yaml

  closuretest:
    filterseed  : 0   # Random seed to be used in filtering data partitions
    fakedata    : true
    fakepdf     : MMHT2014nnlo68cl
    fakenoise   : true
    inconsistent_fakedata : true


In order to specify which datasets should be made inconsistent in the closure test
we need to add a new specification.
For instance

.. code:: yaml

    inconsistent_data_settings:

        treatment_names: [MULT]
        names_uncertainties: [CORR, SPECIAL]

        inconsistent_datasets:
            - HERA_NC_318GEV_EM-SIGMARED
            - HERA_NC_251GEV_EP-SIGMARED
            - HERA_NC_300GEV_EP-SIGMARED
            - HERA_NC_318GEV_EP-SIGMARED

        sys_rescaling_factor: 0.0

specifies that the inconsistency should be applied to the datasets
``HERA_NC_318GEV_EM-SIGMARED``, ``HERA_NC_251GEV_EP-SIGMARED``,
``HERA_NC_300GEV_EP-SIGMARED`` and ``HERA_NC_318GEV_EP-SIGMARED`` by multiplying
the MULT, CORR and SPECIAL (correlated intra-datasets) uncertainties by a factor of 0.0.

