================================ 
Hyperoptimization algorithm
================================

Motivation
----------
While the methodology used up to the 3.1 release of NNPDF considerable reduced the dependency on the
functional form of the PDFs, there exists a bias regarding the choice of hyperparameters that define
the NNPDF neural network and optimization strategy.

Of the main advantages introduced by the ``n3fit`` framework with respect to ``nnfit`` is the
possibility of running fits in a fraction of the time. This allow us to reduce the dependence of the
hyperparameters by running a grid scan on the relevant parameters. Together with an appropriate
figure of merit these grid search or *hyperparameter scan* will minimize the bias of the network
finding the best one for each possible situation.

The final goal is for the methodology to be robust enough that a change on the physics
(fitted experiments, choice of basis, choice of constraints, ...) depends only on a new run of the
hyperparameter scan to be functional.


Figure of merit
---------------
Our goal when running a hyperparameter scan is not just to find the hyperparameter combination that
produces the minimal :math:`\chi^2`. In fact, looking for the minimal :math:`\chi^2` is known to
produce overlearning even when optimizing on the validation loss, as can be seen
`here <https://vp.nnpdf.science/yG3XvinBQriLdqqTAHg3Sw==/>`_. 

Despite producing a very good :math:`\chi^2`, the previous fit will fail when challenged with new
non-seen data. This needs to be accounted for by the figure of merit.

The desired features of this figure of merit can be summarized as:

1. Produce a low :math:`\chi^2` for both fitted experiments and non-fitted experiments.
2. Be stable upon random fluctuations.
3. Be reliable even when the number of points is not extremely large.



K-folding cross-validation
--------------------------
A good compromise between all previous points is the usage of the cross validation technique
usually known as k-folds.

When k-folding  we take all datapoints in our fit (in this case one datapoint refers to one dataset)
and break it down into *k* partitions. Now, for every combination of hyperparameter we do *k* fits,
leaving out a different partition each time.
At the end of the fit we test the goodness of the fit with the data that was left out.
In this way all datasets are used exactly once for testing of the hyperparameters
(and *k-1* times for fitting).

For the fit we perform the usual training validation split within each dataset and use it for
stopping.


**Under construction:**
the choice of figure of merit is still under development, but we have several possibilities

1. We can take the combination that produces the best average for the partitions' :math:`\chi^2`

.. math::
    L_{hyperopt} = \frac{1}{N_{k}} \sum \chi^2

2. We can take the combination that produces the *best* *worst* loss

.. math::
    L_{hyperopt} = max(\chi^2)


Implementation
--------------

The hyperparameter scan capabilities are implemented using the hyperopt framework which
systematically scans over a selection of parameter using Bayesian optimization and measures model
performance to select the best architecture. The hyperopt library implements the tree-structured of
Parzen estimator which is a robust sequential model based optimization approach `[SMBO] <https://en.wikipedia.org/wiki/Hyperparameter_optimization>`_.

We optimize on a combination of the best validation loss and stability of the fits. In other words,
we select the architecture which produces the lowest validation loss after we trim those
combinations which are deemed to be unstable.

.. note:: 
    The fits done for hyperoptimization are one-replica fits. We take advantage of the
    stability of the Gradient Descent and of the fact that the difference between set of hyperparameters
    is small. This is a trade-off as we sustain a loss of "accuracy" (as some very ill-behave replicas
    might destroy good sets of parameters) in exchange for being able to test many more parameters in
    the same time. Once a multireplica ``n3fit`` is implemented we can hyperoptimize without having to
    rely in the one-replica proxy without a loss of performance.


Interpretation 
--------------

While performing the hyperparameter scan we found that optimizing only looking at the validation
loss produced results which would usually be considered overfitted: very low training and validation
chi2 and complex replica patterns. Thanks to the high performance of the ``n3fit`` procedure the
usual cross-validation algorithm used in the NNPDF framework was not enough to prevent overlearning
for all architectures.

The cross-validation implemented in NNPDF is successful on avoiding the learning of the noise within
a dataset. However, we observe that this choice is not enough to prevent overfitting due to
correlations within points in a same dataset when using hyperopt with ``n3fit``.

For hyperopt we have implemented k-folding cross-validation.
This method works by refitting with the same set of parameters several times (k times) each time leaving out
a partition of the datasets.
By using this method we reduce the bias associated with a particular choise of the datasets to leave out
while at the same time, refitting with the same set of parameters, allow us to assess the stability of the
particular combgination of hyperparameters.

The partitions can be chosen by adding a ``partitions`` key to the ``hyperscan`` dictionary.

.. code-block:: yaml
    
    kfold:
        partitions:
            - datasets:
                - data_1
                - data_2
            - datasets:
                - data_3
            - datasets:
                - data_4
                - data_5

An example runcard can be found at ``n3fit/runcards/Basic_hyperopt.yml``.

The loss function is currently computed as the average of the loss function over the partition sets.

.. math::
    L_{hyperopt} = \frac{1}{N_{k}} \sum (L_{k})

The hyperoptimization procedure performed in `hep-ph/1907.05075 <https://arxiv.org/abs/1907.05075>`_
used a different approach in order to avoid overfitting, by leaving out a number of datasets to compute
a "testing set". The loss function was then computed as:

.. math::
    L_{hyperopt} = \frac{1}{2} (L_{validation} + L_{testing})

The group of datasets that were left out were:


* NMC 
* BCDMSP 
* BCDMSD 
* HERACOMBNCEP460 
* H1HERAF2B 
* D0ZRap 
* CDFR2KT 
* D0WMASY
* ATLASZHIGHMASS49FB 
* CMSZDIFF12 
* ATLASTTBARTOT
