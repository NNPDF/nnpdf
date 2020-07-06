================================ 
Hyperoptimization algorithm
================================


Idea 
----

The main advantages of the points presented above consists in the possibility to test several models
in a fraction of time in comparison to the ``nnfit`` framework.

This is of key importance for a proper hyper-parameter scan where everything is potentially
interconnected.

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

.. code-block:: yml
    
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

The loss function is then computed as the average of the loss function over the partition sets.

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
