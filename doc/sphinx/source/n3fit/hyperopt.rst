.. _hyperoptimization:

================================
Hyperoptimization algorithm
================================

Motivation
----------
While the methodology used up to the 3.1 release of NNPDF considerably reduced the dependency on the
functional form of the PDFs compared to other collaborations, there existed a bias regarding the choice of hyperparameters that define
the NNPDF neural network and optimization strategy.

The hyperparameters of the model used by ``n3fit`` are determined by running a grid scan on the
relevant parameters. Together with an appropriate
figure of merit this grid search or *hyperparameter scan* will minimize the bias of the network
by finding the best one for each possible situation.

The final goal is for the methodology to be robust enough that a change in the physics
(fitted experiments, choice of basis, choice of constraints, ...) depends only on a new run of the
hyperparameter scan to be functional.

It is important to remember that the best hyperparameter combination is not necessarily the one that
produces the minimal training/validation :math:`\chi^2`. In fact, looking for the minimal :math:`\chi^2` is known to
produce overlearning even when optimizing on the validation loss, as can be seen
`here <https://vp.nnpdf.science/yG3XvinBQriLdqqTAHg3Sw==/>`_.

Despite producing a very good :math:`\chi^2`, the previous fit will fail when challenged with new
unseen data. This needs to be accounted for in the figure of merit of the hyperoptimization.

The desired features of this figure of merit can be summarized as:

1. Produce a low :math:`\chi^2` for both fitted experiments and non-fitted experiments.
2. Be stable upon random fluctuations.
3. Be reliable even when the number of points is not very large.


K-folding cross-validation
--------------------------
A good compromise between all previous points is the usage of the cross-validation technique
usually known as `k-folds <https://web.stanford.edu/~hastie/Papers/ESLII.pdf#page=260>`_.

In its most general form, we take all data points that enter the fit and break them down into *k*
partitions. Then, for every combination of hyperparameters, we do *k* fits leaving out a different
partition each time. We then use this partition to evaluate the goodness of the fit for each of the *k* fits and construct,
with these results, a reward function for the combination of hyperparameters.

In the NNPDF implementation of k-folding, each of the data points can be identified with a dataset.
Note that during the fit we perform the usual training-validation split within each dataset and use it for
stopping.

The choice of this method for selecting the hyperparameters of the NNPDF fitting methodology
has been discussed `in the mailing list <https://lists.cam.ac.uk/mailman/private/ucam-nnpdf/2020-March/msg00066.html>`_.
Some public discussion about the different hyperoptimization techniques that have been used and
tested during the development of ``n3fit`` can be found in `public slides <http://n3pdf.mi.infn.it/wp-content/uploads/2019/10/JCruz-Martinez_Mexico_102019.pdf>`_
as well as in `internal presentations <https://www.wiki.ed.ac.uk/display/nnpdfwiki/Amsterdam+Feb+2020+NNPDF+Collaboration+Meeting+agenda?preview=/432523942/436448892/juanCM.pdf>`_.


The choice of figure of merit is still under development, but we have several possibilities.

1. By default we take the combination that produces the best average for the partitions' :math:`\chi^2`.

.. math::
    L_{\rm hyperopt} = \frac{1}{N_{k}} \sum \chi^2

An example of a DIS fit using this loss function can be found here: [`best average <https://vp.nnpdf.science/iAaUMPgsTKyngsK5haLYMw==>`_]. It can be selected in the runcard using the target ``average``.

2. We can take the combination that produces the *best* *worst* loss.

.. math::
    L_{\rm hyperopt} = max(\chi^2)

An example of a DIS fit using this loss function can be found here: [`best worst <https://vp.nnpdf.science/0sWyhJZGQbuezEc7nMGATQ==>`_]. It can be selected in the runcard using the target ``best_worst``.

3. We can take the most stable combination which gets the loss under a certain threshold.

.. math::
   L_{\rm hyperopt} = \left\{
  \begin{array}{lr}
         std(\{\chi^{2}\}) & \text{  if } avg(\chi^2) < \text{ threshold } \\
         \infty & \text{otherwise}
  \end{array}
  \right.

An example of a DIS fit using this loss function with the threshold :math:`\chi^2` set to 2.0
can be found here: [`best std <https://vp.nnpdf.science/vcPtqM8KSXCVB2GheENd8Q==>`_].
It can be selected in the runcard using the target ``std``.

As observed, for DIS fits we obtain fits of similar quality using these losses.
This is not unexpected but it is a good test of the robustness of the method.

While this method is much more robust that the previously used "test set" (which is
similar to doing the limit :math:`k\rightarrow 1`) we can still find overfitting configurations.
For instance, if one of the metrics gives a much more complicated network structure,
overfitting is expected. Here's an example where, after 10000 hyperparameter trials,
the network structure had an order of magnitude more free parameters than normal,
in the case of the best average loss function:
[`best avg overlearned <https://vp.nnpdf.science/AQpgs2SyRbGlNqSnWWvMJw==>`_].


.. _hyperextrapolation-label:

Creating partitions
-------------------
The K-folding method is based on the creation of several partitions such that we can evaluate
how the fit would behave on completely unseen data.
The choice of this partitions is completely arbitrary, but defines the method completely.
Here we list some important considerations
to be taken into account when constructing these partitions.

- The reward function of the partitions must be comparable.

All loss functions implemented in ``n3fit`` for the optimization of hyperparameters use the reward
of all partitions as if they were equivalent.
When they are not equivalent the ``weight`` flag should be used (see :ref:`hyperoptrc-label`)



- Not all datasets should enter a partition: beware of extrapolation.

Beyond the last dataset that has entered the fit we find ourselves in what is usually known as
the extrapolation region. The behaviour of the fit in this region is not controlled by any data but
rather by the choice of preprocessing exponents (:math:`\alpha` at small x, :math:`\beta` at large x).
For this reason, if a dataset is included in a partition which however falls in the extrapolation region of the fit,
its loss function will be determined by these exponents (which are randomly chosen)
rather than by the hypeparameter combination.

The general rule that we follow is to always include in the fit the lowest-x dataset that determines
each of the PDF functions.
This means that no partition has datasets which falls in the extrapolation regions.
As a practical proxy-rule we can classify the datasets by process type and exclude from the partitioning
the ones that reach the lowest value of x.


Interpretation of results
-------------------------

While performing the hyperparameter scan we found that optimizing by only looking at the validation
loss produced results which would usually be considered overfitted: very low training and validation
:math:`\chi^2` but very complex replica patterns. Thanks to the high performance of the ``n3fit`` procedure the
usual within-dataset cross-validation algorithm used in the NNPDF framework was not enough to prevent overlearning
for all architectures.

The cross-validation implemented in NNPDF is successful in avoiding the learning of the noise within
a dataset. However, we observe that this choice is not enough to prevent overfitting due to
correlations between points in the same dataset when using hyperopt with ``n3fit``.

For hyperopt we have implemented k-folding cross-validation.
This method works by refitting with the same set of parameters several times (k times) each time leaving out
a partition of the datasets.
By using this method we reduce the bias associated with a particular choice of the datasets to leave out,
while at the same time, refitting with the same set of parameters allows us to assess the stability of the
particular combination of hyperparameters.

Implementation in ``n3fit``
---------------------------

The hyperparameter scan capabilities are implemented using the `hyperopt <https://github.com/hyperopt/hyperopt>`_ framework which
systematically scans over a selection of parameter using Bayesian optimization and measures model
performance to select the best architecture.
A `Jupyter Notebook is provided <https://github.com/NNPDF/tutorials/blob/master/hyperparameter%20scan/Hyperparameter%20scan.ipynb>`_
with a practical example of the usage of the hyperopt framework. This example is a simplified version
of the hyperparameter scan used in ``n3fit``.
The hyperopt library implements the tree-structured Parzen estimator algorithm
which is a robust sequential-model-based optimization approach `[SMBO] <https://en.wikipedia.org/wiki/Hyperparameter_optimization>`_.

We optimize on a combination of the best validation loss and the stability of the fits. In other words,
we select the architecture that produces the lowest validation loss after we trim those
combinations which are deemed to be unstable.

.. note::
    The fits done for hyperoptimization are one-replica fits. We take advantage of the
    stability of the Gradient Descent and of the fact that the difference between set of hyperparameters
    is small. This is a trade-off as we sustain a loss of "accuracy" (as some very ill-behaved replicas
    might destroy good sets of parameters) in exchange for being able to test many more parameters in
    the same time. Once a multireplica ``n3fit`` is implemented we can hyperoptimize without having to
    rely on the one-replica proxy and without a loss of performance.


From the fitting point of view, the implementation of the k-folding is done by setting all experimental
data points from the fold to 0 and by masking the respective predictions from the Neural Network to 0.
In the code this means that during the data-reading phase ``n3fit`` also creates one mask per k-fold
per experiment to apply to the experimental data before compiling the Neural Network.
Note that this is not a boolean mask that drops the points but rather it just sets the data to 0.
The reason for doing it in this way is to minimize the number of things that change when doing a
hyperparameter scan with respect to a fit.

Implementation in ``validphys``
-------------------------------

A ``hyperscan`` object is also available from ``validphys`` which behaves as a special case of ``fit``.
It can be accessed and inspected through the validphys API (see :ref:`vpapi`).
The product of a hyperparameter scan are ``tries.json`` files which can be acccessed with the
``tries_files`` attribute.

.. code-block:: python

   from validphys.api import API
   hyperscan = API.hyperscan(hyperscan="test_hyperopt_fit_300621")


It is also possible to access a ``hyperscan`` by using the ``validphys`` loader with:

.. code-block:: python

        from validphys.loader import Loader
        l = Loader()
        hyperscan = l.check_hyperscan("test_hyperopt_fit_300621")


.. _pos-int-hyperopt:

Positivity and integrability
----------------------------

Since positivity is a hard constraint of the fit (i.e., a replica fit will not be marked as good
unless it passes the positivity constraints), it enters the hyperoptimization in a similar way.
There is no threshold, either the replica passes positivity or it doesn't, and if it doesn't
hyperopt will receive a failure instead of a fit (so the run will be discarded).

Integrability instead is implemented as a penalty.
In order to activate it it is necessary to add ``integrability`` to the penalties section
of the hyperoptimization namespace (see below).
In this case the integrability is implemented as an exponential penalty, this means that as
the "integrability number" grows, the test loss will grow as well, favouring replicas with
an "integrability number" below the chosen threshold.
For consistency the threshold used during hyperoptimization is read directly from the ``fitveto.py`` variable.


.. _hyperoptrc-label:

Practical Usage
---------------

.. note::
  An example runcard can be found at ``n3fit/runcards/Basic_hyperopt.yml``.

The partitions can be chosen by adding a ``kfold::partitions`` key to the runcard.

.. code-block:: yaml

    kfold:
        target: average
        verbosity:
            training: True
            kfold: True
        threshold: 5.0
        penalties:
            - saturation
            - patience
            - integrability
        partitions:
            - overfit: True
              datasets:
                - data_1
                - data_2
            - weight: 2.0
              datasets:
                - data_3
            - datasets:
                - data_4
                - data_5

The ``overfit`` flag, when applied to one of the partitions, introduces this partition in the
fitted data, i.e., the training and validation always include that partition and will work normally.
This is useful for very broad scans where we want to find an architecture which is able to
fit, without worrying about things like overlearning which might be a second-order problem.

The ``weight`` (default 1.0) is multiplied with the loss function of the partition for which it is set.
Note that the weight is applied before the threshold check.

The ``threshold_loss`` flag will make the fit stop if any of the partitions produces a loss greater
than the given threshold. This is useful for quickly discarding hyperparameter subspaces without
needing to do all ``k`` fits.

The ``verbosity`` dictionary allows fine control over what to report each 100 epochs. When both ``training``
and ``kfold`` are set to ``False``, nothing is printed until the end of the fit of the fold.
When set to ``True``, the losses for the training (training and validation) and for the partition are printed.

During hyperoptimization we might want to search for specific features, such as quickly fitting
(giving an incentive to quicker runs) or avoiding saturation (increasing the loss for models that
have produce saturation after a fit). New penalties can easily be added in the ``src/n3fit/hyper_optimization/penalties.py`` file.


The target function for minimization can be selected with the ``target`` key.
By default, and if no ``target`` is chosen, ``n3fit`` defaults to
the average of the loss function over the partition sets (``average``).

.. math::
    L_{\rm hyperopt} = \frac{1}{N_{k}} \sum L_{k}

New target functions can be easily added in the ``src/n3fit/hyper_optimization/rewards.py`` file.

The hyperoptimization procedure performed in `hep-ph/1907.05075 <https://arxiv.org/abs/1907.05075>`_
used a slightly different approach in order to avoid overfitting,
by leaving out a number of datasets to compute a "testing set".
The loss function was then computed as

.. math::
    L_{\rm hyperopt} = \frac{1}{2} (L_{\rm validation} + L_{\rm testing})

The group of datasets that were left out followed the algorithm :ref:`mentioned above<hyperextrapolation-label>` with only one fold:


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

These were chosen attending to their `process type` as defined in their :ref:`commondata files <exp_data_files>`.


Changing the hyperoptimization target
-----------------------------------

Beyond the usual :math:`\chi2`-based optimization figures above, it is possible to utilize other measures as the target for hyperoptimization.
One possibility is to use a :ref:`future test<futuretests>`-based metric for which the goal is not to get the minimum :math:`\chi2` but to get the same :math:`\chi2` (with PDF errors considered) for different datasets. The idea is that this way we select models of which the prediction is stable upon variations in the dataset.
In order to obtain the PDF errors used in the figure of merit it is necessary to run multiple replicas, luckily ``n3fit`` provides such a possibility also during hyperoptimization.

Take the following modifications to a normal hyperopt runcard
(note that for convenience we take the trials directly from a previous run, so we don't have to create a new
hyperopt configuration dictionary).

.. code-block:: yaml

        dataset_inputs:
        - {dataset: NMCPD_dw_ite, frac: 0.75}
        - {dataset: NMC, frac: 0.75}
        - {dataset: SLACP_dwsh, frac: 0.75}
        - {dataset: SLACD_dw_ite, frac: 0.75}
        - {dataset: BCDMSP_dwsh, frac: 0.75}
        - {dataset: BCDMSD_dw_ite, frac: 0.75}
        - {dataset: HERACOMBNCEP575, frac: 0.75}
        - {dataset: HERACOMBCCEM, frac: 0.75}
        - {dataset: HERACOMBCCEP, frac: 0.75}

        hyperscan_config:
          use_tries_from: 210508-hyperopt_for_paper

        kfold:
          target: fit_future_tests
          partitions:
          - datasets:
            - HERACOMBCCEP
            - HERACOMBCCEM
            - HERACOMBNCEP575
          - datasets:

        parallel_models: true
        same_trvl_per_replica: true

We can run this hyperparameter scan for 10 parallel replicas for 20 trials with:

.. code-block:: bash

   n3fit runcard.yml 1 -r 10 --hyperopt 20

The above runcard will, for a sample of 20 trials in ``210508-hyperopt_for_paper`` (according to their rewards),
run two fits of 10 replicas each.
The first fit will hide the data from HERA and the second one (an empty fold) will take into consideration all data.
In order to properly set up a future test the last fold (the future) is recommended to be left as an empty fold such that no data is masked out.
The figure of merit will be the difference between the :math:`\chi2` of the second fit to the folded data and the :math:`\chi2` of the first fit to the folded data *including* pdf errors
(the :ref:`future test<futuretests>` :math:`\chi2`).

.. math::
   L_{\rm hyperopt} = \chi^{2}_{(1) \rm pdferr} - \chi^{2}_{(2)}


Restarting hyperoptimization runs
---------------------------------

In addition to the ``tries.json`` files, hyperparameter scans also produce ``tries.pkl`` `pickle <https://docs.python.org/3/library/pickle.html>`_ files,
which are located in the same directory as the corresponding ``tries.json`` file.
The generated ``tries.pkl`` file stores the complete history of a previous hyperoptimization run, making it possible to resume the process using the ``hyperopt`` framework.
To achieve this, you can use the ``--restart`` option within the ``n3fit`` command, e.g.,:

.. code-block:: bash

   n3fit runcard.yml 1 -r 10 --hyperopt 20 --restart

The above command example is effective when the number of saved trials in the ``test_run/nnfit/replica_1/tries.pkl`` is
less than ``20``. If there are ``20`` or more saved trials, ``n3fit`` will simply terminate, displaying the best results.


Running hyperoptimizations in parallel with MongoDB
---------------------------------------------------

In NNPDF, you can effectively run hyperoptimization experiments in parallel using `MongoDB <https://www.mongodb.com>`_.
This functionality is provided by the :class:`~n3fit.hyper_optimization.mongofiletrials.MongoFileTrials` class,
which extends the capabilities of `hyperopt <https://github.com/hyperopt/hyperopt>`_'s `MongoTrials` and enables the
simultaneous evaluation of multiple trials.

To set up and run a parallelized hyperopt search, follow these steps:

 1. **Instantiate the MongoDB database:** Start by setting up the database in your current directory.
 This database is referred to as ``hyperopt-db`` in the following instructions. You can initiate it with the command:

  .. code-block:: bash

    mongod --dbpath ./hyperopt-db

  By default, ``mongod`` uses port ``27017``. This is also the default port for the ``n3fit --db-port`` option.
  If you wish to use a different port, specify it as follows: ``mongod --dbpath ./hyperopt --port YOUR_PORT_NUMBER``.

 2. **Launch NNPDF with MongoDB integration:** Open a new command prompt and run ``n3fit`` with the desired configuration:

  .. code-block:: bash

    n3fit hyper-quickcard.yml 1 -r N_replicas --hyperopt N_trials --parallel-hyperopt --num-mongo-workers N

  Here, ``N`` represents the number of MongoDB workers you wish to launch in parallel.
  Each mongo worker handles one trial in Hyperopt. So, launching more workers allows for the simultaneous calculation of a greater number of trials.
  Note that there is no need to manually launch mongo workers, as the ``hyperopt-mongo-worker`` command is automatically
  executed by the :meth:`~n3fit.hyper_optimization.mongofiletrials.MongoFileTrials.start_mongo_workers` method.
  By default, the ``host`` argument is set to ``localhost``, and the database is named ``hyperopt``.
  If necessary, you can modify these settings using the ``n3fit --db-host`` or ``n3fit --db-name`` options.


.. note::

  Unlike in serial execution, parallel hyperoptimization runs do not generate ``tries.pkl`` files.
  To resume an experiment, simply retain the MongoDB database created during your previous run.
  Then, follow steps 1 and 2 as described above to restart the experiment.
