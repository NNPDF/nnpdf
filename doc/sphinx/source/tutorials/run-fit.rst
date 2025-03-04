.. _n3fit-usage:

How to run a PDF fit
====================

The user should perform the steps documented below in order to obtain a complete
PDF fit using the latest release of the NNPDF fitting code: ``n3fit``.
The fitting methodology is detailed in :ref:`methodology`.

The three main points in this tutorial are:

- :ref:`Preparing a fit runcard <prepare-fits>`
- :ref:`Running the fitting code <run-n3fit-fit>`
- :ref:`Upload and analyse the fit <upload-fit>`
- :ref:`advance-run-fit`


.. _prepare-fits:

Preparing a fit runcard
-----------------------

The runcard is written in YAML. The runcard is the unique identifier of a fit
and contains all required information to perform and reproduce a fit, which includes the
experimental data, the theory setup and the fitting setup.
A detailed explanation on the parameters accepted by the ``n3fit`` runcards
can be found in the :ref:`n3fit detailed guide <runcard-detailed>`.

For newcomers, it is recommended to start from an already existing runcard,
example runcards (and runcard used in NNPDF releases) are available at
`n3fit/runcards <https://github.com/NNPDF/nnpdf/tree/master/n3fit/runcards>`_.

.. note::

  While we aim for the code to be both backwards and forwards compatible with respect to runcards,
  by setting sensible defaults when introducing new features,
  make sure that you are using a runcard tagged with the same version of the code you are using to
  avoid any surprises.

The runcards are mostly self explanatory, see for instance below an
example of the ``parameter`` dictionary that defines the Machine Learning methodology.

.. code:: yaml

  # runcard example
  ...
  parameters:
    nodes_per_layer: [15, 10, 8]
    activation_per_layer: ['sigmoid', 'sigmoid', 'linear']
    initializer: 'glorot_normal'
    optimizer:
      optimizer_name: 'RMSprop'
      learning_rate: 0.01
      clipnorm: 1.0
    epochs: 900
    positivity:
      multiplier: 1.05
      threshold: 1e-5
    stopping_patience: 0.30 # Ratio of the number of epochs
    layer_type: 'dense'
    dropout: 0.0
  ...

The runcard system is designed such that the user can utilize the program
without having to tinker with the codebase.
One can simply modify the options in ``parameters`` to specify the
desired architecture of the Neural Network as well as the settings for the optimization algorithm.


.. _run-n3fit-fit:

Running the fitting code
------------------------

After successfully installing the ``n3fit`` package and preparing a runcard
following the points presented above you can proceed with a fit.

1.  Prepare the fit using ``vp-setupfit``. This command will generate a
    folder with the same name as the runcard (minus the file extension) in the
    current directory, which will contain a copy of the original YAML runcard.
    The required resources will be downloaded, which includes:

      - The t0 PDF set (an LHAPDF object).
      - The FastKernel tables in the form of a ``theory_xxx.tgz`` file
      - The postfit evolution operator ``EKO.tar``.

    If the runcard requires to precompute some heavy objects shared among replicas,
    such as the theory covariance matrix, it will be done during this step.

::

  vp-setupfit <runcard>.yml

2.  Run the fit using ``n3fit``. The ``n3fit`` program takes a ``runcard.yml`` as input and a replica number, e.g.
    :code:`n3fit runcard.yml replica` where ``replica`` goes from 1-n where n is the
    maximum number of desired replicas. Note that if you desire, for example, a 100
    replica fit you should launch more than 100 replicas (e.g. 120) since not all of them will necessarily converge.
    While by default the code runs each replica separately, it is possible to run many replicas in parallel, see :ref:`parallel-label`.

::

  for i in {1..120} ; do
    n3fit <runcard>.yml $i
  done

3.  Once all replicas have finished, you need to run the ``evolven3fit`` program in order to
    evolve the PDF from the fitting scale to the whole range of scales needed to create an LHAPDF grid.
    This is done using the EKO library to perform DGLAP evolution.

::

  evolven3fit evolve <runcard>

4.  Finally, use ``postfit`` to finalize the PDF set by applying post selection criteria and compute the central replica.
    This will produce a set of ``number_of_replicas`` error replicas and one mean replica for a total of ``number_of_replicas+1``.
    The number of replicas should be that which you desire in the final fit (e.g., 100).
    Note that the standard behaviour of ``postfit`` can be modified by using various flags.
    More information can be found at :ref:`Processing a fit <postfit>`.

::

  postfit <number of desired replicas> <runcard>



Output of the fit
~~~~~~~~~~~~~~~~~
Every time a replica is finalized, the output is saved to the ```runcard/nnfit/replica_$replica```
folder, which contains a number of files:

- ``chi2exps.log``: a json log file with the χ² of the training every 100 epochs.
- ``runcard.exportgrid``: a file containing the PDF grid.
- ``runcard.json``: Includes information about the fit (metadata, parameters, times) in json format.

.. note::

  The reported χ² refers always to the actual χ², i.e., without positivity loss or other penalty terms.


.. _upload-fit:

Upload and analyse the fit
--------------------------
After obtaining the fit you can proceed with the fit upload and analysis by:

1.  *For members of NNPDF*, it is possible to upload the results to the nnpdf server using ``vp-upload runcard_folder`` then install the fitted set with ``vp-get fit fit_name``. Otherwise, copy or link to the results to the ``share/NNPDF/results/`` folder (usually under ``~/.local/share`` or ``${CONDA_PREFIX}/share/``.

2.  Analysing the results with ``validphys``, see the :ref:`vp-guide <vp-index>`.
    Consider using the ``vp-comparefits`` tool.


.. _advance-run-fit:

Advanced topics
---------------


Fit performance
~~~~~~~~~~~~~~~
The ``n3fit`` framework is currently based on `Keras <https://keras.io/>`_
and it is tested to run with the `Tensorflow <https://www.tensorflow.org/>`_
and `pytorch <https://pytorch.org>`_ backends.
This also means that anything that make any of these packages faster will also
make ``n3fit`` faster.
Note that at the time of writing, ``TensorFlow`` is approximately 4 times faster than ``pytorch``.

The default backend for ``keras`` is ``tensorflow``.
In order to change the backend, the environment variable ``KERAS_BACKENDD`` need to be set (e.g., ``KERAS_BACKEND=torch``).

The best results are obtained with ``tensorflow[and-cuda]`` installed from pip
and running ``n3fit`` in GPU, see :ref:`parallel-label`.

Iterate the fit
~~~~~~~~~~~~~~~

It may be desirable to iterate a fit to achieve a higher degree of convergence/stability in the fit.
To read more about this, see :ref:`How to run an iterated fit <run-iterated-fit>`.

QED fit
~~~~~~~

In order to run a QED fit see :ref:`How to run a QED fit <run-qed-fit>`.


Hyperparameter optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~

An important feature of ``n3fit`` is the ability to perform :ref:`hyperparameter scans <hyperoptimization>`,
for this we have also introduced a ``hyperscan_config`` key which specifies
the trial ranges for the hyperparameter scan procedure.
See the following self-explanatory example:

.. code:: yaml

  hyperscan_config:
      stopping: # setup for stopping scan
          min_epochs: 5e2  # minimum number of epochs
          max_epochs: 40e2 # maximum number of epochs
          min_patience: 0.10 # minimum stop patience
          max_patience: 0.40 # maximum stop patience
      positivity: # setup for the positivity scan
          min_multiplier: 1.04 # minimum lagrange multiplier coeff.
          max_multiplier: 1.1 # maximum lagrange multiplier coeff.
          min_initial: 1.0 # minimum initial penalty
          max_initial: 5.0 # maximum initial penalty
      optimizer: # setup for the optimizer scan
          - optimizer_name: 'Adadelta'
            learning_rate:
              min: 0.5
              max: 1.5
          - optimizer_name: 'Adam'
            learning_rate:
              min: 0.5
              max: 1.5
      architecture: # setup for the architecture scan
          initializers: 'ALL' # Use all implemented initializers from keras
          max_drop: 0.15 # maximum dropout probability
          n_layers: [2,3,4] # number of layers
          min_units: 5 # minimum number of nodes
          max_units: 50 # maximum number of nodes
          activations: ['sigmoid', 'tanh'] # list of activation functions

It is also possible to take the configuration of the hyperparameter scan from a previous
run in the NNPDF server by using the key ``from_hyperscan``:

.. code:: yaml

  hyperscan_config:
    from_hyperscan: 'some_previous_hyperscan'

or to directly take the trials from said hyperscan:

.. code:: yaml

  hyperscan_config:
    use_tries_from: 'some_previous_hyperscan'


If you are planning to perform a hyperparameter scan just perform exactly the
same steps as in :ref:`run-n3fit-fit` by adding the ``--hyperopt number_of_trials`` argument to ``n3fit``,
where ``number_of_trials`` is the maximum allowed value of trials required by the
fit. Usually when running hyperparameter scan we switch-off the MC replica
generation so different replicas will correspond to different initial points for
the scan, this approach provides faster results. We provide the ``vp-hyperoptplot``
script to analyse the output of the hyperparameter scan.
