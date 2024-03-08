```eval_rst
.. _n3fit-usage:
```

How to run a PDF fit
====================


The user should perform the steps documented below in order to obtain a complete
PDF fit using the latest release of the NNPDF fitting code: ``n3fit``.
The fitting methodology is detailed in the [Methodology](methodology) page.

- [Preparing a fit runcard](#preparing-a-fit-runcard)
- [Running the fitting code](#running-the-fitting-code)
- [Upload and analyse the fit](#upload-and-analyse-the-fit)


```eval_rst
.. _prepare-fits:
```

Preparing a fit runcard
-----------------------

The runcard is written in YAML. The runcard is the unique identifier of a fit
and contains all required information to perform a fit, which includes the
experimental data, the theory setup and the fitting setup.

A detailed explanation on the parameters accepted by the ``n3fit`` runcards
can be found in the [detailed guide](runcard-detailed).

For newcomers, it is recommended to start from an already existing runcard,
example runcards (and runcard used in NNPDF releases) are available at
[``n3fit/runcards``](https://github.com/NNPDF/nnpdf/tree/master/n3fit/runcards).
The runcards are mostly self explanatory, see for instance below an
example of the ``parameter`` dictionary that defines the Machine Learning framework.

```yaml
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
```

The runcard system is designed such that the user can utilize the program without having to
tinker with the codebase.
One can simply modify the options in ``parameters`` to specify the
desired architecture of the Neural Network as well as the settings for the optimization algorithm.

An important feature of ``n3fit`` is the ability to perform [hyperparameter scans](hyperoptimization),
for this we have also introduced a ``hyperscan_config`` key which specifies
the trial ranges for the hyperparameter scan procedure.
See the following self-explanatory example:
```yaml
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
```

It is also possible to take the configuration of the hyperparameter scan from a previous
run in the NNPDF server by using the key `from_hyperscan`:
```yaml
hyperscan_config:
  from_hyperscan: 'some_previous_hyperscan'
```

or to directly take the trials from said hyperscan:
```yaml
hyperscan_config:
  use_tries_from: 'some_previous_hyperscan'
```


```eval_rst
.. _run-n3fit-fit:
```

Running the fitting code
------------------------

After successfully installing the ``n3fit`` package and preparing a runcard
following the points presented above you can proceed with a fit.

1. Prepare the fit: ``vp-setupfit runcard.yml``. This command will generate a
    folder with the same name as the runcard (minus the file extension) in the
    current directory, which will contain a copy of the original YAML runcard.
    The required resources (such as the theory and t0 PDF set) will be
    downloaded automatically. Alternatively they can be obtained with the
    ``vp-get`` tool.

2. The ``n3fit`` program takes a ``runcard.yml`` as input and a replica number, e.g.
```n3fit runcard.yml replica``` where ``replica`` goes from 1-n where n is the
maximum number of desired replicas. Note that if you desire, for example, a 100
replica fit you should launch more than 100 replicas (e.g. 130) because not
all of the replicas will pass the checks in ``postfit``
([see here](postfit-selection-criteria) for more info).

3. Wait until you have fit results. Then run the ``evolven3fit`` program once to
evolve all replicas using DGLAP. The arguments are ``evolven3fit evolve
runcard_folder``.

4. Wait until you have results, then use ``postfit number_of_replicas
runcard_folder`` to finalize the PDF set by applying post selection criteria.
This will produce a set of ``number_of_replicas + 1`` replicas. This time the
number of replicas should be that which you desire in the final fit (100 in the
above example). Note that the
standard behaviour of ``postfit`` can be modified by using various flags.
More information can be found at [Processing a fit](postfit).

It is possible to run more than one replica in one single run of ``n3fit`` by
using the ``--replica_range`` option. Running ``n3fit`` in this way increases the
memory usage as all replicas need to be stored in memory but decreases disk load
as the reading of the datasets and fktables is only done once for all replicas.


If you are planning to perform a hyperparameter scan just perform exactly the
same steps by adding the ``--hyperopt number_of_trials`` argument to ``n3fit``,
where ``number_of_trials`` is the maximum allowed value of trials required by the
fit. Usually when running hyperparameter scan we switch-off the MC replica
generation so different replicas will correspond to different initial points for
the scan, this approach provides faster results. We provide the ``vp-hyperoptplot``
script to analyse the output of the hyperparameter scan.


Output of the fit
-----------------
Every time a replica is finalized, the output is saved to the ```runcard/nnfit/replica_$replica```
folder, which contains a number of files:

- ``chi2exps.log``: a json log file with the χ² of the training every 100 epochs.
- ``runcard.exportgrid``: a file containing the PDF grid.
- ``runcard.json``: Includes information about the fit (metadata, parameters, times) in json format.

``` note:: The reported χ² refers always to the actual χ², i.e., without positivity loss or other penalty terms.
```



```eval_rst
.. _upload-fit:
```

Upload and analyse the fit
--------------------------
After obtaining the fit you can proceed with the fit upload and analisis by:

1. Uploading the results using ``vp-upload runcard_folder`` then install the
fitted set with ``vp-get fit fit_name``.

2. Analysing the results with ``validphys``, see the [vp-guide](../vp/index).
Consider using the ``vp-comparefits`` tool.



Performance of the fit
----------------------
The ``n3fit`` framework is currently based on [Tensorflow](https://www.tensorflow.org/) and as such, to
first approximation, anything that makes Tensorflow faster will also make ``n3fit`` faster.

``` note:: Tensorflow only supports the installation via pip. Note, however, that the TensorFlow pip package has been known to break third party packages. Install it at your own risk. Only the conda tensorflow-eigen package is tested by our CI systems.
```

When you install the nnpdf conda package, you get the [tensorflow-eigen](https://anaconda.org/anaconda/tensorflow-eigen) package, which is not the default.
This is due to a memory explosion found in some of the conda mkl builds.

If you want to disable MKL without installing ``tensorflow-eigen`` you can always set the environment variable ``TF_DISABLE_MKL=1`` before running ``n3fit``.
When running ``n3fit`` all versions of the package show similar performance.


When using the MKL version of tensorflow you gain more control of the way Tensorflow will use
the multithreading capabilities of the machine by using the following environment variables:

```bash

KMP_BLOCKTIME=0
KMP_AFFINITY=granularity=fine,verbose,compact,1,0
```

These are the best values found for ``n3fit`` when using the mkl version of Tensorflow from conda
and were found for TF 2.1 as the default values were suboptimal.
For a more detailed explanation on the effects of ``KMP_AFFINITY`` on the performance of
the code please see [here](https://software.intel.com/content/www/us/en/develop/documentation/cpp-compiler-developer-guide-and-reference/top/optimization-and-programming-guide/openmp-support/openmp-library-support/thread-affinity-interface-linux-and-windows.html).

By default, ``n3fit`` will try to use as many cores as possible, but this behaviour can be overriden
from the runcard with the ``maxcores`` parameter. In our tests the point of diminishing returns is found
at ``maxcores=4``.

Note that everything stated above is machine dependent so the best parameters for you might be
very different. When testing, it is useful to set the environmental variable ``KMP_SETTINGS`` to 1
to obtain detailed information about the current variables being used by OpenMP.

Below we present a benchmark that have been run for the Global NNPDF 3.1 case, as found in the
example runcards [folder](https://github.com/NNPDF/nnpdf/tree/master/n3fit/runcards).

Settings of the benchmark:
  - TF version: tensorflow-eigen from conda, TF 2.2
  - NNPDF commit: [f878fc95a4f32e8c3b4c454fc12d438cbb87ea80](https://github.com/NNPDF/nnpdf/commit/f878fc95a4f32e8c3b4c454fc12d438cbb87ea80)
  - Number of epochs: 5000
  - maxcores: 4
  - no early stopping

Hardware:
  - Intel(R) Core(TM) i7-6700 CPU @ 4.00GHz
  - 16 GB RAM 3000 MHz DDR4

Timing for a fit:
  - Walltime: 397s
  - CPUtime: 1729s

Iterate the fit
---------------

It may be desirable to iterate a fit to achieve a higher degree of convergence/stability in the fit.
To read more about this, see [How to run an iterated fit](run-iterated-fit).

QED fit
-------

In order to run a QED fit see [How to run a QED fit](run-qed-fit)
