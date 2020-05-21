```eval_rst
.. _n3fit-usage:
```

How to run a PDF fit
====================


The user should perform the steps documented below in order to obtain a complete
PDF fit:

- [Preparing a fit runcard](#preparing-a-fit-runcard)
- [Running the fitting code](#running-the-fitting-code)
- [Upload and analyse the fit](#upload-and-analyse-the-fit)


These steps are similar to those required to run the
[`nnfit` code](nnfit-usage) as `n3fit` is designed to maintain close
compatibility with it.

Preparing a fit runcard
-----------------------

The runcard is written in YAML. The runcard is the unique identifier of a fit
and contains all required information to perform a fit, which includes the
experimental data, the theory setup and the fitting setup.

The `n3fit` code accepts the same YAML keys used by `nnfit` and extra keys
required by the new techniques introduced in [Methodology](methodology).

In particular we introduce new keys in the `fitting` section such as:
- different seeds for the training and validation split (`trvlseed`), neural
network intialization (`nnseed`) and the MC data replica generation (`mcseed`).
- some debug flags to store or load the model in/from hd5 files (`save`,
`savefile`, `load`, `loadfile`, `plot`)
- a new `parameters` dictionnary with the model specifications based on the
arguments described in [Methodology](methodology).

See, as an example, the following self explanatory runcard fragment:
```yaml
# runcard example
...
fitting:
  trvlseed: 1
  nnseed: 2
  mcseed: 3
  save: False
  savefile: 'weights.hd5'
  load: False
  loadfile: 'weights.hd5'
  plot: False

  parameters: # This defines the parameter dictionary that is passed to the Model Trainer
    nodes_per_layer: [15, 10, 8]
    activation_per_layer: ['sigmoid', 'sigmoid', 'linear']
    initializer: 'glorot_normal'
    learning_rate: 0.01
    optimizer: 'RMSprop'
    epochs: 900
    pos_multiplier: 1.05
    pos_initial:  # believe the pos_lambda below
    stopping_patience: 0.30 # percentage of the number of epochs
    layer_type: 'dense'
    dropout: 0.0
...
```

On the other hand, we have also introduced a new `hyperscan` key which specifies
the trial ranges for the hyperparameter scan procedure. See the following self
explanatory example:
```yaml
hyperscan:
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


Finally, complete working examples of DIS-only and global fits are available at
the git repository in
[n3fit/runcards](https://github.com/NNPDF/nnpdf/tree/master/n3fit/runcards).

For a more detailed explanation on the parameters that are specific for the
`n3fit` runcard see the [detailed guide](runcard_detailed).

Running the fitting code
------------------------

After successfully installing the `n3fit` package and preparing a runcard
following the points presented above you can proceed with a fit.

1. Prepare the fit: `vp-setupfit runcard.yml` this command will generate a
`runcard_folder` folder in the current directory with a copy of the original
YAML runcard.  The required resources (such as the theory and t0 PDF) will be
downloaded automatically. Alternatively they can be obtained with the `vp-get`
tool.

2. The `n3fit` program takes a `runcard.yml` as input and a replica number, e.g.
```n3fit runcard.yml replica``` where `replica` goes from 1-n where n is the
maximum number of desired replicas.

3. Wait until you have fit results. Then run the `evolven3fit` program once to
evolve all replicas using DGLAP. The arguments are `evolven3fit runcard_folder
number_of_replicas`.

4. Wait until you have results, then use `postfit number_of_replicas
runcard_folder` to finalize the PDF set by applying post selection criteria.
This will produce a set of `number_of_replicas + 1` replicas.

It is possible to run more than one replica in one single run of `n3fit` by
using the ``--replica_range`` option. Running `n3fit` in this way increases the
memory usage as all replicas need to be stored in memory but decreases disk load
as the reading of the datasets and fktables is only done once for all replicas.


If you are planning to perform a hyperparameter scan just perform exactly the
same steps by adding the `--hyperopt number_of_trials` argument to `n3fit`,
where `number_of_trials` is the maximum allowed value of trials required by the
fit. Usually when running hyperparameter scan we switch-off the MC replica
generation so different replicas will correspond to different initial points for
the scan, this approach provides faster results. We provide the `n3Hyperplot`
script to analyse the output of the hyperparameter scan.


Output of the fit
-----------------
In the same fashion as `nnfit`, every time a replica is finalized a folder is
created in ```runcard/nnfit/replica_$replica```. This folder contains several
files which follow the same structure as `nnfit`:

- `runcard.exportgrid`: a file containing the PDF grid.
- `chi2exps.log`: a log file with the χ² of the training every 100 epochs.
- `runcard.preproc`: Empty file.
- `runcard.fitinfo`: Includes information about the fit. The first line
contains, in this order, the number of epochs, the validation χ², training
χ², experimental χ² and the state of the positivity. The second line the
arclength for each flavour.
- `runcard.time`: Includes the total time the fit took in CPU time and walltime.
The times are separated by the time of the actual fit and the time of the data
load.

``` note:: The reported χ² refers always to the actual χ², i.e., without positivity loss or other penalty terms.
```


Upload and analyse the fit
--------------------------
After obtaining the fit you can proceed with the fit upload and analisis by:

1. Uploading the results using `vp-uploadfit runcard_folder` then install the
fitted set with `vp-get fit fit_name`.

2. Analysing the results with `validphys`, see the [vp-guide](../vp/index).
Consider using the `vp-comparefits` tool.



Performance of the fit
----------------------
The `n3fit` framework is currently based on [Tensorflow](https://www.tensorflow.org/) and as such, to
first approximation, anything that makes Tensorflow faster will also make ``n3fit`` faster.

In our tests the bests results are obtained using the MKL-compiled version of Tensorflow as found by
default in Conda.

When using the MKL version the following environmental variables are relevant to control the
performance of Tensorflow.

```bash

KMP_BLOCKTIME=0
KMP_AFFINITY=granularity=fine,verbose,compact,1,0

```

When these variables are not set, `n3fit` will default to the values shown above.
For more detailed explanation on the effects of `KMP_AFFINITY` on the performance of
the code please see [here](https://software.intel.com/content/www/us/en/develop/documentation/cpp-compiler-developer-guide-and-reference/top/optimization-and-programming-guide/openmp-support/openmp-library-support/thread-affinity-interface-linux-and-windows.html)

By default, `n3fit` will try to use as many cores as possible, but this behaviour can be overriden
from the runcard with the `maxcores` parameter. In our tests the point of diminishing returns is found
at `maxcores=4`.

Note that everything stated above is machine dependent so the best parameters for you might be
very different. When testing it is useful to set the environmental variable `KMP_SETTINGS` to 1
to obtain detailed information about the current variables being used by OpenMP.

Below we list some benchmarks that have been run for the Global NNPDF 3.1 case, as found in the
example runcards [folder](https://github.com/NNPDF/nnpdf/tree/master/n3fit/runcards).
