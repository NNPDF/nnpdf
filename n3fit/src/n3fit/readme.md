Ciao everyone,

As promised here is the "beta" version of the new fitting code for everyone to have fun with it. We put all our code in the folder `n3fit` (i.e., everything below is referred to this folder).

# n3fit.py

This version is what we promised to make public in Varenna next August, in the end we decided to push also the hyperopt code since it was already ready rather than have a stripped-down version for Easter.


## Look before you leap

There are a few things to note before running this code.

### Differences with NNPDF

- Preprocessing:

This is perhaps one of the biggest differences with respect to NNPDF. The preprocessing variables are being fitted (within the parameters given in the runcard). There is room for new studies on this.

- The network:

We have constructed the network as a single net, i.e., all layers are fully connected. 

- Positivity and stopping:

Positivity is enforced and replicas are only accepted once they pass some positivity criteria. This is related to stopping in that the stop only occurs once the validation loss doesn't improve anymore **_and_** the positivity criteria is fulfilled.


- Seeds:

Everything is seeded at the level of the experiments, i.e., if you were to change the order of two experiments the results _should_ still be the same. 

### Performance, hardware usage

- CPU usage:

As said in Amsterdam, the code is quite fast and able to do a fit for a single replica in the order of 10 minutes / 1 hour for DIS / Global.  However the usage of the CPU power is not the most efficient

By default tensorflow will try to generate as many processes as logical threads your computer has (x2), however most of them are sitting idle most of the time, this is due to the amount of data we have available, it is actually not enough to actually give something to do to every core (_if_ you have a big enough number of powerful enough cores).

We have limited the code to 8-16 threads. The limitation is done in line 34,35 of `backend/keras_backend/debug_mode.py`. 

- Memory usage:

In our tests we find an usage of 2 GB of memory for DIS and 16 GB of memory for Global. I would personally recommend to be very generous with the memory requirements since otherwise you might end up swapping, which will greatly hinder performance.

**Warning:** If you were to run this code in a cluster with not very good sandboxing you will be still running in 8 processors even if you only asked for 2. Therefore if you are going to ask for a small number of processors I would recommend to change those settings in the code as well. You don't want to make your IT department sad.

_WIP:_ Eventually the resources management will be part of the runcard, but we are studying how to best reduce the impact (specially on memory) so it is still work in progress and so run it at your own risk.

### Caveats, things to be wary of 

The gradient descent algorithms are quite powerful and the possibility of running much faster means you can run much deeper networks. However, as a wise man once said, _with great power comes great responsability_. In particular in our tests we found **it is not difficult to run into overfitting issues**.

Overfitting is most easily done in DIS where the number of data points is small enough that is not difficult to generate a network with a number of parameters comparable to the number of data points. In this situation it is possible to find overfitting even with cross validation measures. We uploaded a report where the overfitting is quite obvious as an example: https://vp.nnpdf.science/EIA8y1JXQVmT9RYBLlUeJw==/


For truly reproducible runs you should set `debug: true` in the runcard as, by default, we allow tensorflow freedom to play with its own internal variables. Debug mode turns off hyperthreading and fixes tensorflow's internal variables.

The folder structure is not final, specially the files found in the top level folder. There is still some refactoring to do here so the explanation for the folders below might be outdated in a few weeks.

The code works as a validphys app but it is not integrated in validphys. It is not installed when you run `conda install nnpdf` or `make install`. At the moment this could be taken as "Release Candidate" but not as a final product.

At the moment since this is a beta version the final touchs of development are being in the n3pdf repository since there are a few pull requests open. Once we have merged all pull requests there we will move development to this one.


### Requirements:

- **NNPDF:** Make sure NNPDF is installed and that python is able to import validphys.
- **Python packages:** usually having NNPDF installed and working means you will have everything you need, possible pitfalls: ruamel.yaml
- **Keras & tensorflow:** if you installed nnpdf with conda I would suggest installing them also with conda to ensure compatibility between python versions (3.6/3.7 brings some changes that can affect tensorflow). Note: tensorflow-opt includes some optimizations for running on CPU so it is the recommended version.
- **Hyperopt:** in order to use hyperopt you need to install the hyperopt python package, pandas and seaborn.
- **Operative System:** we have run this code in many several different Linux machines (with wildly different Linux flavours) in both python 3.6 and 3.7. We haven't tested Mac OS but there is no reason it should not work.

## Folder Structure:

### Backend 

All backend dependent functions and classes are in the `backends` folder. At the moment there is only one backend implemented: keras (`backend/keras_backend`).

In order to implement a new different backend it should be enough to ensure it provides the same classes and methods.

### Layers

The custom layers that define the different pieces of the Neural Network are all in the `layers/` folder. Each layer is a different class extending `backends/?/MetaLayer.py`.

Note for instance the DIS layer `layers/DIS.py`. This layer is an operation implementing the convolution of one PDF with one FKTable (according to the active flavours).
The plan is to replace some of the operations in these layers by custom operations as we believe this will have an impact on the memory usage of the code.


### Base folder

The base folder is (temporarily) where the rest of the code lies, but let me break it down.

#### Main code:

- n3fit.py

Validphys app definition.

- fit.py

Driver of the fitting code, the hyperoptimization and the output of the code.

- constrains.py

Implements restrictions such as the momentum sum rule.

- model_gen.py

Implements functions dedicated to generate the different independent pieces of the Neural Network such as the Observable layers or the PDF layers.

- reader.py

Implements the reading of the datasets, experiments and positivity sets into python objects. This is the only place where libnnpdf is used. The hope is that most of the logic in this file will be substituted by PR #404.

- writer.py

Implements the functions dedicated to the output of the fit (crucially, .exportgrid)

- statistics.py

Class that saves information about the fitting such a number of epochs, losses of the different pieces, etc. It is also manages the stopping criteria.

- ModelTrainer.py

This class defines the behaviour of the fit, creates the Neural Network and performs the fit.


#### Hyperopt: 

These are files related to the hyperoptimization mode of the code. Below in the basic usage section there are some notes on how to use it in case someone wants to give a try.

- HyperAlgorith.py

- HyperScanner.py

- filetrials.py

- plotting.py

- create_testset.py


#### Runcards and report examples

We provide three runcards (in the `runcards` folder):

- Basic: Basic_runcard.yml, for quick testing with just two datasets
- DIS: PN3_DIS_example.yml (report example: https://vp.nnpdf.science/fpzXjidYR7SMrFBXpSacEA==)
- Global: PN3_Global_example.yml(report example: https://vp.nnpdf.science/BT98HzUSTxqg6_sDkkUx0A== )


#### Evolution code

Instead of using apfel replica per replica the code found in the folder `evolven3fit' will run the evolution for all replicas at once.
This code doesn't do any kind of sanity check or similar, as long as it finds a .gridinfo file, it will be happy. 

Usage:

    make
    ./evolven3fit run_folder number_of_replicas
    

### Basic Usage:

The most basic usage is running a fit for a given runcard for a given replica:
    
    python3 n3fit.py runcard.yml replica_number

This will generate a folder structure akin to the one nnfit currently generates. It is not necessary to run vp-setupfit (if you do it, it will not hurt, you will need it anyway for generating a report at some point in the pipeline).

Example:

    python3 n3fit.py runcards/NNPDF31_nlo_as_0118_basic.yml 1

will generater a folder `NNPDF31_nlo_as_0118_basic/nnfit/replica_1` containing a number of files with arclengths, chi2, etc.


#### Running multiple replicas sequentially:
With the option -r is it possible to run more replicas without the need for loading the data each time. This is just a convenience feature

    python3 n3fit.py runcards/NNPDF31_nlo_as_0118_basic.yml 1 -r 4

will run replicas 1 to 4 one after the other.

#### Hyperoptimization:
The option --hyperopt will activate the hyperoptimizer

    python3 n3fit.py runcards/NNPDF31_nlo_as_0118_basic.yml 1 --hyperopt 100

will run 100 iterations of the hyperoptimizer. We also include a script that will plot the results as a `scan.pdf' file.

    python3 plotting.py NNPDF31_nlo_as_0118_basic

As mentioned before, it is important to avoid overfitting, relying blindly in hyperopt, however, can result in models prone to overfitting. For hyperopt a model that overfits is very good, the most straightforward way of solving this issue is to create a set of datasets that don't enter the fitting or the validation. We include a script to do this in a systematic way with the following option:

    python3 n3fit.py runcards/NNPDF31_nlo_as_0118_basic.yml 1 --create_test_card

This will generate a TEST runcard with a TEST set included which can be use for hyperoptimization.
    

#### Example: from a fit to a validphys report


In order to generate a fit it is necessary to run the code for a number of replicas. Let's say we want a fit of the runcard DIS.yml with 20 replicas:

```
for replica in {1..20}
do
    ./n3fit.py DIS.yml ${replica}
done
```

After this we'll have a `DIS/nnfit' folder with 20 replicas. Now in order to generate the fit we use the following workflow:

First some legacy code, it's only needed for the other scripts not to complain

    vp-setupfit DIS.yml

Now we can run evolven3fit with the folder were the replicas are included and the maximum number of replicas to find.

    ./evolven3fit DIS 20`

From this point onwards the steps are not different from before, for completeness they are also included here:

After running the evolution code we run postfit to generate the final PDF, for that we give the number of members we want to find and the folder.

    postfit 15 DIS

After postfit has finished you will want the results to the NNPDF share folder, assuming a conda installation:  

    cp -RL DIS $CONDA_PREFIX/share/NNPDF/results/

Finally, for generating the plots you need to use comparefits and follow the instructions, for the DIS runcard we provide you can use 181101-001-sc as the reference fit.

    vp-comparefits DIS 181101-001-sc --title comparison --author $USER --keywords n3fit -o my_DIS_report


### Going forward:

This is a list of ideas we have so people can start working on things if they feel like it. Some of them are necessary for the full integration with NNPDF (in particular testing and implementing the new fkparser) other can be interesting or fun side projects.

- Testing the fkparser.

Modifying reader to use the new fkparsers implemented in PR #404 instead of using the C functions.
All calls to libnnpdf are concentrated in this file so that the changes are easy to implement.

- Closure tests

We have not run any closure test or implemented a way of doing it in a systematic way. It could be an interesting exercise.


- Implement new backends

As said before, it should be straightforward (a straight paths can be very long, complicated and full of holes, but the path itself has not many curves).

- Run with other algorithms

At the moment we have run with some of the algorithms included in the default keras installation, but the space of algorithms for machine learning is quite big. We have also limited ourselves to Gradient Descent family, but that's not the end of the story.


- Study different strategies for preprocessing:

At the moment we are taking a very naive approach to the fitting of the processing: we let them vary freely within the strict limits set by the runcard. This may well be suboptimal.

#### Some other things

Some stuff that people is working on and is close to completion:

- Implementation of custom operators and gradients to improve performance.
- An visualizer for the fitting of the PDF in the form of an animation.


--------

Let us know if anything is unclear and feel free to open any issues you find!


                 Juan & Stefano
