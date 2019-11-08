Methodology overview
====================

The goal of this document is to summarise from a conceptual point of view the main points which are different in comparison to the latest NNPDF (i.e. [NNPDF3.1](https://arxiv.org/abs/1706.00428)) methodology.

``` warning:: The default implementation of the concepts presented here are implemented with Keras and Tensorflow. The n3fit code inherits its features, so in this document we avoid the discussion of specific details which can be found in the Keras Documentation
```


``` note:: The final setup used in n3fit fits can be extracted from the runcards stored in nnpdf/n3fit/runcards.
```

This document contains a more specific discussion about the choices currently implemented in the `n3fit` package and discussed for the first time in [hep-ph/1907.05075](https://arxiv.org/abs/1907.05075). The Keras documentation can be found [here](https://keras.io/).

**Table of contents:**
- [Introduction](#introduction)
- [Neural network architecture](#neural-network-architecture)
- [Preprocessing](#preprocessing)
- [Optimizer](#optimizer)
- [Stopping](#stopping-algorithm)
- [Positivity](#positivity)
- [Hyperoptimization](#hyperoptimization-algorithm)

Introduction
------------

The approach presented here inherits the technology developed by the NNPDF collaboration in terms of fit pipeline but extends the possibility to test and improve fitting performance with modern techniques inspired by the deep learning community.

The `n3fit` code is designed in python and one of its main goals is to replace the `nnfit` program. It provides a simple abstraction layer which simplifies the life of developers when considering the possibility of adding new fitting algorithms.

In the following table we list some of the differences between both codes:

```eval_rst
+--------------------+---------------------------------+-------------------------------------+
| Component          | nnfit                           | n3fit                               |
+====================+=================================+=====================================+
| MC data generation | single seed                     | multi seed                          |
+--------------------+---------------------------------+-------------------------------------+
| Data management    | libnnpdf                        | same as nnfit                       |
+--------------------+---------------------------------+-------------------------------------+
| Neural net         | fixed architecture, per flavour |**single net, flexible architecture**|
+--------------------+---------------------------------+-------------------------------------+
| Preprocessing      | random fixed                    | **fitted in range**                 |
+--------------------+---------------------------------+-------------------------------------+
| Integration        | a posteriori per iteration      | **buildin in the model**            |
+--------------------+---------------------------------+-------------------------------------+
| Optimizer          | genetic optimizer               | **gradient descent**                |
+--------------------+---------------------------------+-------------------------------------+
| Stopping           | lookback                        | **patience**                        |
+--------------------+---------------------------------+-------------------------------------+
| Positivity         | penalty and threshold           | penalty, **PDF must be positive**   |
+--------------------+---------------------------------+-------------------------------------+
| Postfit            | 4-sigma chi2 and arclenght      | same as nnfit                       |
+--------------------+---------------------------------+-------------------------------------+
| Fine tuning        | manual                          | **semi-automatic**                  |
+--------------------+---------------------------------+-------------------------------------+
| Model selection    | closure test                    | closure test, **hyper optimization**|
+--------------------+---------------------------------+-------------------------------------+
```


``` note:: In the next sections we focus on the n3fit specifics marked in **bold**.
```

Neural network architecture
---------------------------

The main advantage of using a modern deep learning backend such as Keras/Tensorflow consists in the possibility to change the neural network architecture quickly as the developer is not forced to fine tune the code in order to achieve efficient memory management and PDF convolution performance.

The current `n3fit` code supports sequential dense networks with custom number of layers, nodes, activation functions and initializers from [Keras](https://keras.io/).

A big difference in comparison to `nnfit` is the number of neural networks involved in the fit. Here we use a **single neural network** model which maps the input (x, log x) to 8 outputs, nominally they correspond exactly the 8 PDF flavours defined in NNPDF3.1. The choice of using a single model is justified by the interest in preserving cross-correlations between flavours.

``` image:: figures/nn.png
```

The momentum sum rules are implemented as a **neural network layer** which computes the normalization coefficients for each flavour as a sum over a fixed grid of x points. The number and density of points in x is selected in such way that the final quality of the integrals are at least permille level in comparison to 1D integration algorithms.

The network initialization relies on modern deep learning techniques such as glorot uniform and normal (see [Keras initializers](https://keras.io/initializers/)), which have demonstrated to provide a faster convergence to the solution.

``` important:: Parameters like the number of layers, nodes, activation functions are hyper-paramters that require tuning.
```

Preprocessing
-------------
Preprocessing has been modified from fixed random range selection to fitted preprocessing in a **bounded range** (via gradient clipping). The preprocessing ranges are defined in the  the same from NNPDF3.1 and are defined in the `fitting:basis` parameter in the nnpdf runcard.


The old behaviour, in which the preprocessing is fixed randomly at the beginning of the fit, can be recovered by setting the `trainable` flag to false. See the [detailed runcard guide](runcard_detailed.html#preprocessing) for more information on how to define the preprocessing.


Optimizer
---------

In `n3fit` the genetic algorithm optimizer is replaced by modern stochastic gradient descent algorithms such as RMS propagation, Adam, Adagrad, among others provided by [Keras](https://keras.io/).

Following the gradient descent approach the training is performed in iteration steps where:
- for each data point the neural network is evaluated (forward propagation)
- the accumulated errors of each parameter is computed using the backward propagation algorithm, where starting from the analytical gradient of the loss function as a function of the neural network parameters the errors for each parameter is estimated.
- each parameter is updated accordingly to its weight, the gradient direction and the gradient descent update scheme (which controls the convergence step size and speed).

The gradient descent schemes are usually controlled by the **learning rate**, and the total **number of iterations**.
Examples of fits using the `n3fit` methodology are available here:
- DIS-only fit based on NNPDF3.1 NNLO setup, [view](https://vp.nnpdf.science/KTzrle5FQGuuBdcigkDKnQ==/)
- Global fit based on NNPDF3.1 NNLO setup: [view](https://vp.nnpdf.science/qtXzt-BbQZGkV6P4pf9-UA==/)

``` important:: The gradient descent scheme (RMSprop, Adagrad, etc.), the learning rate, the number of iteractions are hyper-parameters that require tuning.
```

Stopping algorithm
------------------

`n3fit` implements a patience algorithm which, together with the [positivity](#positivity) constraints, define when a fit is allowed to stop:

``` image:: figures/stopping.png
```

Following the diagram presented in the figure above, we then train the network until the validation stops improving. From that point onwards, and to avoid false positives, we enable a patience algorithm which waits for a number of iterations before actually considering the fit finished.

If the patience is set to a ratio 1.0 (i.e., wait until all epochs are finished) this algorithm is equal to that used in `nnfit`.

The look-back approach implemented in `nnfit` is not required by `n3fit` due to its less stochastic/random path towards the solution.

``` important:: The patience and the lagrange multipliers are hyper-parameters of the fit which require specific fine tuning.
```

Positivity
----------

In NNPDF3.1 there were a number of datasets added in order to constraint positivity based on DIS and fixed-target Drell-Yan processes. A similar technology and methodology is implemented in `n3fit` based on a penalty term controlled by a **positivity multiplier**.

The main difference to `nnfit` is that in `n3fit` a hard threshold is set such that no replicas generating negative values for the positivity sets are generated.

``` important:: The positivity multiplier is a hyper-parameter of the fit which require specific fine tuning.
```

Hyperoptimization algorithm
---------------------------

The main advantanges of the points presented above consists in the possibility to test several models in a fraction of time in comparison to the `nnfit` framework.

This is of key importance for a proper hyper-parameter scan where everything is potentially interconnected.

The hyperparameter scan capabilities are implemented using the hyperopt framework which systematically scans over a selection of parameter using Bayesian optimization and measures model performance to select the best architecture. The hyperopt library implements the tree-structured of Parzen estimator which is a robust sequential model based optimization approach ([SMBO](https://en.wikipedia.org/wiki/Hyperparameter_optimization)).

We optimize on a combination of the best validation loss and stability of the fits. In other words, we select the architecture which produces the lowest validation loss after we trim those combinations which are deemed to be unstable.

While performing the hyperparameter scan we found that optimizing only looking at the validation loss produced results which would usually be considered overfitted: very low training and validation chi2 and complex replica patterns. Thanks to the high performance of the `n3fit` procedure the usual cross-validation algorithm used in the NNPDF framework was not enough to prevent overlearning for all architectures.

The cross-validation implemented in NNPDF is successful on avoiding the learning of the noise within a dataset. However, we observe that this choice is not enough to prevent overfitting due to correlations within points in a same dataset when using hyperopt with `n3fit`.

In order to eliminate architectures that allowed overlearning we proceed by including a testing set where the model generalization power is tested. This is a set of datasets where none of the points are used in the fitting either for training or validation. Defining the best appropriate test dataset for PDF fits is particularly challenging due to the nature of the model regression through convolutions. For the present results the test set is defined by removing from the training setdatasets with duplicate process type and smaller leading-order kinematic range coverage. We call the loss produced by the removed datasets "testing loss" and we use it as a third criterion (beyond stability and combined with the value of the validation loss) to discard combinations of hyperparameters. With this procedure we are able to find combinations of hyperparameters which produce good fits for which we are confident no obvious overfitting can be generated.
