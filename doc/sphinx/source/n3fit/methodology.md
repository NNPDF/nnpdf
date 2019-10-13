Methodology overview
====================

In this document we summarise from a conceptual point of view the main points which are different in comparison to the latest NNPDF (i.e. [NNPDF3.1](https://arxiv.org/abs/1706.00428)) methodology.

This document contains a more specific discussion about the choices currently implemented in the `n3fit` package and discussed for the first time in [hep-ph/1907.05075](https://arxiv.org/abs/1907.05075).

**Table of contents:**
- [Introduction](#introduction)
- [Neural network architecture](#neural-network-architecture)
- [Optimizer](#optimizer)
- [Stopping](#stopping-algorithm)
- [Positivity](#positivity)
- [Hyperoptimization](#hyperoptimization-algorithm)

``` warning:: The final setup used in n3fit fits can be extracted from the runcards stored in nnpdf/n3fit/runcards.
```

Introduction
------------

The approach presented here inherits the technology develop by the NNPDF collaboration in terms of fit pipeline but extends the possibility to test and improve fitting performance with modern techniques inspired by the deep learning community.

The `n3fit` code is designed in python and one of its main goals is to replace the `nnfit` program. It provides a simple abstraction layer which simplifies the live of developers when considering the possibility to add new fitting algorithms.

In the following table we list some of the differences between both codes:

```eval_rst
+--------------------+---------------------------------+-----------------------------------+
| Component          | nnfit                           | n3fit                             |
+====================+=================================+===================================+
| MC data generation | single seed                     | multi seed                        |
+--------------------+---------------------------------+-----------------------------------+
| Data management    | libnnpdf                        | same as nnfit                     |
+--------------------+---------------------------------+-----------------------------------+
| Neural net         | fixed architecture, per flavour | single net, flexible architecture |
+--------------------+---------------------------------+-----------------------------------+
| Preprocessing      | random fixed                    | fitted in range                   |
+--------------------+---------------------------------+-----------------------------------+
| Integration        | a posteriori per iteration      | buildin in the model              |
+--------------------+---------------------------------+-----------------------------------+
| Optimizer          | genetic optimizer               | gradient descent                  |
+--------------------+---------------------------------+-----------------------------------+
| Stopping           | lookback                        | patience                          |
+--------------------+---------------------------------+-----------------------------------+
| Positivity         | fixed lagrange multiplier       | dynamic multiplier                |
+--------------------+---------------------------------+-----------------------------------+
| Postfit            | 4-sigma chi2 and arclenght      | same as nnfit                     |
+--------------------+---------------------------------+-----------------------------------+
| Fine tuning        | manual                          | automatic                         |
+--------------------+---------------------------------+-----------------------------------+
```


``` note:: In the next sections we will focus on n3fit specifics.
```

Neural network architecture
---------------------------

The main advantage of using a modern deep learning backend such as keras/tensorflow consists in the possibility to change neural network architecture quickly because the developer can avoid fine tuning the code in order to achieve efficient memory management and PDF convolution performance.

The current `n3fit` code supports sequential dense networks with custom number of layers, nodes, activation functions and initializers.

A big difference in comparison to `nnfit` is the number of neural networks involved in the fit. Here we use a **single neural network** model which maps (x, log x) in 8 outputs, nominally exactly the 8 PDF flavours defined in NNPDF3.1. The choice of using a single model is justified by the interest in preserving cross-correlations between flavours.

Preprocessing has been modified from fixed random range selection to fitted preprocessing in a **bounded range** (via gradient clipping). The preprocessing ranges are the same from NNPDF3.1.

The momentum sum rules are implemented as a **neural network layer** which computes the normalization coefficients for each flavour as a sum over a fixed grid of x points. The number and density of points in x is select in such way that the final quality of the integrals are at least permille level in comparison to 1D integration algorithms.

Optimizer
---------

In `n3fit` the genetic optimizer is replaced by modern stochastic gradient descent algorithms such as RMS propagation, Adam, Adagrad, etc.

In the gradient descent approach the training is performed in iteration where:
- for each data point the neural network is evaluated (forward propagation)
- the accumulated errors of each parameter is computed using the backward propagation algorithm, where starting from the analytical gradient of the loss function as a function of the neural network parameters the errors for each parameter is estimated.
- each parameter is updated accordingly to its weight, the gradient direction and the select scheme (which controls the convergence step size and speed).

Effects in NNPDF
(links to reports in the wiki)


Stopping algorithm
------------------

`n3fit` implements a patience algorithm which, together with the [positivity](./positivity) constraint define when a fit is allowed to stop.

(add picture from the paper)

If the patience is set to a ratio 1.0 (i.e., wait until all epochs are finished) this algorithm is equal to that used in `nnfit`.

(links from reports)


Positivity
----------

In NNPDF 3.1 there were a number of datasets added in order to constraint positivity.

In `n3fit` instead a hard threshold is set such that no replicas generating negative values for the positivity sets are generated.

(link to reports)

Hyperoptimization algorithm
---------------------------

The main advantanges of the points presented above consists in the possibility to test several models in a fraction of time in comparison to the `nnfit` framework.

The  NNPDF  approach  aimed  to  reduce  the  bias  intro-duced in the determination of the functional form utilizedto  parametrize  the  PDFs  [22].  Neural  networks  provideuniversal  function  approximators  [23]  which  reduce  sys-tematic biases introduced by the choice of specific func-tional  forms.  Neural  networks  themselves,  however,  arenot unique and the space of hyperparameters is big enoughthat finding the best choice becomes a overwhelming task.In this work we aim to improve over the previous it-eration of the NNPDF methodology by boxing the entireframework  under  hyperparameter  optimization  routines,there are several key points which allow us to do this now.Firstly, the new design of the code exposes all parametersof the fit including (but not restricted to) the neural net-work architecture. This is of key importance for a properhyperparameter  scan  where  everything  is  potentially  in-terconnected. Secondly, the new methodology has such asmaller impact on computing resources that we can per-form more fits on a scale of orders of magnitude, in otherwords, for each fit using the old methodology we can nowtest hundreds of architectures.The hyperparameter scan capabilities are implementedusing  the  hyperopt  framework  [17]  which  systematicallyscans over a selection of parameter using Bayesian opti-mization [24] and measures model performance to selectthe best architecture.As a proof of concept, for this paper we make a firstselection of parameters on which to scan, shown in Table 2separated between the parameters which define the NeuralNetwork  architecture  and  those  which  define  the  fittingprocedure.In  this  study  we  apply  the  framework  to  both  theglobal  and  DIS  only  setup  and  in  order  to  achieve  thebest model configuration we limit the data input to theexperimental central values instead of using artificial repli-cas. We optimize on a combination of the best validationloss and stability of the fits. In other words, we select thearchitecture which produces the lowest validation loss af-ter we trim those combinations which are deemed to beunstable.In Fig. 3 we show an example of DIS only scan. Wepresent eight examples of those shown in Table 2In this scan we find the Adadelta optimizer, for whichno  learning  rate  is  used,  to  be  more  stable  and  system-atically produce better results than RMSprop and Adamwith a wide choice of learning rates. The initializers, onceunstable options such as a random uniform initializationhave been removed, seem to provide similar qualities witha  slight  preference  for  the  “glorotnormal”  initializationprocedure described in [25].Concerning the parameters related to the stopping cri-teria, we observe that when the number of epochs is verysmall the fit can be certainly unstable, however after a cer-tain threshold no big differences are observed. The stop-ping patience shows a very similar pattern, stopping tooearly  can  be  disadvantageous  but  stopping  to  late  doesseem to make a big difference. The positivity multiplier,however, shows a clear preference for bigger numbers.Finally, concerning the neural network architecture weobserve  that  a  small  number  of  layers  seem  to  produceslightly  better  absolute  results,  however,  one  single  hid-den layer seem to be also very inconsistent. The activationfunctions present with a slight preference for the hyper-bolic tangent. Once we have a acceptable hyperparametersetup we ran again for fine tuning as some of the choicescould have been biased by a bad combination on the otherparameters.The main take away from this scan is the implemen-tation of a semi-automatic methodology able to find thebest  hyperparameter  combination  as  the  setup  changes,e.g.with new experimental data, new algorithms or tech-nologies.

While performing the hyperparameter scan we found thatoptimizing  only  looking  at  the  validation  loss  producedresults which would usually be considered overfitted: verylow training and validationχ2and complex replica pat-terns. Thanks to the high performance of then3fitpro-cedure  the  usual  cross-validation  algorithm  used  in  theNNPDF framework was not enough to prevent overlearn-ing for all architectures.

The  cross-validation  implemented  in  NNPDF  is  suc-cessful on avoiding the learning of the noise within a dataset.However, we observe that this choice is not enough to pre-vent overfitting due to correlations within points in a samedataset when using hyperopt withn3fit. In order to elim-inate architectures that allowed overlearning we proceedby including a testing set where the model generalizationpower  is  tested.  This  is  a  set  of  datasets  where  none  ofthe  points  are  used  in  the  fitting  either  for  training  orvalidation.Defining the best appropriate test dataset for PDF fitsis particularly challenging due to the nature of the modelregression  through  convolutions.  For  the  present  resultsthe test set is defined by removing from the training setdatasets with duplicate process type and smaller leading-order kinematic range coverage. We call the loss producedby the removed datasets “testing loss” and we use it asa third criterion (beyond stability and combined with thevalue  of  the  validation  loss)  to  discard  combinations  ofhyperparameters. With this procedure we are able to findcombinations of hyperparameters which produce good fitsfor which we are confident no obvious overfitting can begenerated. In Table 3 we list the datasets which have beenused as test set for this study.In Fig. 4 we show an example of a PDF produced bytwo  very  different  architectures,  both  of  which  are  gen-erated by the hyperoptimization procedure. We observe amuch more unstable behaviour in the fit in which we allowfor overtraining which in turn translates for aχ2on thetesting set of more than twice the value of the training set.We believe the issue of hidden correlations in experimentalmeasurements as well as its impact on PDF fits requiresa much deeper study outside the scope of this paper.