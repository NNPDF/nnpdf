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

Optimizer
---------

What is Gradient Descent
- Forward Propagation
- Backward Propagation
- Deterministic

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