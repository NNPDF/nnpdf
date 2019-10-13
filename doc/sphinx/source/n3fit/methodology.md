Methodology overview
====================

Optimizer
---------

## What is Genetic Algorithm
(reference to somewhere else in this documentation I HOPE)

## What is Gradient Descent
- Forward Propagation
- Backward Propagation
- Deterministic

## Effects in NNPDF
(links to reports in the wiki)


Stopping algorithm
---------

`n3fit` implements a patience algorithm which, together with the [positivity](./positivity) constraint define when a fit is allowed to stop.

(add picture from the paper)

If the patience is set to a ratio 1.0 (i.e., wait until all epochs are finished) this algorithm is equal to that used in `nnfit`.

(links from reports)


Positivity
---------

In NNPDF 3.1 there were a number of datasets added in order to constraint positivity.

In `n3fit` instead a hard threshold is set such that no replicas generating negative values for the positivity sets are generated.

(link to reports)

Hyperoptimization algorithm
---------

Differences with nnfit
---------