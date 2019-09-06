Stopping algorithm
=============================

`n3fit` implements a patience algorithm which, together with the [positivity](./positivity) constraint define when a fit is allowed to stop.

(add picture from the paper)

If the patience is set to a ratio 1.0 (i.e., wait until all epochs are finished) this algorithm is equal to that used in `nnfit`.

(links from reports)
