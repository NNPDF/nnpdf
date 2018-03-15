%TODOS for validphys

Urgent for gallery
==================

Kinematics plot
---------------

Way to recover the coordinates from the PLOTTING specification.
Hopefully some `f(k1,k2,k3, process_type)` after the kinematics are
transformed should suffice.

Otherwise, we would add `kin_x`, `kin_q` to the plotting specification
for nondefault cases.


Luminosity
----------

Just need to hook up the proper calls to  `scipy.integrate` and logic
to sum over flavors. Should be based on an already existing code.

From there, we can make the usual lumi plots, as well as the Karlsruhe
plot.

SMPDF plots
----------

Done. Needs to be compared to the SMPDF. Some minor improvements could
be:

 - Put the colorbar at the bottom.
 - Reduce the space between the title and the axis.
 - In HTML, why is it so narrow? The width of the file is the same as
   for the rest.

Obs-Obs correlations
--------------------

Need to decide how to group them We don't want a 4000x4000 table.
Maybe by default take `percentile(abs(rho), 90)` or similar.

Distance plots
--------------

Trivial to do the PDF distances. Do we want a better measure of
distance such as Kolmogorov?


Functional improvements
=======================

Smart tools
-----------

 - Implement filters for pdf replicas, data points, datasets,
	 experiments.
º	
 - Make the relevant plots aware of the filters and teach them how to
	 display the filtered points.

The UI for using them should look like:

```yaml
replica_filters:
   - bad_dataset:
       threshold: 2

dataset_input:
   dataset: CMSDY2D12

theoryid: 65

pdf: NNPDF30_nlo_as_0118

fit: myfit

use_cuts:True

actions_:
  - - plot_pdfreplicas
    - plot_training_validation
```
and  what this does is mark all the replicas that are above the
threshold in all the plots with replicas.

Similar to filters, there should be an interface for colors, so that
we can e.g. do the kinematic plot with a custom color scale.

HTML+CSS (+Python)
--------

Important missing feature? Decide after some more experience with
reports produced with the existing features. Maybe more functionality
to group figures.


Important for completeness
==========================

Basis transformations
----------------------

  > To late to hook with C++?

Otherwise implement the same logic in Python, which would surely be
faster.

Training statistics
-------------------

Length and so on.

Statistical estimators
----------------------

$\phi$ and so on.
