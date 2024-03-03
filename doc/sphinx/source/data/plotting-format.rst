.. _plotting-format:

===============
Plotting format
===============

The ``plotting`` dictionary within the metadata of a dataset
defines a set of options that are used for analysis
and representation purposes, particularly to determine how datasets
should be represented in plots.

.. warning:: the information in this page is not up to date

Format
======

The plotting file specifies the variable in which the data
is to be plotted (in the  *x* axis) as well as the variables
in which the data will be split in different lines in the
same figure or in different figures. The possible variables
('*kinematic labels*') are described below.

The format also allows the control of several plotting properties, such
as whether to use log scale, or the axes labels.

Kinematic labels
================

.. note:: very outdated information that only applies to legacy data

When a dataset has been ported from the old implementation and thus
it has no well defined kinematic variables (but instead just k1, k2, k3)
the default kinematic variables are inferred from the *process type*
declared in the commondata files (more specifically from
a substring). Currently they are:

.. code-block:: python

  'DIS': ('$x$', '$Q^2 (GeV^2)$', '$y$'),
  'DYP': ('$y$', '$M^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
  'EWJ_JPT': ('$p_T (GeV)$', '$M^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
  'EWJ_JRAP': ('$\\eta/y$', '$M^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
  'EWJ_MLL': ('$M_{ll} (GeV)$', '$M_{ll}^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
  'EWJ_PT': ('$p_T (GeV)$', '$M^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
  'EWJ_PTRAP': ('$\\eta/y$', '$p_T^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
  'EWJ_RAP': ('$\\eta/y$', '$M^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
  'EWK_MLL': ('$M_{ll} (GeV)$', '$M_{ll}^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
  'EWK_PT': ('$p_T$ (GeV)', '$M^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
  'EWK_PTRAP': ('$\\eta/y$', '$p_T^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
  'EWK_RAP': ('$\\eta/y$', '$M^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
  'HIG_RAP': ('$y$', '$M_H^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
  'HQP_MQQ': ('$M^{QQ} (GeV)$', '$\\mu^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
  'HQP_PTQ': ('$p_T^Q (GeV)$', '$\\mu^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
  'HQP_PTQQ': ('$p_T^{QQ} (GeV)$', '$\\mu^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
  'HQP_YQ': ('$y^Q$', '$\\mu^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
  'HQP_YQQ': ('$y^{QQ} (GeV)$', '$\\mu^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
  'INC': ('$0$', '$\\mu^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
  'JET': ('$\\eta$', '$p_T^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
  'PHT': ('$\\eta_\\gamma$', '$E_{T,\\gamma}^2 (GeV^2)$', '$\\sqrt{s} (GeV)$'),
  'SIA': ('$z$', '$Q^2 (GeV^2)$', '$y$')


The three kinematic variables are referred to as `k1`, `k2` and `k3`
in the plotting files. For example, for DIS processes, `k1` refers to `x`,
`k2` to `Q`, and `k3` to `y`.

These kinematic values can be overridden by some transformation of
them. For that purpose, it is possible to define
a `kinematics_override` key.  The value must be a class defined
in: `validphys2/src/validphys/plotoptions/kintransforms.py`

The class must have a `__call__` method that takes three parameters:
`(k1, k2 k3)` as defined in the dataset implementation, and returns
three new values `('k1', 'k2', k3')` which are the "transformed"
kinematical variables, which will be used for plotting purposes every
time the kinematic variables `k1`, `k2` and `k3` are referred to.
Additionally, the class must implement a `new_labels` method, that
takes the old labels and returns the new ones, and an `xq2map`
function that takes the kinematic variables and returns a tuple of (x,
Q²) with some approximate values. An example of such transform is:

.. code-block:: python

  class dis_sqrt_scale:
      def __call__(self, k1, k2, k3):
          ecm = sqrt(k2/(k1*k3))
          return k1, sqrt(k2), ceil(ecm)

      def new_labels(self, *old_labels):
          return ('$x$', '$Q$ (GeV)', r'$\sqrt{s} (GeV)$')

      def xq2map(self, k1, k2, k3, **extra_labels):
          return k1, k2*k2


Additional labels
=================
Additional labels can be specified by declaring an **extra_labels**
key in the plotting file, and specifying for each new label a value
for each point in the dataset.

For example:

..  code-block:: yaml
  
   extra_labels:
      idat2bin:  [0, 0, 0, 0, 0, 0, 0, 0, 100, 100, 100, 100, 100, 200, 200, 200, 300, 300, 300, 400, 400, 400, 500, 500, 600, 600, 700, 700, 800, 800, 900, 1000, 1000, 1100]

defines one label where the values for each of the datapoints are
given in the list. Note that the name of the extra_label (in this case
`idat2bin` is completely arbitrary, and will be used for plotting
purposes (LaTeX math syntax is allowed as well). However, adding labels
manually for each point can be tedious. This should only be reserved
for information that cannot be recovered from the kinematics as
defined in the CommonData file. Instead, new labels can be generated
programmatically: every function defined in `validphys2/src/validphys/plotoptions/labelers.py`
is a valid label. These functions take as keyword arguments the
(possibly transformed) kinematical variables, as well as any extra
label declared in the plotting file. For example, one might declare:

.. code-block:: python

  def high_xq(k1, k2, k3, **kwargs):
      return k1 > 1e-2 and k2 > 1000


Note that it is convenient to always declare the `**kwargs`
parameter so that the code doesn't crash when the function is called
with extra arguments. Similarly to the kinematics transforms, it is
possible to decorate them with a `@label` describing a nicer latex
label than the function name. For example:

.. code-block:: python

  @label(r"$I(x>10^{-2})\times I(Q > 1000 GeV)$")
  def high_xq(k1, k2, k3, **kwargs):
      return (k1 > 1e-2) & (k2 > 1000)


Plotting and grouping
=====================

The variable in which the data is plotted is simply
declared as

..  code-block:: yaml
    
  x: <label>

For example:

..  code-block:: yaml
    
  x: k1

If a `line_by` key is specified, variables with different values for
each of the labels listed, will be represented as different lines. For
example,


..  code-block:: yaml
    
  line_by:
    - k2

for DIS would mean that the data in the same Q bin is plotted in the
same line.

Similarly, it is possible to define a `figure_by` key: Points
with different values for the listed keys will be split across
separated figures. For example:

..  code-block:: yaml
    
  figure_by:
    - idat2bin
    - high_xq


Transforming the result
=======================

.. note:: very outdated information that only applies to legacy data

By default the *y* axis represents the central value and error. However,
it is possible to define a results_transform in the plotting file:

..  code-block:: yaml
    
  result_transform: qbinexp

The value must be a function declared in
`validphys2/src/validphys/plotoptions/results_transform.py`
taking the error, the central value, as well as all the labels, and
returning a new error and central value. For example:

..  code-block:: python
    
  def qbinexp(cv, error, **labels):
      q = labels['k2']
      qbin = bins(q)
      return 10**qbin*cv, 10**qbin*error

Plotting options
================

Several plotting options can be specified.
These include

 - x/y_scale: 'linear' or 'log'.
 - x/y_label: Any string, possibly latex formatted. Note that the
	 x_label will be deduced automatically.

Overriding configuration for normalized plots
=============================================

When the results are to be plotted as a ratio, it may be convenient to
alter the configuration of the plots, for example by changing the
`line_by` labels into `figure_by` (because otherwise the points would
overlap), or by changing the scale from log to linear. To do so, we
specify the options we want to override in a `normalize` key.
Everything defined inside will take precedence when we produce a ratio
plot and will be ignored for absolute value plots. For example:

..  code-block:: yaml
    
  x: k1

  x_label: '$\left\|\eta/y\right|$'

  y_label: '$d\sigma/dy$ (fb)'

  line_by:
    - Boson

  normalize:
      figure_by:
          - Boson

  extra_labels:
    Boson:  ["$W^+$","$W^+$","$W^+$","$W^+$","$W^+$","$W^+$","$W^+$","$W^+$","$W^+$","$W^+$","$W^+$","$W^-$","$W^-$","$W^-$","$W^-$","$W^-$","$W^-$","$W^-$","$W^-$","$W^-$","$W^-$","$W^-$","$Z$","$Z$","$Z$","$Z$","$Z$","$Z$","$Z$","$Z$"]

Here, we would split the data by different figure files for each
unique value of the key `Boson` (which is defined explicitly as an
`extra_label`), but only one plot with the three bosons split across
different lines will be produced in absolute value plots.

Example
=======

..  code-block:: yaml

  dataset_label: "Some hypothetical dataset"
  x: k3
  x_scale: log
  kinematics_override: dummy_transform #defined in transforms.py
  line_by:
    - k2

  figure_by:
    - idat2bin #defined below
    - high_xq  #defined in labelers.py

  normalize: # Change the scale for ratio plots
      x_scale: linear

  extra_labels:
      idat2bin:  [0, 0, 0, 0, 0, 0, 0, 0, 100, 100, 100, 100, 100, 200, 200, 200, 300, 300, 300, 400, 400, 400, 500, 500, 600, 600, 700, 700, 800, 800, 900, 1000, 1000, 1100]
