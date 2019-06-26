============================================================
**Documentation for theory covariance module in validphys2**
============================================================

:Author: Rosalyn Pearson

.. raw:: latex

   \maketitle

.. raw:: latex

   \tableofcontents

Summary
=======

-  The module of ``validphys2`` which deals with computation and
   interpretation of theoretical covariance matrices can be found in
   ``nnpdf/validphys2/src/validphys/theorycovariance/``, which consists
   of three files:

   #. ``construction.py``: deals with construction of covariance
      matrices and associated quantities

   #. ``output.py``: plots and tables

   #. ``tests.py``: actions for testing the covariance matrices against
      the NNLO-NLO shift

-  Theoretical covariance matrices are built according to the
   prescriptions found in Richard’s notes.

-  As input you need theories at the relevant scale combinations which
   correspond to the prescription.

-  These must be ordered in a specific way in the runcard.

-  The prescription is chosen based on the number of input theories,
   which must be one of :math:`\{3,5,7,9\}`.

-  In the case of 5 theories, you must further specify whether the 5 or
   :math:`\bar{5}` prescription is required. You can do this by
   allocating the flag ``fivetheories`` to ``nobar`` or ``bar`` in the
   runcard.

-  Currently the renormalisation scales are correlated within each
   process type. These process types are categorised as {DIS CC, DIS NC,
   Drell-Yan, Jets, Top}. However, this needs future consideration - if
   we add more experiments we need to be careful which process type they
   are allocated as.

-  Outputs include tables and heat plots of theoretical and combined
   (theoretical + experimental) covariance matrices, comparisons of
   theoretical and experimental errors, and plots and tables of
   :math:`\chi^2` values.

Important information about runcard layout
==========================================

Below is an example runcard (template in the next section) to calculate
the theoretical covariance matrix for 9-point variations at NLO and
display a range of outputs.

-  **IMPORTANT:** In runcards, theories must be listed according to
   points :math:`(\mu_0, \mu_i)` in the following order: :math:`(0,0)`,
   :math:`(+,0)`, :math:`(-,0)`, :math:`(0,+)`, :math:`(0,-)`,
   :math:`(+,+)`, :math:`(-,-)`, :math:`(+,-)`, :math:`(-,+)`. Here "+"
   refers to doubled scale, "-" to halved scale and "0" to central
   scale.

-  In terms of ``theoryids``, at NLO this corresponds to: 163, 177, 176,
   179, 174, 180, 173, 175, 178. This ensures that the prescriptions
   will be implemented correctly. The ``theoryids`` to input for each
   prescription are listed explicitly in the section on point
   prescriptions below. If an incorrect configuration of ``theoryids``
   is given, you should receive an error message.

-  The flag ``fivetheories`` specifies the choice of 5 or
   :math:`\bar{5}` prescription for the case of 5 input theories. You
   must assign a value ``nobar`` or ``bar`` correspondingly. If you do
   not do this, ``validphys`` will give an error.

-  As a fit you must provide the fit for central scales at the relevant
   order.

-  When calculating at NNLO, you must also provide c-factors.

-  Experiments will appear in tables and plots in the order they are
   listed in the runcard.

Outputs
=======

Below is the template corresponding to the runcard in the previous
section. This will produce a report with the important features of the
theory covariance module. These fall into three broad categories:

#. **Matrices and plots of matrices**. Heatmap plots, diagonal element
   plots comparing experimental and theoretical errors.

#. **:math:`\chi^2` values**. Total :math:`\chi^2`\ s and per dataset
   and per experiment. Bar chart comparing before and after adding
   theory errors for the per dataset case.

#. **Scale variations as a function of kinematics.** Plots of the
   theoretical predictions for each data point for the different scale
   choices.

Point prescriptions for theory covariance matrices
==================================================

The equations below display the different point prescriptions, as they
appear in ``validphys2``.

3 points
--------

theoryids: 163, 180, 173

.. math:: s_{11} = \frac{1}{2}\bigg\{ \Delta_1(+,+)^2 + \Delta_1(-,-)^2 \bigg\}

.. math:: s_{12} = \frac{1}{4}\bigg\{\bigg(\Delta_1(+,+) + \Delta_1(-,-) \bigg) \bigg(\Delta_2(+,+) + \Delta_2(-,-) \bigg) \bigg\}

.. _points-1:

 5 points
---------

theoryids: 163, 177, 176, 179, 174

.. math:: s_{11} = \frac{1}{2}\bigg\{ \Delta_1(+,0)^2 + \Delta_1(-,0)^2 + \Delta_1(0,+)^2 + \Delta_1(0,-)^2 \bigg\}

.. math::

   \begin{split}
       s_{12} = \frac{1}{2}\bigg\{ &\Delta_1(+,0)\Delta_2(+,0) + \Delta_1(-,0)\Delta_2(-,0) \bigg\} \\
               + \frac{1}{4}\bigg\{ &\bigg(\Delta_1(0,+) + \Delta_1(0,-) \bigg)\bigg(\Delta_2(0,+) + \Delta_2(0,-)\bigg)\bigg\}
   \end{split}

:math:`\mathbf{\overline{5}}` points
------------------------------------

theoryids: 163, 180, 173, 175, 178

.. math:: s_{11} = \frac{1}{2}\bigg\{ \Delta_1(+,+)^2 + \Delta_1(-,-)^2 + \Delta_1(+,-)^2 + \Delta_1(-,+)^2 \bigg\}

.. math::

   \begin{split}
       s_{12} = \frac{1}{4}\bigg\{ &\bigg(\Delta_1(+,+) + \Delta_1(+,-)\bigg) \bigg(\Delta_2(+,+) + \Delta_2(+,-) \bigg) \\
       + &\bigg(\Delta_1(-,+) + \Delta_1(-,-)\bigg) \bigg(\Delta_2(-,+) + \Delta_2(-,-) \bigg) \bigg\}
   \end{split}

7 points - original
-------------------

| Specify in the runcard ``seventheories: original``
| theoryids: 163, 177, 176, 179, 174, 180, 173

  .. math::

     \begin{split}
         s_{11} = \frac{1}{3}\bigg\{ &\Delta_1(+,0)^2 + \Delta_1(-,0)^2 + \Delta_1(0,+)^2 + \Delta_1(0,-)^2  \\                                 + &\Delta_1(+,+)^2 + \Delta_1(-,-)^2 \bigg\}
     \end{split}

.. math::

   \begin{split}
       s_{12} = \frac{1}{6}\bigg\{ &\bigg(\Delta_1(+,0) + \Delta_1(+,+) \bigg) \bigg(\Delta_2(+,0) + \Delta_2(+,+) \bigg) \\
               + &\bigg(\Delta_1(-,0)+\Delta_1(-,-)\bigg) \bigg(\Delta_2(-,0) + \Delta_2(-,-) \bigg) \\
               + &\bigg(\Delta_1(0,+)+\Delta_1(0,-)\bigg)\bigg(\Delta_2(0,+) + \Delta_2(0,-) \bigg)\bigg\}
   \end{split}

7 points - Gavin (default)
--------------------------

theoryids: 163, 177, 176, 179, 174, 180, 173

.. math::

   \begin{split}
       s_{11} = \frac{1}{3}\bigg\{ &\Delta_1(+,0)^2 + \Delta_1(-,0)^2 + \Delta_1(0,+)^2 + \Delta_1(0,-)^2  \\                                 + &\Delta_1(+,+)^2 + \Delta_1(-,-)^2 \bigg\}
   \end{split}

.. math::

   \begin{split}
       s_{12} = \frac{1}{6}\bigg\{ &2\bigg(\Delta_1(+,0)\Delta_2(+,0) + \Delta_1(-,0)\Delta_2(-,0) \bigg) \\
               + &\bigg(\Delta_1(0,+)+\Delta_1(0,-)\bigg) \bigg(\Delta_2(0,+) + \Delta_2(0,-) \bigg) \\
               + &\bigg(\Delta_1(+,+)+\Delta_1(-,-)\bigg)\bigg(\Delta_2(+,+) + \Delta_2(-,-) \bigg)\bigg\}
   \end{split}

.. _points-2:

9 points
--------

theoryids: 163, 177, 176, 179, 174, 180, 173, 175, 178

.. math::

   \begin{split}
       s_{11} = \frac{1}{4}\bigg\{ &\Delta_1(+,0)^2 + \Delta_1(-,0)^2
                               + \Delta_1(0,+)^2 + \Delta_1(0,-)^2 \\
                               + &\Delta_1(+,+)^2 + \Delta_1(+,-)^2 
                               + \Delta_1(-,+)^2 + \Delta_1(-,-)^2 \bigg\}
   \end{split}

.. math::

   \begin{split}
       s_{12} = \frac{1}{12}\bigg\{&\bigg(\Delta_1(+,0)+\Delta_1(+,+) + \Delta_1(+,-)\bigg) \bigg(\Delta_2(+,0) + \Delta_2(+,+) + \Delta_2(+,-) \bigg) \\
               + &\bigg(\Delta_1(-,0) + \Delta_1(-,+) + \Delta_1(-,-)\bigg)\bigg(\Delta_2(-,0) + \Delta_2(-,+) + \Delta_2(-,-) \bigg) \bigg\}\\
               + \frac{1}{8}&\bigg(\Delta_1(0,+)+ \Delta_1(0,-)\bigg)\bigg(\Delta_2(0,+) + \Delta_2(0,-) \bigg)
   \end{split}

Tests
=====

.. _projection:

Projection onto the subspace of nonzero eigenvalues
---------------------------------------------------

**Notation overview:**

-  :math:`D` is the space of data, indexed by
   :math:`i,j = 1,...,N_{dat}`

-  :math:`E \in D` is the subspace of non-zero eigenvalues, indexed by
   :math:`\alpha, \beta = 1,..., N_{eval}`

-  non-dashed variables are in :math:`D`, dashed variables are in
   :math:`E`

The covariance matrix :math:`S` can be expressed as a sum of outer
products of vectors :math:`X_{\alpha}`, generated by taking the
difference between observables at the central scale and at a given
varied scale:

.. math::

   \label{covmat}
       S = \sum_{\alpha, \beta} \lambda_{\alpha \beta}X_{\alpha}X_{\beta}^T.

 :math:`S` satisfies the eigenvalue equation :math:`Sv = wv` for
:math:`N_{eval}` non-zero eigenvalues. :math:`E` is spanned by
:math:`\{X_{\alpha}\}`. Generate an orthonormal basis
:math:`\{Y_{\alpha}\}` from :math:`\{X_{\alpha}\}` by the Gram-Schmidt
process:

.. math:: X_n^\prime = X_n - \sum_{m \leq n} (Y_m.X_n)Y_m,

 normalising to give :math:`Y_n = X^\prime_n/|X^\prime_n|`. Then define

.. math::

   P =
     \begin{bmatrix}
       Y_1 & ... & Y_{N_{eval}} \\
       \downarrow & ... & \downarrow 
     \end{bmatrix}

 and project :math:`S` into :math:`E` using
:math:`S^\prime_{\alpha \beta} = P^{-1}_{i \alpha} S_{ij} P_{j \beta}`.
Solve to find the eigenvalues :math:`w_{\alpha}` and their corresponding
eigenvectors :math:`v^\prime_{\alpha}`.

Consider

.. math::

   \begin{split}
           S v_\alpha &= w_\alpha v_\alpha \\
           PS^\prime P^{-1} v_\alpha &= w_\alpha v_\alpha \\
           P^{-1}PS^\prime P^{-1} v_\alpha &= P^{-1}w_\alpha v_\alpha \\
           S^\prime (P^{-1} v_\alpha) &= w_\alpha (P^{-1}v_\alpha).
       \end{split}

 So :math:`v^\prime_\alpha = P^{-1}v_\alpha`, and the eigenvectors in
:math:`D` can be found using

.. math::

   \begin{split}
       (v_{\alpha})_i &= P_{i \beta}(v^\prime_\alpha)_\beta \\
                      &= \sum_\beta (v^\prime_\alpha)_\beta Y_{i\beta}.
   \end{split}

Code implementation
-------------------

The code for testing theory covariance matrices against the observed
NNLO-NLO shift, :math:`\delta`, is contained in ``tests.py``. In order
to compare these, we need to ensure that cuts are matched between the
NNPDF3.1 theories 52 and 53 (NLO and NNLO respectively), and the
scale-varied theories.

The action ``evals_nonzero_basis`` takes the matched theory covariance
matrix and projects it from the data space into the basis of non-zero
eigenvalues, dependent on point prescription. It then returns the
eigenvalues and the data-space eigenvectors. These are taken as inputs
by ``theory_shift_test``, which compares them with the NNLO-NLO shift to
assess the missing fraction :math:`\delta_{miss}`, ``fmiss``, of the
shift vector which is not covered by the theory covariance matrix, and
the projections of the shift vector onto each of the eigenvectors
(``projectors``). The various outputs are:

#. ``theory_covmat_eigenvalues``: returns a table of
   :math:`s = \sqrt{eval}`, the projector and the ratio of the two,
   ordered by largest eigenvalue

#. ``efficiency``: returns the efficiency with which the theory
   covariance matrix encapsulates the NNLO-NLO shift,
   :math:`\epsilon = 1-\frac{\delta_{miss}}{\delta}`.

#. ``validation_theory_chi2``: returns the theory :math:`\chi^2`,
   defined as
   :math:`\frac{1}{N_{eval}}\sum_a \bigg(\frac{\delta_a}{s_a}\bigg)^2`.

#. ``projector_eigenvalue_ratio``: produces a plot of the ratio between
   the projectors and the square roots of the corresponding eigenvalues.

#. ``shift_diag_cov_comparison``: produces a plot of the NLO-NNLO shift
   compared to an envelope given by the diagonal elements of the theory
   covariance matrix.

``evals_nonzero_basis``
~~~~~~~~~~~~~~~~~~~~~~~

This section details the way that ``evals_nonzero_basis`` determines the
non-zero eigenvalues and eigenvectors.

First we construct the differences between the theory results for
shifted scales and for the central scale, ``diffs``. These are ordered
by dataset, so processes are jumbled up. We split these up into a series
of vectors, ``splitdiffs``, of which there are :math:`p` for each of the
``diffs``. They each correspond to the values for one process from one
of the ``diffs``, and zeroes everywhere else, but where the entries are
ordered such that each process occupies a different space in the vector,
and is gathered together. i.e.

.. math::

   \Delta(+;+) \to 
       \begin{pmatrix}
       \Delta^1(+;+) \\
       0 \\
       ...\\
       0
       \end{pmatrix},
       \begin{pmatrix}
       0 \\
       \Delta^2(+;+) \\
       ... \\
       0
       \end{pmatrix}, ...,
       \begin{pmatrix}
       0 \\
       0 \\
       ... \\
       \Delta^p(+;+)
       \end{pmatrix},

 where the superscript labels the process number. This allows us to
decorrelate scale variations between different processes.

The next stage is to use the ``splitdiffs`` to construct the linearly
independent vectors corresponding to each prescription, which are those
detailed in Sec. `5.3 <#vectors>`__. Each of these is implemented in a
separate function called ``vectors_#pt``, where # is the number of
points in the prescription, which is called as appropriate. The next
stage is to orthonormalise these vectors (labelled ``xs``) using the
Gram-Schmidt method to produce ``ys``. The covariance matrix is then
projected into the space of these, and the eigenvalues and eigenvectors
are found. The eigenvectors in the dataspace are then calculated by
rotating back, as detailed in Sec. `5.1 <#projection>`__.

.. _vectors:

Linearly independent vectors for each prescription
--------------------------------------------------

3 point
~~~~~~~

For :math:`(\mu_1, \mu_2, ..., \mu_p)`, where :math:`\mu_0` is
correlated with :math:`\mu_i` variation for each process.
:math:`(+, +, +, ...)` :math:`(-, +, +, ...)` + cyclic :math:`p+1`
vectors

.. _point-1:

5 point
~~~~~~~

For :math:`(\mu_0; \mu_1, \mu_2, ..., \mu_p)`.
:math:`(\pm; 0, 0, 0, ...)` :math:`(0; +, +, +, ...)`
:math:`(0; -, +, +, ...)` + cyclic :math:`p+3` vectors

:math:`\bar{5}` point
~~~~~~~~~~~~~~~~~~~~~

For :math:`(\mu_0; \mu_1, \mu_2, ..., \mu_p)`.
:math:`(\pm; +, +, +, ...)` :math:`(\pm; -, +, +, ...)` + cyclic
:math:`2p+2` vectors

.. _point-2:

7 point
~~~~~~~

Combine 3 and 5 point :math:`2p+4` vectors

.. _point-3:

9 point
~~~~~~~

For :math:`(\mu_0; \mu_1, \mu_2, ..., \mu_p)`.
:math:`(\pm; +, +, +, ...)` :math:`(0; +, +, +, ...)`
:math:`(\pm; -, +, +, ...)` + cyclic :math:`(0; -, +, +, ...)` + cyclic
:math:`(\pm; 0, +, +, ...)` + cyclic :math:`5p+3` vectors
