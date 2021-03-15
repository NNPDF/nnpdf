 .. _vptheorycov-tests:
 
Tests
=====

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

In this section we summarise for reference the linearly independent vectors
contributing to the covariance matrix in each prescription.

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
