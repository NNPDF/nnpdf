.. _fktables:

============================================================
Fast Kernel (FK) tables
============================================================

.. raw:: latex

   \maketitle

.. raw:: latex

   \tableofcontents

Here we discuss the numerical implementation of the calculations of the DIS structure functions.

In the framework of collinear QCD factorization, the :math:`F_2` structure function
can be decomposed in terms of hard-scattering coefficient functions and PDFs as

.. math::

    F_2(x,Q^2) = \sum_i^{n_f} C_i(x,Q^2) \otimes f_i(x,Q^2),

where :math:`C_i(x,Q^2)` are the process-dependent coefficient functions which
can be computed perturbatively as an expansion in the QCD and QED
couplings, :math:`f_i(x,Q^2)` are the PDFs at the scale :math:`Q^2`,
and :math:`\otimes` denotes the Mellin convolution.
The sum over flavours runs over the :math:`n_f` active quark flavours and
their antiquarks at the relevant scale, as well as over the gluon.

The PDFs at :math:`Q^2` are obtained from those at the initial parameterization
scale :math:`Q_0^2` using the evolution operator
:math:`E_{ij}(Q^2 \leftarrow Q_0^2)`, determined by the DGLAP equations.
For a detailed discussion, see `Solving DGLAP in the EKO documentation
<https://eko.readthedocs.io/en/latest/theory/DGLAP.html>`_.

The direct calculation of the above equation during the PDF fit is not practical
since it requires first solving the DGLAP evolution equation for each new boundary
condition at :math:`Q_0` and then convoluting with the coefficient
functions.

Constructing FK tables
----------------------

To evaluate the observable efficiently, we precompute the perturbative
coefficient functions :math:`C_i`, using `Yadism
<https://yadism.readthedocs.io/en/v0.13.8/overview/tutorials/getting_started.html#dumping-predictions-into-fast-interpolation-grids>`_
for DIS, and the evolution operators :math:`E_{ij}`, using EKO.
We choose a suitable interpolation basis as discussed in the `EKO documentation
<https://eko.readthedocs.io/en/latest/theory/Interpolation.html>`_.

`Pineko <https://github.com/NNPDF/pineko>`_ combines the coefficient-function
grids in `PineAPPL <https://pineappl.readthedocs.io/en/latest/>`_ format with
the EKO evolution operators to produce Fast Kernel (FK) tables.
These encode the perturbative coefficients and PDF evolution, leaving only the
PDFs at the input scale to be supplied during the fit:

.. math::

    \begin{aligned}
    F_2(x,Q^2)
    &= \sum_{i,j}^{n_f} C_i(x,Q^2) \otimes E_{ij}(Q^2 \leftarrow Q_0^2)
       \otimes f_j(x,Q_0^2) \\
    &= \sum_j^{n_f} \sum_{\alpha=1}^{n_x}
       FK_{j,\alpha}(x,x_{\alpha},Q^2,Q_0^2)\,
       f_j(x_{\alpha},Q_0^2).
    \end{aligned}

Here :math:`j` runs over the input parton flavours and :math:`\alpha` labels
the :math:`n_x` momentum-fraction grid points at the input scale.

The convolutions are thus reduced to a matrix multiplication of the precomputed FK table
with the PDFs evaluated at the input scale. This avoids repeating the
perturbative calculation for each trial PDF during the fit, increasing the numerical calculation speed of the 
DIS structure functions by up to several orders of magnitude.

.. note::

   Hadronic processes (proton-proton collisions), such as Drell-Yan production,
   follow the same logic. Their FK tables have a slightly more complicated
   structure to account for the evolution of two PDFs instead of one.
