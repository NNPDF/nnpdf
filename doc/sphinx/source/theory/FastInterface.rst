.. _fktables:

============================================================
Fast Interface (FK tables)
============================================================

.. raw:: latex

   \maketitle

.. raw:: latex

   \tableofcontents

Here we discuss the numerical implementation of the calculations of the DIS structure functions.

In the framework of collinear QCD factorization, the :math:`F_2` structure function
can be decomposed in terms of hard-scattering coefficient functions and PDFs as,

.. math::

    \begin{align} 
    \label{eq:ev} 
    F_2(x,Q^2) &= \sum_i^{n_f} C_i(x,Q^2) \otimes f_i(x,Q^2) \nonumber \\
    &= \sum_{i,j}^{n_f} C_i(x,Q^2) \otimes \Gamma_{ij}(Q^2,Q_0^2) \otimes f_j(x,Q_0^2),
    \end{align}

where :math:`C_i(x,Q^2)` are the process-dependent coefficient functions which
can be computed perturbatively as an expansion in the QCD and QED
couplings;  :math:`\Gamma_{ij}(Q^2,Q_0^2)` is an evolution operator, determined by the
solutions of the DGLAP equations, which evolves the PDF from the initial
parameterization scale :math:`Q_0^2` into the hard-scattering scale :math:`Q^2`,
:math:`f_i(x,Q^2_0)` are the PDFs at the parameterization scale, and
:math:`\otimes` denotes the Mellin convolution.

The sum over flavors :math:`i,j` runs over the :math:`n_f` active quarks and antiquarks flavors at a given
scale :math:`Q`, as well as over the gluon.

The direct calculation of the above equation during the PDF fit is not practical
since it requires first solving the DGLAP evolution equation for each new boundary
condition at :math:`Q_0` and then convoluting with the coefficient
functions.

To evaluate the observable in a more computationally efficient way, it is better 
to precompute all the perturbative information, i.e. the coefficient functions :math:`C_i`
and the evolution operators :math:`\Gamma_{ij}`, with a suitable
interpolation basis.

Several of these approaches have been made available in the context of
PDF fits.
Here we use the `APFELgrid` to precompute the perturbative
information of the DIS structure functions provided by the `APFEL`.

Within this approach, we can factorize the dependence on the PDFs at the input scale :math:`Q_0` as follows.

First, we introduce an expansion over a set of interpolating functions :math:`\{ I_{\beta}\}` spanning both :math:`Q^2` and :math:`x` such that

.. math::

    \begin{equation}
    f_i(x,Q^2) = \sum_{\beta} \sum_{\tau} f_{i,\beta \tau} I_{\beta}(x) I_{\tau}(Q^2) \, ,
    \end{equation}

where the PDFs are now tabulated
in a grid in the :math:`(x,Q^2)` plane, :math:`f_{i,\beta \tau}\equiv f_i(x_\beta,Q^2_{\tau})`.

We can express this result in terms of the PDFs at the input evolution scale
using the (interpolated) DGLAP evolution operators,

.. math::

    \begin{equation}
    f_{i,\beta \tau} = \sum_j \sum_{\alpha} \Gamma^{\tau}_{ij,\alpha \beta}\,f_j(x_{\alpha},Q_0^2) \, ,
    \end{equation}
so that the nuclear DIS structure function can be
evaluated as

.. math::

    \begin{equation}
    F_2(x,Q^2) = \sum_i^{n_f} C_i(x,Q^2) \otimes \left[
    \sum_{\alpha,\beta,\tau} \sum_j \Gamma^{\tau}_{ij,\alpha \beta}\,f_j(x_{\alpha},Q_0^2) I_{\beta}(x) I_{\tau}(Q^2)\right]\, .
    \end{equation}

This can be rearranged to give

.. math::

    \begin{align}
    \label{eq:ev_interp}
    F_2(x,Q^2) &= \sum_i^{n_f} \sum_{\alpha}^{n_x} FK_{i,\alpha}(x,x_{\alpha},Q^2,Q^2_0) \, f_i(x_{\alpha},Q_0^2) 
    \end{align}

where all of the information about the partonic cross-sections and the DGLAP
evolution operators is now encoded into the so-called FK table, :math:`FK_{i,\alpha}`.

Therefore, with the `APFELgrid` method we are able to
express the series of convolutions by a matrix
multiplication, increasing the numerical 
calculation speed of the DIS structure functions by up to several orders
of magnitude.

**Note: The hadronic processes (proton-proton collision) like Drell-Yan for example follows the same logic, only with slighly more complicated structure of FK tables to take into account the evolution of two PDFs instead of one**
