Chi square figures of merit
================================================================================

Within the NNPDF methodology various figures of merit are used, each of which
can be used in different situations. To avoid confusion, it is important to
understand the differences between the various figures of merit, and to
understand which definition we are referring to in a given context. In
particular, it is worth stressing that whenever a figure of merit is discussed,
the :math:`t_0` method (discussed below) applies.

Here we we provide an overview of the different figures of merit, and discuss
when each of them is used.


The basis of the loss functions: 𝜒²
--------------------------------------------------------------------------------
The :math:`\chi^2` figures of merit used in the NNPDF methodology are all
based on the chi square statistic:

.. math::
    \chi^{2}=\sum_{i, j}^{N_{\text {dat }}}(D-P)_{i} C_{i j}^{-1}(D-P)_{j},

where :math:`D_i` is the :math:`i`-th datapoint, :math:`P_i` is the prediction
of the corresponding datapoint calculated from the convolution product
between the :ref:`FastKernel tables<fktables>` for point :math:`i` and the PDF
model, and :math:`C_{ij}` is the covariance between datapoints :math:`i`
and :math:`j`.

The covariance matrix accounts for correlated systematic uncertainties,
normalization uncertainties, and statistical uncertainties as provided by the
experimental collaborations.

We refer to this figure of merit as *experimental* :math:`\chi^2`.

.. note::
    This definition of :math:`\chi^2` is not used as a figure of merit
    anywhere in NNDPF fits. Instead, variations discussed below
    are used.


Avoiding bias: t₀ method
~~~~~~~~~~~~~~~~~~~~~~~~
The :math:`t_0` method introduced in
:cite:p:`Ball:2009qv` aims to
remove systematic biases as a result of a naive treatment of multiplicative
uncertainties. This is done by redefining the covariance matrix in the
definition of :math:`\chi^2`, resulting in a :math:`t_0` covariance matrix
:math:`C_{t_0}` and a corresponding figure of merit sometimes denoted by
:math:`\chi^2_{t_0}`. The new covariance matrix is constructed by replacing the
central value of the data (which is used as reference for multiplicative
uncertainties) with the theory predictions computed using some existing PDF
set, which needs to be specified.

.. note::
    From NNPDF2.0 onwards the t₀ formalism has been used to define the figure of
    merit used during the fitting of the PDFs.

.. note::

    The :math:`t_0` method is **not** used by default in other :ref:`validphys
    applications <vp-index>`, and instead the default is to compute the
    experimental :math:`\chi^2`. To compute :math:`\chi^2_{t_0}`, users need to
    specify

    .. code-block::  yaml

        use_t0: True
        t0pdfset: <Some LHAPDF set>

    in the relevant :ref:`namespace <namespaces>`. This will instruct actions
    such as :py:func:`validphys.results.dataset_chi2_table` to compute the
    :math:`t_0` estimator.


Missing higher order uncertainties
--------------------------------------------------------------------------------
Another source of uncertainties that we may want to include in the covariance
matrix are theoretical uncertainties, particularly missing higher order
uncertainties estimated through scale variations. These uncertainties can be
considered in the figure of merit through the implementation of a 'theory
covariance matrix'. A paper discussing the formalism can be found here:
:cite:p:`AbdulKhalek:2019bux`. For a tutorial see
:ref:`How to include a theory covariance matrix in a fit <thcov_tutorial>`.


Future test: including PDF errors
--------------------------------------------------------------------------------
To test the generalization power of the NNPDF fitting framework in the region
where PDFs are not constrained by data, the 'future test' has been developed.
The figure of merit considered in a future test is again the :math:`\chi^2`,
however, in this case the covariance matrix is not only the covariance matrix
corresponding to the datasets, but it is instead the sum of the covariance
matrix describing the data uncertainties and the covariance matrix describing
the PDF uncertainties.

For a more detailed discussion of the future test formalism see e.g.
:cite:p:`Cruz-Martinez:2021rgy`, or learn
:ref:`How to run a Future Test <futuretests>`


.. _covmat-reg:

Regularized covariance matrices
--------------------------------------------------------------------------------
Information about the accuracy of the experimental uncertainty is generally not
available, nevertheless inaccuracies in an experimental covariance matrix can
lead to problems during optimization. Simply making a conservative estimate of
the correlations does not always guarantee this problem is avoided and this is
where the regularized covariance matrix comes in: it aims to provide a matrix
which is closely related to the original experimental covariance matrix while
avoiding the problems during optimization.

The stability characteristic for a given dataset can be computed using the
:py:func:`validphys.covmats.covmat_stability_characteristic`. All the dataset
covariance matrices can be altered so that their stability characteristic is
less than a given value by specifying such value as a `norm_threshold`
parameter in the runcard. Adding it in an analysis results in computing a
regularized :math:`\chi^2` that is less sensitive to inaccuracies in the
correlation model. Adding it in a :ref:`fit runcard <runcard-detailed>` results
in a fit with regularized covariance matrices.

.. note::
    There is currently no support for displaying regularized :math:`\chi^2`
    values in :ref:`vp-comparefits <compare-fits>`

A more detailed discussion of regularization procedure, and how it is used
within NNPDF can be found in sections 4.2 and 8.7 of the NNPDF4.0 paper
:cite:p:`nnpdf40`.



The weighted fit method
--------------------------------------------------------------------------------
To determine whether a specific dataset shows inconsistencies with the
global dataset, one can produce a PDF determination in which that measurement
is given an increased weight (usually equal to the combined weight of the other
datasets). The idea being that if -- in oder to accommodate the dataset under
investigation -- the agreement to the other datasets deteriorates, this dataset
is likely inconsistent with the global dataset.

When performing a weighted fit the figure of merit is hence redefined as

.. math::
    \chi^{2}=\frac{1}{N_{\text {dat }}-N_{\text {dat }}^{(j)}}
    \sum_{i \neq j}^{n_{\text {exp }}}N_{\text {dat }}^{(i)}\chi_{i}^{2}
    +\omega^{(j)} \chi_{j}^{2}

with :math:`w^{(j)}=N_{\rm dat}/N^{(j)}_{\rm dat}`.

A dataset can be given an additional weight by explictitly writing a weight key
for a given dataset in the :ref:`n3fit runcard <runcard-detailed>`. For example,
while the default weight is 1, one can set the weight of the
HERACOMB_SIGMARED_C dataset to 100 by adding the following to the runcard:

.. code-block:: yaml

    dataset_inputs:
        - {dataset: HERACOMB_SIGMARED_C, frac: 0.75, weight: 100}


Experimental, validation, and training 𝜒²
--------------------------------------------------------------------------------
When performing a PDF fit we generally distinguish three different definitions
of the :math:`\chi^2` loss function, namely the experimental loss
:math:`\chi^2_{\rm exp}`, the training loss :math:`\chi^2_{rm tr}` and the
validation loss :math:`\chi^2_{val}`, all of which are defined using the
:math:`t_0` method. Here the experimental loss is calculated with respect to the
experimental covariance matrix and corresponding central values, while the
training and validation losses are defined with respect to the central values
of the psuedodata replicas.

The training and validation losses are used for cross-correlation in the
early stopping algorithm, and can further be adjusted to ensure positivity and
integrability of the resulting PDFs after the fit by adding a component to the
loss function (see :ref:`below <lagrange-multipliers>`).

More details of these loss functions and the role they play within the training
of the neural network can be found in the :ref:`methodology overview
<methodology>`.


.. _lagrange-multipliers:
Positivity and integrability: Lagrange multipliers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Generally in an NNPDF fit we will want to ensure positivity and integrability of
the resulting PDFs. This is enforced by means of Lagrange multipliers, which
provide an additional contribution to the definition of the chi squared
loss function.

For an discussion of how exactly the loss function is adjusted upon including
the Lagrange multipliers, see sections 3.1.3 and 3.1.4 of the NNPDF4.0 paper
:cite:p:`nnpdf40`.

An explanation of how the runcard should be adjusted to include the additional
positivity Lagrange multiplier can be found :ref:`here <positivity-label>`,
while the analogous information for integrability can be found 
:ref:`here <integrability-label>`.


Hyperoptimized figure of merit
--------------------------------------------------------------------------------
To test the generalization power of a given methodology (a specific set of
hyperparameter values), we employ hyperoptimization, specifically we use
K-folds cross-validation. The idea of K-folds cross-validation is to create
subsets of data representative of the global dataset, and then perform a
fit to :math:`K-1` subsets while using the :math:`K^{\rm th}` subset as a test
set to check the generalization performance after the neural network has been
trained. The figure of merit that is minimized during the hyperoptimization
routine is obtained by summing over all :math:`K` test losses that are obtained
after performing :math:`K` fits to each possible combination of :math:`K-1`
datasets.

For a more detailed description of the hyperoptimization loss see the
documentation of the :ref:`hyperoptimization algorithm<hyperoptimization>`.
