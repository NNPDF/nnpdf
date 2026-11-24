.. _vptheorycov-pc:

Power corrections
=================

Power corrections (also referred to as higher twist corrections for DIS-like
processes) model contributions from non-perturbative effects that scale as
inverse powers of the hard scale. They are implemented in the
``theorycovariance`` module and can be included as a theory covariance matrix in
a fit. Power corrections for jets and higher twists for DIS data have been
determined in :cite:p:`Ball:2025xtj`, based on NNPDF4.0, where the reader can
find further details on the methodology and phenomenological implications.

The implementation is in
`higher_twist_functions.py <https://github.com/NNPDF/nnpdf/tree/master/validphys2/src/validphys/theorycovariance/higher_twist_functions.py>`_
and
`construction.py <https://github.com/NNPDF/nnpdf/tree/master/validphys2/src/validphys/theorycovariance/construction.py>`_.


Overview
--------

In NNPDF, power corrections modify theoretical predictions by introducing multiplicative
shifts. For a generic observable :math:`O`, the corrected prediction is

.. math:: O \to O \times (1 + \mathrm{PC}),

where :math:`\mathrm{PC}` is the power correction. The shift to the prediction
is therefore

.. math:: \Delta O = O \times \mathrm{PC}.

Different functional forms for the power correction are used depending on the
process type:

- **DIS** (neutral current and charged current): the correction depends on
  Bjorken-:math:`x` and :math:`Q^2`, and scales as :math:`1/Q^2`.
- **Single-inclusive jets**: the correction depends on rapidity and transverse
  momentum :math:`p_T`, and scales as :math:`1/p_T`.
- **Dijets**: the correction depends on a rapidity variable and the dijet
  invariant mass :math:`m_{jj}`, and scales as :math:`1/m_{jj}`.


Parametrisation
---------------

Power corrections are parametrised using a piecewise-linear interpolation
between a set of nodes. The node positions (``nodes``) and the function values
at each node (``yshift``) are specified in the runcard.

The interpolation is constructed as a sum of triangular basis functions: each
node :math:`i` is associated with a triangle that peaks at the node position
with value ``yshift[i]`` and drops linearly to zero at the two neighbouring
nodes. The resulting function is continuous and piecewise-linear.

For DIS processes, the nodes are placed in Bjorken-:math:`x` and the power
correction for a data point at :math:`(x, Q^2)` is

.. math:: \mathrm{PC}(x, Q^2) = \frac{h(x)}{Q^2},

where :math:`h(x)` is the piecewise-linear interpolation.

For jet processes, the nodes are placed in rapidity and the correction at
:math:`(y, p_T)` is

.. math:: \mathrm{PC}(y, p_T) = \frac{h(y)}{p_T}.

For dijets, the same functional form is used but the suppression scale is the
dijet invariant mass :math:`m_{jj}`.


Dataset routing
---------------

Each dataset is mapped to one or more power correction parameter keys via the
function ``get_pc_type``. The mapping depends on the process type and dataset
name:

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Process type
     - PC type key
     - Datasets
   * - DIS NC (proton :math:`F_2`)
     - ``f2p``
     - SLAC, BCDMS proton :math:`F_2`; NMC, HERA :math:`\sigma_{\mathrm{red}}`
   * - DIS NC (deuteron :math:`F_2`)
     - ``f2d``
     - SLAC, BCDMS deuteron :math:`F_2`
   * - DIS NC (NMC ratio :math:`F_2^d / F_2^p`)
     - ``(f2p, f2d)``
     - NMC ratio dataset
   * - DIS CC
     - ``dis_cc``
     - CHORUS, NuTeV, HERA CC
   * - Jets
     - ``Hj``
     - Single-inclusive jet datasets
   * - Dijets (ATLAS)
     - ``H2j_ystar``
     - ATLAS dijet datasets
   * - Dijets (CMS)
     - ``H2j_ymax``
     - CMS dijet datasets
   * - Factorised dijets
     - ``(H2j_ystar, H2j_yb)``
     - Triple-differential dijets and their one-dimensional projections


Special case: NMC ratio
~~~~~~~~~~~~~~~~~~~~~~~~

The NMC ratio dataset :math:`F_2^d / F_2^p` receives contributions from both
the proton and deuteron power corrections. The corrected ratio is

.. math::

   \frac{F_2^d}{F_2^p} \to \frac{F_2^d \,(1 + \mathrm{PC}_d)}{F_2^p \,(1 + \mathrm{PC}_p)},

and the shift is

.. math::

   \Delta\!\left(\frac{F_2^d}{F_2^p}\right) = \frac{F_2^d}{F_2^p} \,
   \frac{\mathrm{PC}_d - \mathrm{PC}_p}{1 + \mathrm{PC}_p}.


Covariance matrix construction
------------------------------

The theory covariance matrix is constructed from the shifts :math:`\Delta O` by
taking outer products. For each combination of power correction parameters, a
shift vector is computed per dataset. The sub-matrix between datasets :math:`i`
and :math:`j` is then

.. math::

   S_{ij} = \sum_k \Delta_i^{(k)} \otimes \Delta_j^{(k)},

where :math:`k` runs over all parameter combinations (one non-zero ``yshift``
entry at a time, with all others set to zero). This corresponds to the
``covmat_power_corrections`` function in ``construction.py``.


Runcard configuration
---------------------

Power corrections are included via the ``theorycovmatconfig`` section of the
runcard. The key ``"power corrections"`` must be added to the
``point_prescriptions`` list, alongside any scale variation prescriptions.

The following keys are used:

- ``pc_parameters``: a dictionary mapping PC type keys to their parametrisation
  (``yshift`` and ``nodes`` arrays). The length of ``yshift`` must match the
  length of ``nodes``.
- ``pc_included_procs``: list of process types to which power corrections
  are applied (e.g. ``["DIS NC", "DIS CC", "JETS", "DIJET"]``).
- ``pc_excluded_datasets``: list of dataset names to exclude from power corrections
  even if their process type is included.
- ``pdf``: the PDF used for computing the theory predictions that enter the
  multiplicative shifts.

Example
~~~~~~~

.. code:: yaml

    theorycovmatconfig:
      point_prescriptions: ["9 point", "power corrections"]
      pc_parameters:
        f2p:
          yshift: [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0]
          nodes:  [0.0, 0.001, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
        f2d:
          yshift: [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0]
          nodes:  [0.0, 0.001, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
        dis_cc:
          yshift: [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0]
          nodes:  [0.0, 0.001, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
        Hj:
          yshift: [2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
          nodes:  [0.25, 0.75, 1.25, 1.75, 2.25, 2.75]
        H2j_ystar:
          yshift: [2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
          nodes:  [0.25, 0.75, 1.25, 1.75, 2.25, 2.75]
        H2j_ymax:
          yshift: [2.0, 2.0, 2.0, 2.0, 2.0]
          nodes:  [0.25, 0.75, 1.25, 1.75, 2.25]
      pc_included_procs: ["JETS", "DIJET", "DIS NC", "DIS CC"]
      pc_excluded_datasets:
        - HERA_NC_318GEV_EAVG_CHARM-SIGMARED
        - HERA_NC_318GEV_EAVG_BOTTOM-SIGMARED
      pdf: NNPDF40_nnlo_as_01180
      use_thcovmat_in_fitting: true
      use_thcovmat_in_sampling: true

.. warning::
   The lengths of ``yshift`` and ``nodes`` must be equal for each PC type.
   A mismatch will raise an error at initialisation time.

.. note::
   Power corrections can be combined with scale variation prescriptions.
   Both contributions are summed into a single theory covariance matrix.
   See the tutorial on :ref:`including a theory covmat in a fit <thcov_tutorial>`.


.. _pc-posterior:

Posterior analysis after fitting
--------------------------------

The posterior extraction is available directly in ``validphys``. It requires a
completed, postprocessed Monte Carlo PDF fit with ``"power corrections"`` in
``theorycovmatconfig.point_prescriptions`` and
``use_thcovmat_in_fitting: true``. Run the fit with ``savepseudodata: true`` and
retain its saved pseudodata and ``tables`` directory. The prior PDF specified by
``theorycovmatconfig.pdf`` must also be available.

Create a runcard such as
``validphys2/examples/higher_twist_posterior.yaml``, replacing the fit name:

.. code:: yaml

    fit: MY_PC_PRIOR_FIT
    pc_cfactor_label: PC

    actions_:
      - pc_fit_context pc_posterior_table
      - pc_fit_context pc_posterior_covmat
      - pc_fit_context pc_posterior_corrmat
      - pc_fit_context plot_pc_posterior
      - pc_fit_context plot_pc_posterior_corrmat
      - pc_fit_context pc_posterior_cfactors

Then run:

.. code:: bash

    validphys posterior.yaml -o posterior

``pc_fit_context`` supplies the fit's dataset inputs, theory, cuts, fitted PDF,
prior PDF and PC settings. There is no need to copy these settings into the
analysis runcard. The fitted PDF supplies the replica predictions; the prior
PDF supplies the shifts used to construct the original theory covariance.
The process selection uses ``nnpdf31_process`` names, as in the prior
construction. The older ``pc_excluded_exps`` spelling is also read when
``pc_excluded_datasets`` is absent.

The calculation
~~~~~~~~~~~~~~~

Let :math:`B_{ia}` be the shift of data point :math:`i` under prior combination
:math:`a`, computed by ``compute_deltas_pc``, and let :math:`A_{\alpha a}` be the
corresponding shift of node parameter :math:`\alpha`. Define

.. math::

    S = BB^T, \qquad \widetilde S = AA^T, \qquad \widehat S = AB^T,
    \qquad K = C + S.

Here :math:`C` contains the experimental covariance and all non-PC theory
contributions saved with the fit. The experimental part uses data-normalised
multiplicative uncertainties
(``use_t0: false``). It is built by experiment before reordering by process,
preserving experimental correlations between different process groups.
The PC shifts are checked against the saved PC covariance
before conditioning; a mismatch indicates that the original prior cannot be
reconstructed with the current inputs.

With :math:`\overline T` the mean prediction over the postfit PDF replicas,
:math:`\overline D` their mean saved pseudodata, and
:math:`X=N_{\rm rep}^{-1}\sum_r(T_r-\overline T)(T_r-\overline T)^T`, the
zero-centred prior gives

.. math::

    \overline h &= \widehat S K^{-1}(\overline D-\overline T), \\
    P &= \widetilde S-\widehat S K^{-1}\widehat S^T
          +\widehat S K^{-1}XK^{-1}\widehat S^T.

The implementation uses linear solves and projects replica deviations into
parameter space, avoiding an explicit inverse and the full matrix :math:`X`.
For diagonal-basis fits, saved pseudodata are rotated back to the data basis
and aligned with predictions by dataset and point ID. Fixed nodes with zero
prior width are retained.

For factorised dijets (both ``H2j_ystar`` and ``H2j_yb`` present), each product
combination shifts both participating factor nodes. Consequently
:math:`AA^T` is not a diagonal prior with entries ``yshift**2`` and can be
singular. The tables describe these factor-node parameters. The plots instead
show the GeV-valued marginal corrections obtained by averaging over the other
rapidity coordinate. Their means and variances use the Gaussian product
moments, including :math:`\langle ab\rangle=\langle a\rangle\langle b\rangle+
\operatorname{Cov}(a,b)`.

Outputs and c-factors
~~~~~~~~~~~~~~~~~~~~~

The actions write standard validphys tables under ``posterior/tables`` and
figures under ``posterior/figures``:

* ``pc_posterior_table``: node position, posterior mean, posterior standard
  deviation and prior standard deviation, indexed by PC type and node label.
* ``pc_posterior_covmat``: the full parameter covariance, including correlations
  between different PC types. Use this matrix when propagating uncertainties.
* ``pc_posterior_corrmat``: the corresponding correlations. Fixed nodes have
  zero rows and columns.
* ``plot_pc_posterior``: prior and posterior one-standard-deviation bands for
  each PC type. The bands propagate correlated node uncertainties through the
  piecewise-linear interpolation, without Monte Carlo sampling.
* ``plot_pc_posterior_corrmat``: the parameter correlation heatmap.
* ``pc_posterior_cfactors``: a manifest of the c-factor files, written under
  ``posterior/cfactors/MY_PC_PRIOR_FIT``.
* ``pc_prediction_mean_table``: the mean correction in prediction units.
  Factorised products include their node cross-covariance; ratio corrections
  are evaluated at the node means.
* ``pc_shift_covmat``: the per-fit second-moment prescription
  :math:`f C_{\rm post}+\Delta T\Delta T^T`, with
  ``pc_covariance_factor: 2.0`` by default. The factor scales only the
  posterior covariance, not the mean outer product.

The c-factor files are named ``CF_PC_<FK-table-name>.dat`` (replace ``PC`` using
``pc_cfactor_label``). They cover the full FK-table binning, including bins
removed by the fit cuts. Split tables and metadata bin shifts are handled
per FK operand. For the NMC deuteron/proton ratio, separate factors are written
for numerator and denominator. Factorised-dijet central factors include the
posterior cross-covariance term.

Export requires PineAPPL FK tables. Copy the files to the appropriate theory's
``cfactor`` directory and add ``PC`` to the dataset's ``cfac`` list to apply
the central correction in a subsequent calculation. The second column is zero:
these files apply central shifts and do not encode the
correlated posterior uncertainty. That uncertainty is supplied by
``pc_posterior_covmat``. For a ratio observable, applying the two central
operand factors gives their ratio; this is not the expectation of a ratio of
Gaussian random variables.

These actions can also be used inside reports or run over ``fits`` by prefixing
each action with ``fits pc_fit_context``. The numerical provider is
``pc_posterior`` in
``validphys.theorycovariance.higher_twist_posterior``; it returns the mean,
posterior covariance and prior covariance as labelled pandas objects.


Report and prediction-space uncertainties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a complete HTML report, set the fit name in
``validphys2/examples/higher_twist_posterior_report.yaml`` and run:

.. code:: bash

    validphys validphys2/examples/higher_twist_posterior_report.yaml -o pc-report

Open ``pc-report/index.html``. The report includes the posterior function bands,
the node table and parameter correlation heatmap, followed by:

* ``plot_pc_prediction_uncertainties``: prior and posterior PC standard deviations
  at fitted data points, shown as percentages of the absolute mean PDF-replica
  prediction. Each process has a separate plot, with dataset boundaries marked.
* ``plot_pc_prediction_experimental_uncertainties``: posterior PC versus
  experimental standard deviations, both shown as percentages of the absolute
  measured central value. This compares the uncertainty on the correction,
  rather than the size of its central value.
  This uses :math:`\sqrt{(C_{\rm post})_{ii}}` and
  :math:`\sqrt{(C_{\rm exp})_{ii}}`; the posterior is uninflated and excludes
  the mean-shift outer product. Plots are grouped by process with dataset
  boundaries, using the same fitted points as the prior/posterior comparison.
* ``plot_pc_prediction_corrmat``: correlations of the PC shifts between fitted
  points, including correlations between different datasets and processes.
  Points with zero prior PC uncertainty are omitted from these plots.

The per-point values are saved by ``pc_prediction_uncertainties``. Relative
uncertainties are undefined (``NaN``) where the mean prediction vanishes.
These are uncertainties of the PC contribution; they do not represent the total
uncertainty of the PDF prediction plus its correction.

``pc_prediction_experimental_uncertainties`` exports the absolute errors,
data-normalised relative errors and the ratio
:math:`\sqrt{(C_{\rm post})_{ii}/(C_{\rm exp})_{ii}}`. The experimental diagonal
is aligned by dataset and original point ID, using the same data-normalised
experimental covariance as the posterior extraction (without theory or PDF
errors added to it). Relative errors are ``NaN`` at zero measured values;
the ratio is ``NaN`` at zero experimental error. Zero-PC points remain in the
table, while only points with a nonzero PC prior are plotted. This comparison
reuses the existing posterior diagonal calculation and does not construct
another full prediction-space covariance matrix.

The prediction-space covariance :math:`C_{\rm post}` is a projection of the
posterior combination covariance. In terms of the posterior
covariance :math:`V` of the independent, unit-variance prior combinations,

.. math::

    V &= I-B^TK^{-1}B+B^TK^{-1}XK^{-1}B, \\
    P &= AVA^T, \\
    C_{\rm post} &= BVB^T
       = S-SK^{-1}S+SK^{-1}XK^{-1}S, \qquad S=BB^T.

Both projections share ``pc_combination_posterior``, so conditioning is performed
once. The prediction-space calculation keeps the full correlations and works
for singular parameter bases, including factorised dijets. The pointwise
uncertainty plots evaluate only the diagonal and do not require a full
data-space matrix. Here :math:`S` is reconstructed from the prior shifts;
comparison with an analysis using the saved, float32-rounded prior covariance
can therefore show small rounding differences.

``pc_prediction_covmat`` supplies :math:`C_{\rm post}` itself. The second-moment
prescription :math:`2C_{\rm post}+\Delta T\Delta T^T` additionally needs the mean
shift and the chosen covariance inflation factor; neither is folded into this
provider's output.

Here the shifts are those of the fitted prior-combination model. For nonlinear
ratio or factorised-product corrections, this covariance follows the
shift-model prescription; it is not the exact nonlinear propagation of a
Gaussian distribution over the reported node parameters through the c-factor
functions.

The full matrices cover all fitted points, with zero rows and columns for
unaffected points. To export them, add either or both actions to the report
runcard (the global-fit tables can be large):

.. code:: yaml

    actions_:
      - pc_fit_context pc_posterior_covmat
      - pc_fit_context pc_prediction_covmat_table
      - pc_fit_context pc_prediction_corrmat_table
      - report(main=true)

The underlying providers ``pc_prediction_covmat`` and ``pc_prediction_corrmat``
return labelled dataframes for use in other validphys analyses. Central
c-factors can also be exported alongside the report by adding
``pc_fit_context pc_posterior_cfactors``.


Jet correction size and kinematic cuts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To identify where the fitted correction is larger than the experimental error,
use the criterion

.. math::

    |C_{{\rm PC},i}-1| > \frac{\sqrt{(C_{\rm exp})_{ii}}}{|D_i|}.

The left-hand side is the **central correction**, not its posterior standard
deviation. A large, well-determined correction can satisfy this criterion even
when its posterior uncertainty is much smaller than the experimental error.
For example, ``260710-ac-01-nnpdf41-pc`` gives a central correction of 21.86%
for ``ATLAS_1JET_8TEV_R06_PTY`` at :math:`p_T=77.5` GeV and
:math:`|y|=0.25`, compared with an experimental error of 9.53%.
The data-normalised posterior width at that point is only 2.72%.

The posterior report includes ``plot_pc_jet_correction_comparison`` and downloads
of ``pc_jet_correction_comparison`` and ``pc_jet_correction_exceedance``.
For a report over several fits, edit ``fits`` in
``validphys2/examples/higher_twist_jet_correction_report.yaml`` and run:

.. code:: bash

    validphys validphys2/examples/higher_twist_jet_correction_report.yaml -o pc-jet-report

The example selects the five original ``260710-ac-{01,02,03,04,05}-nnpdf41-pc``
fits. Each selected fit supplies its own posterior, prior settings, dataset
variants and fitted cuts. Installed PC c-factor files are not read. Fits are
resolved sequentially, retaining only compact comparison tables between fits.

The plots show signed :math:`C_{\rm PC}-1`, its posterior band, and the positive
experimental relative error, with separate rapidity bins and :math:`p_T` or
:math:`m_{jj}` ticks in GeV. Triple-differential dijets have one row per
:math:`y_b` bin, separated into :math:`y^*` bins. The criterion always uses the
absolute correction, including when the plotted central correction is negative.
Factorised-dijet means include the covariance contribution to the mean product.

The bands propagate the posterior combination covariance directly into
dimensionless c-factor units. If :math:`R_{ia}` is the prior c-factor response
minus one for combination :math:`a`, the displayed width is
:math:`\sqrt{(RVR^T)_{ii}}`. The prediction-space uncertainty comparison instead
divides the width by the measured value; for a multiplicative jet correction this
introduces an additional factor :math:`|T_{{\rm prior},i}|/|D_i|`. Neither covariance inflation
nor the mean outer product enters these bands. For factorised dijets they use
the fitted shift-combination prescription, not exact nonlinear product propagation.

The per-point table preserves dataset names and original point IDs, with
``correction``, ``correction_std`` and ``experimental_relative_std`` expressed
as fractions. ``exceeds_experiment`` is undefined at zero measured values.
The summary reports ``ndata``, ``n_compared``, ``n_exceeding`` and
``max_scale_exceeding`` for every rapidity bin, along with its fitted scale
range. An empty maximum means no compared point satisfies the criterion.
The multi-fit providers prefix these tables with a ``fit`` index:
``pc_fits_jet_correction_comparison`` and ``pc_fits_jet_correction_exceedance``.

``max_scale_exceeding`` is the highest observed exceeding bin-centre scale, not an
interpolated crossing or an automatic cut recommendation. Exceedance need not
be monotonic, and points excluded by the fit are not tested. Use the pointwise
table and rapidity-resolved plots to examine candidate cuts within this coverage.


.. _pc-posterior-merging:

Merging partitioned fits
~~~~~~~~~~~~~~~~~~~~~~~~

``validphys2/examples/higher_twist_merge.yaml`` merges several PC fits into one
labelled data-space matrix. Each dataset has one owner. By default, the first
fit is the reference and owns all shared datasets; a dataset absent from the
reference belongs to the first fit in the list containing it. Thus shared DIS
is counted once, while each distinct jet sector comes from its own fit.

What the two prescriptions do
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Both methods treat the correction as an uncertainty and return
:math:`C_{\rm merged}=fC_{\rm post}+\Delta T\Delta T^T`, where
:math:`\Delta T` is the merged mean correction and :math:`f` scales its
posterior covariance. The mean outer product is kept across **all** sectors.

**Independent merge** (``pc_merge_method: independent``).
Keep each owner's mean correction and posterior covariance block unchanged,
and set posterior covariances between different fits to zero:

.. math::

    C_{\rm post}=\operatorname{blockdiag}
    (C_{{\rm post},1},\ldots,C_{{\rm post},K}).

Correlations within an owner's block are preserved. Different owners are
correlated only through :math:`\Delta T\Delta T^T`; the DIS--jet posterior
correlations measured in the non-reference fits are discarded.

**DIS-conditional merge** (``pc_merge_method: conditional``).
Use one common DIS posterior from the reference fit. For every jet sector,
keep the relation learned in its own fit between its parameters and DIS,
together with the residual uncertainty left when DIS is fixed. Conceptually,
draw DIS once from the reference posterior, then draw each sector from its
own conditional distribution given that same DIS draw. Sector residuals
are assumed independent given DIS.

This induces DIS--jet correlations and jet--jet correlations mediated by
DIS. The reference fit's block is preserved, while the other sectors' means
and marginal uncertainties generally change to follow the reference DIS
determination. It assumes compatible DIS determinations and transferable
conditional relations; it does not reconstruct correlations from a joint
fit of all sectors.

In short, the independent prescription preserves the separate determinations;
the conditional prescription makes them share one fluctuating DIS
determination. Both include the coherent mean-shift outer product and give
positive-semidefinite matrices for valid input covariances and :math:`f\geq0`.

Running and configuring a merge
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: yaml

    fits:
      - FIT_DIS_AND_LARGE_R_JETS
      - FIT_DIS_AND_SMALL_R_JETS
      - FIT_DIS_AND_DIJETS

    pc_merge_method: independent  # or conditional
    pc_covariance_factor: 2.0
    pc_merge_reference: FIT_DIS_AND_LARGE_R_JETS  # optional; defaults to first fit

    actions_:
      - pc_merged_covmat
      - pc_merged_mean_shift
      - pc_merge_dis_comparison
      - pc_merge_dis_covmats

Run with ``validphys higher_twist_merge.yaml -o pc-merge``. The final table is
``pc-merge/tables/pc_merged_covmat.csv``. All matrices use
``(group, dataset, id)`` labels. The companion ``pc_merged_mean_shift`` table
records both the correction and the owner of every point.

To display the merged matrix in an HTML report, use
``validphys2/examples/higher_twist_merge_report.yaml``:

.. code:: bash

    validphys validphys2/examples/higher_twist_merge_report.yaml -o pc-merge-report

Open ``pc-merge-report/index.html``. The runcard accepts the same method,
factor, reference, ownership and subset settings as the table-only example.
The report includes:

* ``plot_pc_merged_covmat``: the final
  :math:`f C_{\rm post}+\Delta T\Delta T^T`, in products of prediction units,
  with a symmetric logarithmic colour scale.
* ``plot_pc_merged_corrmat``: the correlation matrix of that final error model,
  including the mean outer product. This permits comparison across observables
  with different units; it is not the correlation matrix of :math:`C_{\rm post}`
  alone.
* The input DIS comparison table and downloads of the final covariance,
  mean shifts with point owners, and input DIS covariance matrices.

Both heatmaps use dataset names as axis labels and mark dataset boundaries,
preserving the point order of the exported matrix.
Zero-variance points are omitted from the plots unless the entire matrix
vanishes; the exported matrices always retain every selected point.

``pc_covariance_factor`` may be any finite non-negative number. Both methods
apply it only after forming the merged posterior covariance. The optional
``pc_merged_prediction_covmat`` action exports that uninflated covariance.

To override an owner, supply a mapping such as
``pc_dataset_owners: {DATASET_NAME: FIT_NAME}``. To select only part of the
union, supply ``pc_merge_datasets: [DATASET_NAME, ANOTHER_DATASET]``. Output
order is reference first, then the remaining fits in list order, preserving
each owner's point order. The subset list does not reorder the output.
Only the owner's fitted points are exported: cuts from different fits are
not united, and excluded points are not extrapolated. Zero-PC points are
retained with zero shifts and covariance rows. When using the table in a
subsequent fit, its selected PC-affected points must be covered by these
owner cuts. Distinct dataset names referring to overlapping measurements
must be resolved by the chosen dataset subset; ownership cannot detect
experimental overlap from names alone.

Compatibility and nonlinear corrections
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For ``conditional``, all selected DIS datasets must belong to the reference.
The DIS parameter types, node coordinates and exactly fixed nodes must agree
across contributing fits. Fixed nodes are removed before the Cholesky solve;
the DIS covariance must be positive definite on the remaining nodes.
``pc_merge_dis_comparison`` and ``pc_merge_dis_covmats`` expose the original
DIS means, standard deviations and full covariances for comparison. The
mean differences normalised by combined marginal errors are descriptive,
not independent-fit significance tests: the fits share DIS data. There is
no automatic statistical compatibility threshold.

The conditional calculation operates in the original shift-combination
space and anchors it to the physical DIS nodes. This is equivalent to the
node-parameter prescription for linear corrections, and preserves the
fitted covariance for factorised dijets without inverting their
singular node-to-combination map. The correction means use the same evaluation
as ``pc_prediction_mean`` after transporting the
node posterior. For products this includes their cross-covariance; for
ratios it is a plug-in evaluation. The covariance remains the fitted
shift-model covariance, not exact nonlinear propagation through these
corrections.

Fits are resolved sequentially through a single provider, keeping compact
combination posteriors and response matrices. This avoids retaining five
full experimental-covariance dependency graphs in memory at once. The
per-fit resolution happens during execution, so ``--dry`` checks the outer
merge graph rather than every inner fit computation.

Module reference
----------------

``higher_twist_functions.py`` provides the following public functions:

- ``get_pc_type(exp_name, process_type, experiment, pc_dict)``:
  determines which PC type key(s) apply to a given dataset.
- ``linear_bin_function(a, y_shift, bin_edges)``:
  evaluates the piecewise-linear triangular interpolation at points ``a``.
- ``dis_pc_func(delta_h, nodes, x, Q2)``:
  computes the DIS power correction :math:`h(x)/Q^2`.
- ``jets_pc_func(delta_h, nodes, pT, rap)``:
  computes the jet power correction :math:`h(y)/p_T`.
- ``mult_dis_pc(nodes, x, q2, dataset_sp, pdf)``:
  returns a function that computes the multiplicative DIS shift given node values.
- ``mult_dis_ratio_pc(p_nodes, d_nodes, x, q2, dataset_sp, pdf)``:
  returns a function that computes the shift for the :math:`F_2^d/F_2^p` ratio.
- ``mult_jet_pc(nodes, pT, rap, dataset_sp, pdf)``:
  returns a function that computes the multiplicative jet shift given node values.
- ``construct_pars_combs(parameters_dict)``:
  builds the list of one-at-a-time parameter combinations used to construct
  the covariance matrix.
- ``compute_deltas_pc(dataset_sp, pdf, pc_dict)``:
  computes the full set of shifts for a single dataset.

``construction.py`` provides:

- ``covmat_power_corrections(deltas1, deltas2)``:
  computes the theory covariance sub-matrix between two datasets from their
  shift dictionaries.
- ``covs_pt_prescrip_pc(combine_by_type, point_prescription, pdf, pc_parameters, pc_included_procs, pc_excluded_datasets)``:
  assembles the full power correction covariance matrix across all datasets.
