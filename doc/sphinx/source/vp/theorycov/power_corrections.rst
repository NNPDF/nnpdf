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
     - ``H2j_ATLAS``
     - ATLAS dijet datasets (falls back to ``H2j`` if key absent)
   * - Dijets (CMS)
     - ``H2j_CMS``
     - CMS dijet datasets (falls back to ``H2j`` if key absent)


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
- ``pc_excluded_exps``: list of dataset names to exclude from power corrections
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
        H2j_ATLAS:
          yshift: [2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
          nodes:  [0.25, 0.75, 1.25, 1.75, 2.25, 2.75]
        H2j_CMS:
          yshift: [2.0, 2.0, 2.0, 2.0, 2.0]
          nodes:  [0.25, 0.75, 1.25, 1.75, 2.25]
      pc_included_procs: ["JETS", "DIJET", "DIS NC", "DIS CC"]
      pc_excluded_exps:
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
- ``covs_pt_prescrip_pc(combine_by_type, point_prescription, pdf, pc_parameters, pc_included_procs, pc_excluded_exps)``:
  assembles the full power correction covariance matrix across all datasets.
