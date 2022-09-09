.. _thcov_tutorial:

How to include a theory covariance matrix in a fit
==================================================

This section details how to include :ref:`scale variation covariance matrices (covmats) <vptheorycov-index>`
in a PDF fit. At the present time this can only be done at next-to-leading order (NLO), for which the
central theory is :ref:`theory 163 <theory-indexes>`.

First, decide which theory covmat you want
------------------------------------------
- Choose the desired point-prescription listed :ref:`here <prescrips>`.
- Each prescription comes with a ``point_prescription`` flag to include in
  the runcard, one of ["3 point", "5 point", "5bar point", "7 point", "9 point"]

Next, add necessary flags to the runcard
----------------------------------------
- Remember to list the required datasets using ``dataset_inputs`` (see :ref:`data_specification`).
- Add ``theorycovmatconfig`` to the runcard. An example is in the following code snippet:

.. code:: yaml

	############################################################
	theory:
	  theoryid: 163        # database id

	theorycovmatconfig:
	  point_prescription: "3 point"
	  theoryids:
   	    from_: scale_variation_theories
	  pdf: NNPDF31_nlo_as_0118
	  use_thcovmat_in_fitting: true
	  use_thcovmat_in_sampling: true

	############################################################

- ``pdf`` is the PDF used to generate the scale varied predictions which
  construct the theory covmat. Choose something close to the PDF you are
  trying to fit, such as a previous iteration if available.
-  ``theoryids`` are necessary for the construction of the theory covmat.
   To avoid user error in entering them in the correct configuration and order,
   this is handled by the ``produce_scale_variation_theories`` action in
   `config <https://github.com/NNPDF/nnpdf/tree/master/validphys2/src/validphys/config.py>`_,
   using the information in
   `the scalevariations module <https://github.com/NNPDF/nnpdf/tree/master/validphys2/src/validphys/scalevariations>`_.
-  The flags ``use_thcovmat_in_fitting`` and ``use_thcovmat_in_sampling`` specify
   where to use the theory covmat in the code. There are two possible places:
   the fitting (i.e. :math:`\chi^2` minimiser) and the sampling (i.e. pseudodata
   generation). The default is ``True`` for both.

    .. warning::

          Changing either of these to ``False`` will affect the fit outcome and should
          be avoided unless you know what you are doing.


If you want to compare data to another fit
------------------------------------------
-  Sometimes we want to compare data to another fit for validation, for example
   we might want to compare predictions for the NLO fit with MHOUs to the known
   NNLO fit (see :ref:`vptheorycov-tests`).
-  To make sure the cuts match between these two fits, edit the ``datacuts``
   section of the runcard to include the following

.. code:: yaml

	  use_cuts: fromintersection
	  cuts_intersection_spec:
	  - theoryid: 163
	  - theoryid: 53

-  This ensures that the cuts on the data are the intersection of the cuts in
   theory 53 (default NNLO) and theory 163 (central scale variation NLO). See
   :ref:`here <theory-indexes>` for theory definitions.

Example runcard
---------------
The following is an example runcard for an NLO NNPDF3.1-style fit with a 3 point theory covmat.
It can be found `here <https://github.com/NNPDF/nnpdf/tree/master/validphys2/examples/theory_covariance/fit_with_thcovmat.yaml>`_.

.. code:: yaml

	#
	# Configuration file for NNPDF++
	#
	##########################################################################################
	description: Example runcard for NLO NNPDF3.1 style fit with 3pt theory covariance matrix

	##########################################################################################
	# frac: training fraction
	# ewk: apply ewk k-factors
	# sys: systematics treatment (see systypes)
	dataset_inputs:
	  - {dataset: NMCPD, frac: 0.5}
	  - {dataset: NMC, frac: 0.5}
	  - {dataset: SLACP, frac: 0.5}
	  - {dataset: SLACD, frac: 0.5}
	  - {dataset: BCDMSP, frac: 0.5}
	  - {dataset: BCDMSD, frac: 0.5}
	  - {dataset: CHORUSNU, frac: 0.5}
	  - {dataset: CHORUSNB, frac: 0.5}
	  - {dataset: NTVNUDMN, frac: 0.5}
	  - {dataset: NTVNBDMN, frac: 0.5}
	  - {dataset: HERACOMBNCEM, frac: 0.5}
	  - {dataset: HERACOMBNCEP460, frac: 0.5}
	  - {dataset: HERACOMBNCEP575, frac: 0.5}
	  - {dataset: HERACOMBNCEP820, frac: 0.5}
	  - {dataset: HERACOMBNCEP920, frac: 0.5}
	  - {dataset: HERACOMBCCEM, frac: 0.5}
	  - {dataset: HERACOMBCCEP, frac: 0.5}
	  - {dataset: HERAF2CHARM, frac: 0.5}
	  - {dataset: CDFZRAP, frac: 1.0}
	  - {dataset: D0ZRAP, frac: 1.0}
	  - {dataset: D0WEASY, frac: 1.0}
	  - {dataset: D0WMASY, frac: 1.0}
	  - {dataset: ATLASWZRAP36PB, frac: 1.0}
	  - {dataset: ATLASZHIGHMASS49FB, frac: 1.0}
	  - {dataset: ATLASLOMASSDY11EXT, frac: 1.0}
	  - {dataset: ATLASWZRAP11, frac: 0.5}
	  - {dataset: ATLAS1JET11, frac: 0.5}
	  - {dataset: ATLASZPT8TEVMDIST, frac: 0.5}
	  - {dataset: ATLASZPT8TEVYDIST, frac: 0.5}
	  - {dataset: ATLASTTBARTOT, frac: 1.0}
	  - {dataset: ATLASTOPDIFF8TEVTRAPNORM, frac: 1.0}
	  - {dataset: CMSWEASY840PB, frac: 1.0}
	  - {dataset: CMSWMASY47FB, frac: 1.0}
	  - {dataset: CMSDY2D11, frac: 0.5}
	  - {dataset: CMSWMU8TEV, frac: 1.0}
	  - {dataset: CMSZDIFF12, frac: 1.0, cfac: [NRM]}
	  - {dataset: CMSJETS11, frac: 0.5}
	  - {dataset: CMSTTBARTOT, frac: 1.0}
	  - {dataset: CMSTOPDIFF8TEVTTRAPNORM, frac: 1.0}
	  - {dataset: LHCBZ940PB, frac: 1.0}
	  - {dataset: LHCBZEE2FB, frac: 1.0}
	  - {dataset: LHCBWZMU7TEV, frac: 1.0, cfac: [NRM]}
	  - {dataset: LHCBWZMU8TEV, frac: 1.0, cfac: [NRM]}

	############################################################
	datacuts:
	  t0pdfset: 190310-tg-nlo-global                    # PDF set to generate t0 covmat
	  q2min: 13.96                        # Q2 minimum
	  w2min: 12.5                        # W2 minimum
	  combocuts: NNPDF31                 # NNPDF3.0 final kin. cuts
	  jetptcut_tev: 0                    # jet pt cut for tevatron
	  jetptcut_lhc: 0                    # jet pt cut for lhc
	  wptcut_lhc: 30.0                   # Minimum pT for W pT diff distributions
	  jetycut_tev: 1e30                  # jet rap. cut for tevatron
	  jetycut_lhc: 1e30                  # jet rap. cut for lhc
	  dymasscut_min: 0                   # dy inv.mass. min cut
	  dymasscut_max: 1e30                # dy inv.mass. max cut
	  jetcfactcut: 1e30                  # jet cfact. cut
	  use_cuts: fromintersection
	  cuts_intersection_spec:
	  - theoryid: 163
	  - theoryid: 53

	############################################################
	theory:
	  theoryid: 163        # database id

	theorycovmatconfig:
	  point_prescription: "3 point"
	  theoryids:
	   from_: scale_variation_theories
	  fivetheories: None
	  pdf: NNPDF31_nlo_as_0118
	  use_thcovmat_in_fitting: true
	  use_thcovmat_in_sampling: true

	sampling_t0:
	  use_t0: false

	fitting_t0:
	  use_t0: true

	############################################################
	fitting:
	  seed: 65532133530           # set the seed for the random generator
	  genrep: on        # on = generate MC replicas, off = use real data
	  rngalgo: 0        # 0 = ranlux, 1 = cmrg, see randomgenerator.cc
	  fitmethod: NGA    # Minimization algorithm
	  ngen: 30000       # Maximum number of generations
	  nmutants: 80      # Number of mutants for replica
	  paramtype: NN
	  nnodes: [2, 5, 3, 1]

	  # NN23(QED) = sng=0,g=1,v=2,t3=3,ds=4,sp=5,sm=6,(pht=7)
	  # EVOL(QED) = sng=0,g=1,v=2,v3=3,v8=4,t3=5,t8=6,(pht=7)
	  # EVOLS(QED)= sng=0,g=1,v=2,v8=4,t3=4,t8=5,ds=6,(pht=7)
	  # FLVR(QED) = g=0, u=1, ubar=2, d=3, dbar=4, s=5, sbar=6, (pht=7)
	  fitbasis: NN31IC # EVOL (7), EVOLQED (8), etc.
	  basis:
	      # remeber to change the name of PDF accordingly with fitbasis
	      # pos: on for NN squared
	      # mutsize: mutation size
	      # mutprob: mutation probability
	      # smallx, largex: preprocessing ranges
	  - {fl: sng, pos: off, mutsize: [15], mutprob: [0.05], smallx: [1.046, 1.188], largex: [
	      1.437, 2.716]}
	  - {fl: g, pos: off, mutsize: [15], mutprob: [0.05], smallx: [0.9604, 1.23], largex: [
	      0.08459, 6.137]}
	  - {fl: v, pos: off, mutsize: [15], mutprob: [0.05], smallx: [0.5656, 0.7242], largex: [
	      1.153, 2.838]}
	  - {fl: v3, pos: off, mutsize: [15], mutprob: [0.05], smallx: [0.1521, 0.5611], largex: [
	      1.236, 2.976]}
	  - {fl: v8, pos: off, mutsize: [15], mutprob: [0.05], smallx: [0.5264, 0.7246], largex: [
	      0.6919, 3.198]}
	  - {fl: t3, pos: off, mutsize: [15], mutprob: [0.05], smallx: [-0.3687, 1.459], largex: [
	      1.664, 3.373]}
	  - {fl: t8, pos: off, mutsize: [15], mutprob: [0.05], smallx: [0.5357, 1.267], largex: [
	      1.433, 2.866]}
	  - {fl: cp, pos: off, mutsize: [15], mutprob: [0.05], smallx: [-0.09635, 1.204],
	    largex: [1.654, 7.456]}

	############################################################
	stopping:
	  stopmethod: LOOKBACK  # Stopping method
	  lbdelta: 0            # Delta for look-back stopping
	  mingen: 0             # Minimum number of generations
	  window: 500           # Window for moving average
	  minchi2: 3.5          # Minimum chi2
	  minchi2exp: 6.0       # Minimum chi2 for experiments
	  nsmear: 200           # Smear for stopping
	  deltasm: 200          # Delta smear for stopping
	  rv: 2                 # Ratio for validation stopping
	  rt: 0.5               # Ratio for training stopping
	  epsilon: 1e-6         # Gradient epsilon

	############################################################
	positivity:
	  posdatasets:
	  - {dataset: POSF2U, poslambda: 1e6}        # Positivity Lagrange Multiplier
	  - {dataset: POSF2DW, poslambda: 1e6}
	  - {dataset: POSF2S, poslambda: 1e6}
	  - {dataset: POSFLL, poslambda: 1e6}
	  - {dataset: POSDYU, poslambda: 1e10}
	  - {dataset: POSDYD, poslambda: 1e10}
	  - {dataset: POSDYS, poslambda: 1e10}

	############################################################
	closuretest:
	  filterseed: 0     # Random seed to be used in filtering data partitions
	  fakedata: off     # on = to use FAKEPDF to generate pseudo-data
	  fakepdf: MSTW2008nlo68cl      # Theory input for pseudo-data
	  errorsize: 1.0    # uncertainties rescaling
	  fakenoise: off    # on = to add random fluctuations to pseudo-data
	  rancutprob: 1.0   # Fraction of data to be included in the fit
	  rancutmethod: 0   # Method to select rancutprob data fraction
	  rancuttrnval: off # 0(1) to output training(valiation) chi2 in report
	  printpdf4gen: off # To print info on PDFs during minimization

	############################################################
	lhagrid:
	  nx: 150
	  xmin: 1e-9
	  xmed: 0.1
	  xmax: 1.0
	  nq: 50
	  qmax: 1e5

	############################################################
	debug: off
