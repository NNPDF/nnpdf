.. _theory-covmat-examples:

Examples
========

Covariance matrix report
------------------------

The following is an example of how to produce a report detailing the covariance
matrices. In this case, the 3 point prescription is shown, for a global data
at NLO.

You need to provide the central theory under the ``default_theory`` flag, 
corresponding to :math:`(\mu_F, \mu_R) = (0,0)`, which for NLO is theory 163.

You need to provide the required point prescription using the flag in 
:ref:`this section <pointprescrips>`, e.g. ``point_prescription: "3 point"``
in the case below.

``dataspecs`` associates a chosen label (``speclabel``) with each of the theory
choices. This details what scale variation the theory corresponds to.

Here the cuts and PDF are taken from the central NLO scale-varied fit.

You must also list all the datasets you wish to include, along with any 
relevant c-factors. 

.. warning::
	In order to ensure backwards compatibility now that the structure
	of data in runcards has been updated and ``experiments`` is deprecated, you must
	also include ``metadata_group: nnpdf31_process`` in the runcards, so that the
	scale variation prescriptions are done by process rather than by experiment. See
	:ref:`backwards-compatibility` for more details.

.. code-block::  yaml
   :linenos:

   meta:
     author: Andrea Barontini
     keywords: [theory uncertainties, 3-point]
     title: NLO 3-point variations for 5 process types - DIS CC, DIS NC, DY, Top, Jets

   default_theory:
     - theoryid: 163

   theoryid: 163
     point_prescription: '3 point'


   theoryids:
     from_: scale_variation_theories

   dataspecs:
     - theoryid: 163
      speclabel: $(\xi_F,\xi_R)=(1,1)$
     - theoryid: 180
      speclabel: $(\xi_F,\xi_R)=(2,2)$
     - theoryid: 173
      speclabel: $(\xi_F,\xi_R)=(0.5,0.5)$

   normalize_to: 1

   fit: 190315_ern_nlo_central_163_global
   use_cuts: "fromfit"

   pdf:
     from_: fit

   dataset_inputs:
     - dataset: NMCPD
     - dataset: NMC
     - dataset: SLACP
     - dataset: SLACD
     - dataset: BCDMSP
     - dataset: BCDMSD
     - dataset: CHORUSNU
     - dataset: CHORUSNB
     - dataset: NTVNUDMN
     - dataset: NTVNBDMN
     - dataset: HERACOMBNCEM
     - dataset: HERACOMBNCEP460
     - dataset: HERACOMBNCEP575
     - dataset: HERACOMBNCEP820
     - dataset: HERACOMBNCEP920
     - dataset: HERACOMBCCEM
     - dataset: HERACOMBCCEP
     - dataset: HERAF2CHARM
     - dataset: CDFZRAP
     - dataset: D0ZRAP
     - dataset: D0WEASY
     - dataset: D0WMASY
     - dataset: ATLASWZRAP36PB
     - dataset: ATLASZHIGHMASS49FB
     - dataset: ATLASLOMASSDY11EXT
     - dataset: ATLASWZRAP11
     - dataset: ATLAS1JET11
     - dataset: ATLASZPT8TEVMDIST
     - dataset: ATLASZPT8TEVYDIST
     - dataset: ATLASTTBARTOT
     - dataset: ATLASTOPDIFF8TEVTRAPNORM
     - dataset: CMSWEASY840PB
     - dataset: CMSWMASY47FB
     - dataset: CMSDY2D11
     - dataset: CMSWMU8TEV
     - {dataset: CMSZDIFF12, cfac: [NRM]}
     - dataset: CMSJETS11
     - dataset: CMSTTBARTOT
     - dataset: CMSTOPDIFF8TEVTTRAPNORM
     - dataset: LHCBZ940PB
     - dataset: LHCBZEE2FB
     - {dataset: LHCBWZMU7TEV, cfac: [NRM]}
     - {dataset: LHCBWZMU8TEV, cfac: [NRM]}


   template: template_matrix_plots.md

   dataset_report:
     meta: Null
     template_text: |
        ## Scale variations as a function of the kinematics for {@dataset_name@}
        {@plot_fancy_dataspecs@}

   actions_:
     - report(main=true)


 


The corresponding template file is ``template_matrix_plots.md``, shown below. This will produce
a comprehensive set of plots and tables describing the covariance matrices.

.. code-block::  md
   :linenos:

   Covariance matrices
   -------------------
   {@with default_theory@}
      {@plot_normexpcovmat_heatmap@}
      {@plot_normthcovmat_heatmap_custom@}
   {@endwith@}

   Correlation matrices
   --------------------
   {@with default_theory@}
      {@plot_expcorrmat_heatmap@}
      {@plot_thcorrmat_heatmap_custom@}
      {@plot_expplusthcorrmat_heatmap_custom@}
   {@endwith@}

   Diagonal elements of covariance matrices
   ----------------------------------------
   {@with default_theory@}
      {@plot_diag_cov_comparison@}
   {@endwith@}

Validation report
----------------- 

Here is an example of a runcard for a report validating the theory covariance
matrix against the NNLO-NLO shift. In this case the 5bar point prescription is chosen,
and only a single Drell-Yan dataset is considered.

The flag ``orthonormalisation`` corresponds to the method used to orthonormalise 
the basis vectors of the theory covariance matrix. There are three choices:

#. QR decomposition (choose this by default), with the flag ``qr``

#. Singular value decomposition, with the flag ``svd``

#. An in-built Gram-Schmidt orthonormalisation, with the flag ``gs``.

``_datasets_list_nlo`` is a list of all the experiments to be included at NLO.
Defining them as a list here avoids the need to repeat the same block of text
many times later on for each theory.

The remainder of the runcard is divided into two namespaces, ``shiftconfig`` and
``theoryconfig``. The former deals with the information concerning the NNLO-NLO
shift vector, and the latter with the information needed to construct the theory
covariance matrix.

In ``shiftconfig`` we provide an NLO and an NNLO dataspec, so that the shift can
be calculated as the difference between the two. Here we list just the experiments
we wish to consider, e.g. Drell-Yan experiments in this case. Because the experiments
and cuts are matched between ``theoryconfig`` and ``shiftconfig`` this means that
overall only these experiments will be used, even though we can pass the whole
``_experiments_list_nlo`` list to ``theoryconfig``.

In ``theoryconfig`` we again provide the relevant theories, in the correct order.
For each dataspec we can give the ``_experiments_list_nlo``. 

.. code-block::  yaml
   :linenos:

   meta:
     title: Theory shift validation test, 5bar point, single-dataset, QR
     author: Andrea Barontini
     keywords: [example]

   # This is a single dataset comparison. For a full comparison extend
   # the datasets accordingly.

   metadata_group: nnpdf31_process #!IMPORTANT! Do not remove

   orthonormalisation: qr # Choice of orthonormalisation scheme
                       # for finding th covmat basis. Default is qr.
   theoryid: 163

   fit: 190315_ern_nlo_central_163_global

   use_cuts: fromfit

   pdf:
     from_: fit

   _datasets_list_nlo: &datasets_list_nlo
     - dataset: D0WMASY

   shiftconfig:     # For calculating NNLO-NLO shift

   use_cuts: fromfit
   fit: 190315_ern_nlo_central_163_global
  
   theoryid: 163

   dataspecs:
       - theoryid: 163
         pdf:
           from_: fit
         speclabel: "NLO"
         dataset_inputs:
                  - dataset: D0WMASY
       - theoryid: 166
         pdf:
           from_: fit
         speclabel: "NNLO"
         dataset_inputs:
              - { dataset: D0WMASY, cfac: [QCD] }

   theoryconfig:  # For generating theory covariance matrix

   theoryid: 163
   
   theoryids:
      from_: scale_variation_theories
   
   point_prescription: '5bar point'

   use_cuts: fromfit
   fit: 190315_ern_nlo_central_163_global

   pdf:
     from_: fit

   dataspecs:
           - theoryid: 163
             speclabel: $(\xi_F,\xi_R)=(1,1)$
             dataset_inputs: *datasets_list_nlo
           - theoryid: 180
             speclabel: $(\xi_F,\xi_R)=(2,2)$
             dataset_inputs: *datasets_list_nlo
           - theoryid: 173
             speclabel: $(\xi_F,\xi_R)=(0.5,0.5)$
             dataset_inputs: *datasets_list_nlo
           - theoryid: 175
             speclabel: $(\xi_F,\xi_R)=(2,0.5)$
             dataset_inputs: *datasets_list_nlo
           - theoryid: 178
             speclabel: $(\xi_F,\xi_R)=(0.5,2)$
             dataset_inputs: *datasets_list_nlo
     
   template: template_test.md

   dataset_report:
     meta: Null
     template_text: |
        ## Testing 5barpt NLO global covariance matrix against NNLO-NLO shift
   actions_:
     - report(main=true, mathjax=True)


The corresponding file ``template_test.md`` is shown below. This will produce
a range of outputs analysing the theory covariance matrix's performance in 
capturing the NNLO-NLO shift.

.. code-block::  md
   :linenos:

   % Theory shift validation test: 5 pt

   Non-zero eigenvalues
   --------------------

   {@theory_covmat_eigenvalues@}

   Efficiency
   ----------

   {@efficiency@}

   Angle between NNLO-NLO shift vector and its component in the theory subspace
   -----------------------------------------------------------------------------------

   {@theta@} 

   Ratio of projectors to eigenvalues
   ----------------------------------
  
   {@projector_eigenvalue_ratio@}

   Condition number of projected matrix
   ------------------------------------

   {@projected_condition_num@}

   Theory $\chi^2$ 
   ---------------
 
   {@validation_theory_chi2@}

   Comparison of NNLO-NLO shift with theory errors from prescription
   -----------------------------------------------------------------

   {@shift_diag_cov_comparison@}

   Eigenvector plots
   -----------------

   {@eigenvector_plot@}

   $\delta_{miss}$ plot
   --------------------

   {@deltamiss_plot@}
