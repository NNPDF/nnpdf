============================================================
Documentation for theory covariance module in validphys2
============================================================

:Author: Rosalyn Pearson (r.l.pearson@ed.ac.uk)

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

   #. ``tests.py``: actions for validating the covariance matrices against
      the NNLO-NLO shift

-  Theoretical covariance matrices are built according to the various prescriptions.

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
   Drell-Yan, Jets, Top}. 

-  Outputs include tables and heat plots of theoretical and combined
   (theoretical + experimental) covariance matrices, comparisons of
   theoretical and experimental errors, and plots and tables of
   :math:`\chi^2` values.

-  Various validation outputs also exist, including tables of eigenvalues, 
   plots of eigenvectors and shift vs theory comparisons.

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

-  You must also provide the relevant c-factors (EWK, QCD ...).


Outputs
=======

Below is the template corresponding to the runcard in the previous
section. This will produce a report with the important features of the
theory covariance module. These fall into three broad categories:

#. **Matrices and plots of matrices**. Heatmap plots, diagonal element
   plots comparing experimental and theoretical errors.

#. :math:`\chi^2` **values**. Total :math:`\chi^2`\ s and per dataset
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

Examples
========

Covariance matrix report
------------------------

The following is an example of how to produce a report detailing the covariance
matrices. In this case, the 3 point prescription is shown, for a global data
at NLO.

You need to provide the central theory under the ``default_theory`` flag, 
corresponding to :math:`(\mu_F, \mu_R) = (0,0)`,
which for NLO is theory 163. Under ``theoryids`` you need to provide all the
relevant theories, as outlined above.

``dataspecs`` associates a chosen label (``speclabel``) with each of the theory
choices. This details what scale variation the theory corresponds to.

Here the cuts and PDF are taken from the central NLO scale-varied fit.

You must also list all the experiments you wish to include, along with any 
relevant c-factors. 

.. code-block::  yaml
   :linenos:
   
   meta:
      author: Rosalyn Pearson
      keywords: [theory uncertainties, 3-point]
      title: NLO 3-point variations for 5 process types - DIS CC, DIS NC, DY, Top, Jets

   default_theory:
      - theoryid: 163

   theoryids:
      - 163
      - 180
      - 173

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

   experiments:
     - experiment: NMC
       datasets:
         - dataset: NMCPD
         - dataset: NMC
     - experiment: SLAC
       datasets:
         - dataset: SLACP
         - dataset: SLACD
     - experiment: BCDMS
       datasets:
         - dataset: BCDMSP
         - dataset: BCDMSD
     - experiment: NTVDMN
       datasets:
         - dataset: NTVNUDMN
         - dataset: NTVNBDMN
     - experiment: CHORUS
       datasets:
         - dataset: CHORUSNU
         - dataset: CHORUSNB
     - experiment: HERAF2CHARM
       datasets:
         - dataset: HERAF2CHARM
     - experiment: HERACOMB
       datasets:
         - dataset: HERACOMBNCEM 
         - dataset: HERACOMBNCEP460
         - dataset: HERACOMBNCEP575
         - dataset: HERACOMBNCEP820
         - dataset: HERACOMBNCEP920
         - dataset: HERACOMBCCEM 
         - dataset: HERACOMBCCEP 
     - experiment: ATLAS
       datasets:
         - dataset: ATLASWZRAP36PB
         - dataset: ATLASZHIGHMASS49FB
         - dataset: ATLASLOMASSDY11EXT
         - dataset: ATLASWZRAP11
         - dataset: ATLAS1JET11
         - dataset: ATLASZPT8TEVMDIST
         - dataset: ATLASZPT8TEVYDIST
         - dataset: ATLASTTBARTOT
         - dataset: ATLASTOPDIFF8TEVTRAPNORM
     - experiment: CMS
       datasets:
         - dataset: CMSWEASY840PB
         - dataset: CMSWMASY47FB
         - dataset: CMSDY2D11
         - dataset: CMSWMU8TEV
         - { dataset: CMSZDIFF12, cfac: [NRM] }
         - dataset: CMSJETS11
         - dataset: CMSTTBARTOT
         - dataset: CMSTOPDIFF8TEVTTRAPNORM
     - experiment: LHCb
       datasets:
         - dataset: LHCBZ940PB
         - dataset: LHCBZEE2FB
         - { dataset: LHCBWZMU7TEV, cfac: [NRM] }
         - { dataset: LHCBWZMU8TEV, cfac: [NRM] }
     - experiment: CDF
       datasets:
         - dataset: CDFZRAP
     - experiment: D0
       datasets:
         - dataset: D0ZRAP
         - dataset: D0WEASY
         - dataset: D0WMASY

   template: template.md

   dataset_report:
      meta: Null
      template_text: |
         ## Scale variations as a function of the kinematics for {@dataset_name@}
         {@plot_fancy_dataspecs@}

   actions_:
     - report(main=true) 


The corresponding template file is ``template.md``, shown below. This will produce
a comprehensive set of plots and tables describingn the covariance matrices. 

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

   Experimental $\chi^2$
   ---------------------
   {@with default_theory@}
      {@total_experiments_chi2@}

   Total (exp. + th.) $\chi^2$
   ---------------------------
      {@chi2_impact_custom@}

   Experimental $\chi^2$ by dataset
   --------------------------------
      {@experiments_chi2_table@}

   Total (exp. + th.) $\chi^2$ by dataset
   --------------------------------------
      {@experiments_chi2_table_theory@}

   $\chi^2$ including only diagonal theory elements
   ------------------------------------------------
      {@chi2_diag_only@}

   Impact of theory covariance matrix on $\chi^2$s 
   -----------------------------------------------
      {@plot_datasets_chi2_theory@}
   {@endwith@}

   Scale variations as a function of the kinematics
   ------------------------------------------------
   {@with matched_datasets_from_dataspecs@}
      [Plots for {@dataset_name@}]({@dataset_report report@})
   {@endwith@}


Validation report
----------------- 

Here is an example of a runcard for a report validating the theory covariance
matrix against the NNLO-NLO shift. In this case the 5 point prescription is chosen,
and Drell-Yan experiments only are considered.

Note that as we are dealing with 5 theories, we need to set the ``fivetheories``
flag, which in this case is set to ``nobar``. This must be used in conjuction
with the correct ``theoryids`` and ordering of ``theoryids`` in order not to throw 
an error. 

The flag ``orthonormalisation`` corresponds to the method used to orthonormalise 
the basis vectors of the theory covariance matrix. There are three choices:

#. QR decomposition (choose this by default), with the flag ``qr``

#. Singular value decompostion, with the flag ``svd``

#. An in-built Gram-Schmidt orthonormalisation, with the flag ``gs``.

``_experiments_list_nlo`` is a list of all the experiments to be included at NLO.
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

   meta:
       title: Theory shift validation test, 5 point, DY-only, QR
       author: Rosalyn Pearson
       keywords: [test, theory uncertainties, eigenvalues, 5 point]

   fivetheories: nobar

   orthonormalisation: qr

   theoryid: 163

   fit: 190315_ern_nlo_central_163_global

   pdf:
     from_: fit

   _experiments_list_nlo: &experiments_list_nlo
     - experiment: NMC
       datasets:
         - dataset: NMCPD
         - dataset: NMC
     - experiment: SLAC
       datasets:
         - dataset: SLACP
         - dataset: SLACD
     - experiment: BCDMS
       datasets:
         - dataset: BCDMSP
         - dataset: BCDMSD
     - experiment: NTVDMN
       datasets:
         - dataset: NTVNUDMN
         - dataset: NTVNBDMN
     - experiment: CHORUS
       datasets:
         - dataset: CHORUSNU
         - dataset: CHORUSNB
     - experiment: HERAF2CHARM
       datasets:
         - dataset: HERAF2CHARM
     - experiment: HERACOMB
       datasets:
         - dataset: HERACOMBNCEM 
         - dataset: HERACOMBNCEP460
         - dataset: HERACOMBNCEP575
         - dataset: HERACOMBNCEP820
         - dataset: HERACOMBNCEP920
         - dataset: HERACOMBCCEM 
         - dataset: HERACOMBCCEP 
     - experiment: ATLAS
       datasets:
         - dataset: ATLASWZRAP36PB
         - dataset: ATLASZHIGHMASS49FB
         - dataset: ATLASLOMASSDY11EXT
         - dataset: ATLASWZRAP11
         - dataset: ATLAS1JET11
         - dataset: ATLASZPT8TEVMDIST
         - dataset: ATLASZPT8TEVYDIST
         - dataset: ATLASTTBARTOT
         - dataset: ATLASTOPDIFF8TEVTRAPNORM
     - experiment: CMS
       datasets:
         - dataset: CMSWEASY840PB
         - dataset: CMSWMASY47FB
         - dataset: CMSDY2D11
         - dataset: CMSWMU8TEV
         - { dataset: CMSZDIFF12, cfac: [NRM] }
         - dataset: CMSJETS11
         - dataset: CMSTTBARTOT
         - dataset: CMSTOPDIFF8TEVTTRAPNORM
     - experiment: LHCb
       datasets:
         - dataset: LHCBZ940PB
         - dataset: LHCBZEE2FB
         - { dataset: LHCBWZMU7TEV, cfac: [NRM] }
         - { dataset: LHCBWZMU8TEV, cfac: [NRM] }
     - experiment: CDF
       datasets:
         - dataset: CDFZRAP
     - experiment: D0
       datasets:
         - dataset: D0ZRAP
         - dataset: D0WEASY
         - dataset: D0WMASY

   shiftconfig:

      use_cuts: fromfit
      fit: 190315_ern_nlo_central_163_global

      theoryid: 163

      dataspecs:
          - theoryid: 163
            pdf:
              from_: fit
            speclabel: "NLO"
            experiments:
                - experiment: ATLAS
                  datasets:
                     - dataset: ATLASWZRAP36PB
                     - dataset: ATLASZHIGHMASS49FB
                     - dataset: ATLASLOMASSDY11EXT
                     - dataset: ATLASWZRAP11
                     - dataset: ATLASZPT8TEVMDIST
                     - dataset: ATLASZPT8TEVYDIST
                - experiment: CMS
                  datasets:
                     - dataset: CMSWEASY840PB
                     - dataset: CMSWMASY47FB
                     - dataset: CMSDY2D11
                     - dataset: CMSWMU8TEV
                     - { dataset: CMSZDIFF12, cfac: [NRM] }
                - experiment: LHCb
                  datasets:
                     - dataset: LHCBZ940PB
                     - dataset: LHCBZEE2FB
                     - { dataset: LHCBWZMU7TEV, cfac: [NRM] }
                     - { dataset: LHCBWZMU8TEV, cfac: [NRM] }
                - experiment: CDF
                  datasets:
                     - dataset: CDFZRAP
                - experiment: D0
                  datasets:
                     - dataset: D0ZRAP
                     - dataset: D0WEASY
                     - dataset: D0WMASY
          - theoryid: 166
            pdf:
              from_: fit
            speclabel: "NNLO"
            experiments:
                - experiment: ATLAS
                  datasets:
                     - { dataset: ATLASWZRAP36PB, cfac: [QCD]}
                     - { dataset: ATLASZHIGHMASS49FB, cfac: [QCD] }
                     - { dataset: ATLASLOMASSDY11EXT, cfac: [QCD] }
                     - { dataset: ATLASWZRAP11, cfac: [QCD] }
                     - { dataset: ATLASZPT8TEVMDIST, cfac: [QCD], sys: 10 }
                     - { dataset: ATLASZPT8TEVYDIST, cfac: [QCD], sys: 10 }
                - experiment: CMS
                  datasets:
                     - { dataset: CMSWEASY840PB, cfac: [QCD] }
                     - { dataset: CMSWMASY47FB, cfac: [QCD]}
                     - { dataset: CMSDY2D11, cfac: [QCD] }
                     - { dataset: CMSWMU8TEV, cfac: [QCD] }
                     - { dataset: CMSZDIFF12, cfac: [QCD, NRM], sys: 10 }
                - experiment: LHCb
                  datasets:
                     - { dataset: LHCBZ940PB, cfac: [QCD] }
                     - { dataset: LHCBZEE2FB, cfac: [QCD] }
                     - { dataset: LHCBWZMU7TEV, cfac: [QCD, NRM] }
                     - { dataset: LHCBWZMU8TEV, cfac: [QCD, NRM] }
                - experiment: CDF
                  datasets:
                     - { dataset: CDFZRAP, cfac: [QCD] }
                - experiment: D0
                  datasets:
                     - { dataset: D0ZRAP, cfac: [QCD] }
                     - { dataset: D0WEASY, cfac: [QCD] }
                     - { dataset: D0WMASY, cfac: [QCD] }

   theoryconfig:

      theoryids:
         - 163
         - 177
         - 176
         - 179
         - 174

      theoryid: 163

      use_cuts: fromfit
      fit: 190315_ern_nlo_central_163_global

      pdf:
        from_: fit

      dataspecs:
              - theoryid: 163
                speclabel: $(\xi_F,\xi_R)=(1,1)$
                experiments: *experiments_list_nlo
              - theoryid: 177
                speclabel: $(\xi_F,\xi_R)=(2,1)$
                experiments: *experiments_list_nlo
              - theoryid: 176
                speclabel: $(\xi_F,\xi_R)=(0.5,1)$
                experiments: *experiments_list_nlo
              - theoryid: 179
                speclabel: $(\xi_F,\xi_R)=(1,2)$
                experiments: *experiments_list_nlo
              - theoryid: 174
                speclabel: $(\xi_F,\xi_R)=(1,0.5)$
                experiments: *experiments_list_nlo

   template: ../../template_test.md

   dataset_report:
      meta: Null
      template_text: |
         ## Testing 5pt NLO global covariance matrix against NNLO-NLO shift
   actions_:
     - report(main=true, mathjax=True)


The corresponding file ``template_test.md`` is shown below. This will produce
a range of outputs analysing the theory covariance matrix's performance in 
capturing the NNLO-NLO shift.

.. code-block::  md
   :linenos:

   % Theory shift validation test: 9 pts

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

.

