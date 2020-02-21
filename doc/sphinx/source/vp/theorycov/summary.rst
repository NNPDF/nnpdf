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
