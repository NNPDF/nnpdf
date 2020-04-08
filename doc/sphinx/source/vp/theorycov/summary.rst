========
Summary
========

:Author: Rosalyn Pearson (r.l.pearson@ed.ac.uk)

.. raw:: latex

   \maketitle

.. raw:: latex

   \tableofcontents
   
See the `short
<https://arxiv.org/abs/1905.04311>`_  and `long
<https://arxiv.org/abs/1906.10698>`_ papers for reference.

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

-  As input you need theories for the relevant scale combinations which
   correspond to the prescription. This information is taken from the
   ``scalevariations`` module, which consists of two files:

   #. ``pointprescriptions.yaml``: correspondence between each point prescription
      and the scale combinations that are used to construct it

   #. ``scalevariationtheoryids.yaml``: correspondence between each scale combination
      and a theoryid for a given central theoryid

-  The prescription must be one of :math:`\{3,5,7,9\}`.

-  In the case of 5 theories, you must further specify whether the 5 or
   :math:`\bar{5}` prescription is required. You can do this by
   allocating the flag ``fivetheories`` to ``nobar`` or ``bar`` in the
   runcard.

-  In the case of 7 theories, there are two options. The default is the
   modified prescription that Gavin Salam proposed. To use the original
   prescription instead, specify ``seventheories: original`` in the runcard.

-  Currently the renormalisation scales are correlated within each
   process type. These process types are categorised as {DIS CC, DIS NC,
   Drell-Yan, Jets, Top}. 

-  Outputs include tables and heat plots of theoretical and combined
   (theoretical + experimental) covariance matrices, comparisons of
   theoretical and experimental errors, and plots and tables of
   :math:`\chi^2` values.

-  Various validation outputs also exist, including tables of eigenvalues, 
   plots of eigenvectors and shift vs theory comparisons.
