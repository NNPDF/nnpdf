 .. _vptheorycov-index:
 
The `theorycovariance` module
===============================

The ``theorycovariance`` module deals with constructing, testing and 
outputting theory covariance matrices (covmats). Primarily, it is concerned
with scale variation covariance matrices used to model missing higher order
uncertainties. See the `short
<https://arxiv.org/abs/1905.04311>`_  and `long
<https://arxiv.org/abs/1906.10698>`_ NNPDF papers for in-depth information.

Summary
-------

-  The module of ``validphys2`` which deals with computation and
   interpretation of theoretical covariance matrices can be found in
   ``nnpdf/validphys2/src/validphys/theorycovariance/``, which consists
   of three files:

   #. ``construction.py``: deals with construction of covariance
      matrices and associated quantities

   #. ``output.py``: plots and tables

   #. ``tests.py``: actions for validating the covariance matrices against
      the NNLO-NLO shift

-  Theoretical covariance matrices are built according to the various prescriptions
   in :ref:`prescrips`. 
 
-  The prescription must be one of 3 point, 5 point, 5bar point, 7 point, 7original point or 9 point. You can specify
   this using ``point_prescription: "x point"`` in the runcard. The translation of this flag 
   into the relevant ``theoryids`` is handled by the ``scalevariations`` module in ``validphys``.

-  As input you need theories for the relevant scale combinations which
   correspond to the prescription. This information is taken from the
   ``scalevariations`` module, which consists of two files:

   #. ``pointprescriptions.yaml``: correspondence between each point prescription
      and the scale combinations that are used to construct it

   #. ``scalevariationtheoryids.yaml``: correspondence between each scale combination
      and a theoryid for a given central theoryid

-  Renormalisation scales should be correlated within each
   process type. These process types are categorised as {DIS CC, DIS NC,
   Drell-Yan, Jets, Top}. 

-  :ref:`Outputs <thcov_outputs>` includes tables and heat plots of theoretical and combined
   (theoretical + experimental) covariance matrices, comparisons of
   theoretical and experimental errors, and plots and tables of
   :math:`\chi^2` values.

-  Various :ref:`testing <vptheorycov-tests>` outputs also exist, including tables of eigenvalues, 
   plots of eigenvectors and shift vs theory comparisons.
   

More information
-----------------
.. toctree::
   :maxdepth: 1

   ./runcard_layout.rst
   ./outputs.rst
   ./point_prescrip.rst
   ./tests.rst
   ./examples.rst
