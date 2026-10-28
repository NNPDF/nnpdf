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

-  The prescription for scale variations must be one of ``3(f, r) point``,
   ``5(bar) point``, ``7(original) point`` or ``9 point``, see
   :ref:`definitions <prescrips>`. Beside this, the following prescriptions
   are available for IHOU variations and for alphas variations. You can
   specify them using ``point_prescription: "x point"`` in the runcard. The
   translation of this flag into the relevant ``theoryids`` is handled by the
   ``scalevariations`` module in ``validphys``.

   IHOU variations
      ``ad ihou``
         Varies the IHOU Anomalous Dimensions in EKO.

      ``dis ihou``
         Varies the IHOU Coefficient Functions in YADISM.

      ``3pt missing``
         3 point renormalisation scale variations for a subset of data.

      ``3pt hadronic``
         3 point renormalisation scale variations for hadronic datasets.

      ``fhmruvv ihou``
         Varies the FHMRUVV splitting functions in EKO.

   Alphas variations
      ``alphas 0117 0119 0115``
         Varies \alpha_s central at 0.117 from 0.115 to 0.115.

      ``alphas 0118 0120 0116``
         Varies \alpha_s central at 0.118 from 0.120 to 0.116.

      ``alphas 0119 0121 0117``
         Varies \alpha_s central at 0.119 from 0.121 to 0.117.

      ``alphas 0119 0122 0116``
         Varies \alpha_s central at 0.119 from 0.121 to 0.116.

      ``alphas 0119 0123 0115``
         Varies \alpha_s central at 0.119 from 0.123 to 0.115.

      ``alphas 0120 0122 0118``
         Varies \alpha_s central at 0.120 from 0.122 to 0.118.

      ``alphas 0121 0123 0119``
         Varies \alpha_s central at 0.121 from 0.123 to 0.119.

-  You may use multiple prescriptions in the theory covariance matrix, by specifying multiple prescriptions at once. 
   E.g. ``point_prescriptions: ["fhmruvv ihou", "dis ihou", "7 point"]`` will vary the weights using  three prescriptions.

-  As input you need theories for which the relevant scale combinations that
   correspond to the prescription actually exist. This information is taken from the
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
   ./power_corrections.rst
   ./tests.rst
   ./examples.rst
