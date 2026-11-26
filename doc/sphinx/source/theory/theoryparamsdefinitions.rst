.. _th_parameter_definitions:

Theory parameter definitions
============================

Table describing the variables used to specify theory choices in :code:`nnpdf`. Here they are shown
directly as the fields of the runcards in :code:`theory_cards`.
The `Type` column refers directly to data type of the field.

.. csv-table::
   :file: ../theory_params.csv
   :header-rows: 1

QED coupling
------------

The parameter :code:`alphaqed` is the running QED coupling in the :math:`\overline{\rm MS}` scheme,
:math:`\alpha_{\rm em}(Q_{\rm ref})`. It is the starting value of the running of
:math:`\alpha_{\rm em}` in the QED evolution of the PDFs (EKO) and in the computation of the
photon PDF (FiatLux). For :math:`Q_{\rm ref} = M_Z` the value is
:math:`\hat{\alpha}^{(5)}(M_Z)^{-1} = 127.930 \pm 0.008` (for :math:`\alpha_s(M_Z) = 0.1187`), from
Section 10, "Electroweak Model and Constraints on New Physics", of the
`Review of Particle Physics 2024 <https://pdg.lbl.gov/2024/reviews/rpp2024-rev-standard-model.pdf>`_.
The 2026 edition no longer quotes it, because its electroweak fit moved to codes that work in the
on-shell scheme.

:code:`alphaqed` is *not* the coupling used in the computation of the hard cross sections.
These are computed (with :code:`pinefarm`) in the :math:`G_\mu` scheme, where the electroweak
coupling is fixed and derived from :code:`GF`, :code:`MW` and :code:`MZ`,

.. math::

   \alpha_{G_\mu} = \frac{\sqrt{2}}{\pi} G_F M_W^2 \left(1 - \frac{M_W^2}{M_Z^2}\right) \simeq \frac{1}{132.1} \, .

This coupling does not run, and it differs from the running coupling at :math:`M_Z` by about 3%.
As a safeguard, a theory with :code:`QED > 0` whose :code:`alphaqed` is within 2% of
:math:`\alpha_{G_\mu}` is rejected.
