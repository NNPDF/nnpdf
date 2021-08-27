.. _bundled-sets:

Bundle PDFs with :math:`\alpha_s` replicas
==========================================

Using ``validphys`` it is possible to produce a bundled
PDF set which accounts for a combined PDF + :math:`\alpha_s`
uncertainty. The procedure to generate such an LHAPDF set
is to simply take a baseline PDF set and append the central
replica from a list of :math:`\alpha_s` variation fits.

The action to leverage this is
:py:func:`validphys.replica_selector.alpha_s_bundle_pdf`. We
specify the baseline PDF as the ``pdf`` key within the runcard
and the ``pdfs`` list specifies the :math:`\alpha_s` fits that
are to be used. In the following example, the ``NNPDF31_nnlo_as_0118``
PDF set is used as baseline and we append the central replica from
``NNPDF31_nnlo_as_0117`` and ``NNPDF31_nnlo_as_0119``.

.. code-block :: yaml

   pdf: NNPDF31_nnlo_as_0118

   pdfs:
    - NNPDF31_nnlo_as_0117
    - NNPDF31_nnlo_as_0119

   actions_:
    - alpha_s_bundle_pdf

Executing this runcard with ``validphys`` produces the bundled PDF set
in the output folder, which by default will be name the same as the baseline,
except for the ``_pdfas`` suffix being appended:

.. code-block ::

  output/
  ├── NNPDF31_nnlo_as_0118_pdfas
  │   ├── NNPDF31_nnlo_as_0118_pdfas.info
  │   ├── NNPDF31_nnlo_as_0118_pdfas_0000.dat
  │   ├── NNPDF31_nnlo_as_0118_pdfas_0001.dat
  │   ├── NNPDF31_nnlo_as_0118_pdfas_0002.dat

The ``.info`` file now has ``replicas+as`` as the ``ErrorType`` and the
``NumMembers`` key has now been updated to reflect the additional replicas
(in this example it has been incremented by two). The additional replicas
also have the ``AlphaS_MZ`` key prepended at the top of the file so as to
keep track of what PDF set they originated from.

The optional ``target_name`` key can be provided in the runcard so as to
specify the name of the resulting PDF. The following runcard will generate
the same bundled PDF as before, but with the name ``bundled_pdf``:

.. code-block :: yaml

  pdf: NNPDF31_nnlo_as_0118

  pdfs:
    - NNPDF31_nnlo_as_0117
    - NNPDF31_nnlo_as_0119

  target_name: bundled_pdf

  actions_:
    - alpha_s_bundle_pdf


.. code-block ::

  output/
  ├── bundled_pdf
  │   ├── bundled_pdf.info
  │   ├── bundled_pdf_0000.dat
  │   ├── bundled_pdf_0001.dat
  │   ├── bundled_pdf_0002.dat

.. note ::

  Despite adding additional replicas, the central replica in the bundled
  PDF set is **not** recomputed: it is identical to that of the baseline.
