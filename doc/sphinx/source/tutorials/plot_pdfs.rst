.. _pdfplots:

How to plot PDFs, distances and luminosities
============================================

Plotting any number of PDFs can be done using ``validphys``.  There are several
kinds of plots which can be made using the actions in the module ``pdfplots.py``.
The runcards in this section can be found
`here <https://github.com/NNPDF/nnpdf/tree/master/validphys2/examples>`_.

Runcard basics
--------------
First we review key aspects of the runcards.

+----------------------------------------+-------------------------------------------------------------------------------------------------------+
| Element                                | Description                                                                                           |
+========================================+=======================================================================================================+
| ``pdfs``                               | PDFs to be plotted. For each PDF, ``id`` is specified alongside a label chosen by the user.           |
|                                        | If no label is provided, the PDF ``id`` will be displayed on the plots instead.                       |
+----------------------------------------+-------------------------------------------------------------------------------------------------------+
| ``pdfs_noband``                        | One or more PDF sets to be plotted as a central value only (i.e. without an uncertainty band).        |
|                                        | A list of strings corresponding to PDF ``ids``, integers (starting from 1) of the index of the        |
|                                        | PDF in ``pdfs``, or a mixture of the two.                                                             |
+----------------------------------------+-------------------------------------------------------------------------------------------------------+
| ``show_mc_errors (bool)``              | Specifies whether to plot 1\\(\\sigma\\) bands in addition                                            |
|                                        | to 68\% Monte Carlo uncertainties.                                                                    |
+----------------------------------------+-------------------------------------------------------------------------------------------------------+
| ``Q``                                  | Plotting scale.                                                                                       |
+----------------------------------------+-------------------------------------------------------------------------------------------------------+
| ``normalize_to``                       | If specified, PDFs are displayed as ratios to one of the PDFs (``normalize_to: x``, where ``x``       |
|                                        | is the index (starting from 1) of the PDF in ``pdfs`` or is the name of the PDF).                     |
+----------------------------------------+-------------------------------------------------------------------------------------------------------+
| ``basis``                              | The basis for plotting, two popular choices being the ``flavour`` and ``evolution`` bases.            |
|                                        | For information on custom bases see :ref:`here <pdfbases>`.                                           |
+----------------------------------------+-------------------------------------------------------------------------------------------------------+
| ``xscale``                             | The scale on the x-axis, e.g. ``Log`` or ``Linear``.                                                  |
+----------------------------------------+-------------------------------------------------------------------------------------------------------+
| ``flavours``                           | Flag for plotting specific flavours only. The chosen flavours should be provided as a list, where the |
|                                        | name of flavour or PDG value (listed in :ref:`this section <pdgflavs>`) can be used.                  |
+----------------------------------------+-------------------------------------------------------------------------------------------------------+
| ``xmin``, ``xmax``, ``ymin``, ``ymax`` | Plot axes limits.                                                                                     |
+----------------------------------------+-------------------------------------------------------------------------------------------------------+


Plotting PDFs, uncertainties and replicas
-----------------------------------------

.. code:: yaml

    meta:
      title: PDF plot example
      author: Rosalyn Pearson
      keywords: [example]

    pdfs:
        - {id: "NNPDF31_nlo_as_0118", label: "3.1 NLO"}
        - {id: "NNPDF31_nnlo_as_0118", label: "3.1 NNLO"}
        - {id: "NNPDF31_nnlo_as_0118_DISonly", label: "3.1 DIS only NNLO"}

    pdfs_noband: ["NNPDF31_nnlo_as_0118_DISonly"] # Equivalently ["3"]

    show_mc_errors: True

    Q: 10

    PDFnormalize:
        - normtitle: Absolute  # normalize_to default is None
        - normalize_to: 1      # Specify index in list of PDFs or name of PDF
          normtitle: Ratio

    Basespecs:
        - basis: flavour
          basistitle: Flavour basis
        - basis: evolution
          basistitle: Evolution basis

    PDFscalespecs:
        - xscale: log
          xscaletitle: Log
        - xscale: linear
          xscaletitle: Linear

    template_text: |
      {@with PDFscalespecs@}
      {@with Basespecs@}
      {@with PDFnormalize@}
      {@plot_pdfs@}
      {@plot_pdf_uncertainties@}
      {@plot_pdfreplicas@}
      {@endwith@}
      {@endwith@}
      {@endwith@}

    actions_:
      - report(main=True)

- Here the ``Basespecs`` namespace specifies the two bases to be plotted.
- Here the ``PDFscalespecs`` namespace specifies two x-scales to be plotted.
- The actions are called: ``plot_pdfs``, ``plot_uncertainties`` and
  ``plot_pdfreplicas``. Their output is fairly self-evident!

Plotting PDF distances
----------------------

.. code:: yaml

	meta:
	    title: I didn't change the title
	    keywords: [Guilty]
	    author: Lazy Person

	pdfs:
	    - NNPDF31_nlo_as_0118
	    - NNPDF31_nnlo_as_0118_DISonly

	pdf: NNPDF31_nlo_as_0118

	First:
	    Q: 2
	    flavours: [up, down, gluon, 4]

	Second:
	    Q: 100
	    scale: linear
	    flavours: [up, down, gluon, 4]

	normalize_to: 1

	template_text: |

	  Log scale, low Q
	  -----------
	  {@First plot_pdfdistances@}
	  {@First plot_pdfvardistances@}

	  Linear scale, high Q
	  -----------
	  {@Second plot_pdfdistances@}
	  {@Second plot_pdfvardistances@}

	actions_:
	  - report(main=true)

- The actions ``plot_pdfdistances`` and ``plot_pdfvardistances`` plot the
  distances of the PDFs and the variances of these distances with respect to the
  PDF specalised by ``normalize_to``.

Plotting PDF flavours on the same axis
--------------------------------------
.. code:: yaml

	meta:
	  title: PDF flavours plot example
	  author: Rosalyn Pearson
	  keywords: [example]

	pdf:  {id: "NNPDF31_nlo_as_0118", label: "3.1 NLO"}

	Q: 10

	basis: pdg      # [g/10, u_v, d_v, s, ubar, dbar, c] plots well on same axis
	xmin: 0.002

	ymin: 0
	ymax: 0.6

	xscale: log

	template_text: |
	  {@plot_flavours@}

	actions_:
	  - report(main=True)

- ``plot_flavours`` is the action used to plot PDF flavours on the same axes.
- Note that the ``basis`` has been set to ``pdg``, which is a configuration that
  plots the various flavours well on the same axis as the gluon PDF is divided
  by 10. More on PDF bases :ref:`here <pdfbases>`.

Luminosity plots
----------------
.. code:: yaml

	meta:
	  title: Luminosity plot example
	  author: Rosalyn Pearson
	  keywords: [example]

	pdfs:
	  - {id: "NNPDF31_nlo_as_0118", label: "3.1 NLO"}
	  - {id: "NNPDF31_nnlo_as_0118", label: "3.1 NNLO"}
	  - {id: "NNPDF31_nnlo_as_0118_DISonly", label: "3.1 DIS only NNLO"}

	pdf: {id: "NNPDF31_nlo_as_0118", label: "3.1 NLO"}

	sqrts: 13000 # GeV

	lumi_channel: "gg" # one of [gg, gq, qqbar, qq, ddbar, uubar, ssbar,
	                   #         ccbar, bbbar, dubar, udbar, scbar, csbar, pp, gp]

	y_cut: 2.5
			   
	PDFscalespecs:
	  - xscale: log
	    Xscaletitle: Log
	  - xscale: linear
	    Xscaletitle: Linear

	template_text: |
	  {@with PDFscalespecs@}
	  {@Xscaletitle@} scale
	  =====================
	  {@plot_lumi1d@}
	  {@plot_lumi1d_uncertainties@}
	  {@plot_lumi2d@}
	  {@plot_lumi2d_uncertainty@}
	  {@endwith@}

	actions_:
	  - report(main=True)

- Luminosity plots can be made using the actions in the above runcard.
- The 1D plots show the luminosity as a function of the invariant mass of the
  final state *X*, while the 2D plots show it as a function of the invariant
  mass *X* and the rapidity *y* of the final state.
- A choice of ``lumi_channel`` must be provided, which is a string in one of
  [gg, gq, qqbar, qq, ddbar, uubar, ssbar, ccbar, bbbar, dubar, udbar, scbar,
  csbar, pp, gp]. Note that p represents the photon and q and qbar represent all
  six quark flavours, i.e. [u, d, s, c, b, t] or [ubar, dbar, sbar, cbar, bbar,
  tbar], respectively. The exact definitions of the channels in the code can be
  found `here <https://github.com/NNPDF/nnpdf/blob/c20f1892767632f4764ada12bc106c04d5b739d4/validphys2/src/validphys/gridvalues.py>`_.
  Note further that a set of channels can be given by
  specifying ``lumi_channels`` followed by a list of channels.
- The square root of centre of mass energy, √s, in GeV must also
  be provided via ``sqrts``. This is instead of ``Q``.
- The options ``ymin`` and ``ymax`` can be supplied to control the vertical
  range of the plot, while ``mxmin`` and ``mxmax`` (in GeV) control the
  horizontal axis of invariant masses.
- The option ``y_cut`` can be supplied to define a rapidity cut
  ( ``|y|<y_cut``) on the integration range of 1D plots.
- Other options for the band plots, such as ``pdfs_noband`` and
  ``show_mc_errors`` work for ``plot_lumi1d``.
