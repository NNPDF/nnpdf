How to plot PDFs
================

Plotting any number of PDFs can be done using ``validphys``.  The runcard below is an example. 

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

    pdfnormalize:
        - normtitle: Absolute  # normalize_to default is None
        - normalize_to: 1      # Specify index in list of PDFs or name of PDF
          normtitle: Ratio

    basespecs:
        - basis: flavour
          basistitle: Flavour basis
        - basis: evolution
          basistitle: Evolution basis

    pdfscalespecs:
        - xscale: log
          xscaletitle: Log
        - xscale: linear
          xscaletitle: Linear
      
    template_text: |
      {@with pdfscalespecs@}
      {@with basespecs@}
      {@with pdfnormalize@}
      {@plot_pdfs@}
      {@plot_pdf_uncertainties@}
      {@endwith@}
      {@endwith@}
      {@endwith@}
  
    actions_:
      - report(main=True)

	  
- The PDFs to be plotted are listed under ``pdfs``, where the PDF ``id`` is specified alongisde a label chosen by the user. If no label is provided, the PDF ``id`` will be displayed on the plots instead.

- To plot one or more PDF sets as a central value only (i.e. without uncertainty band), include them in the list ``pdfs_noband``, where entries in the list are strings corresponding to PDF ``ìds``, integers (starting from 1) of the index of the PDF in ``pdfs``, or a mixture of the two.

- ``show_mc_errors (bool)`` specifies whether to plot 1\[\sigma\] bands in addition to 68\% Monte Carlo uncertainties.

- ``Q`` specifies the plotting scale

- PDFs can be displayed as absolute values (``normalize_to: None``) or as ratios to one of the PDFs (``normalize_to: x``, where ``x`` is the index of the PDF in ``pdfs`` or is the name of the PDF).

- Here the ``basespecs`` namespace specifies the basis for plotting, two popular choices being the ``flavour`` and ``evolution`` bases. For information on custom bases see :ref:`pdfbases`.

- Here the ``pdfscalespecs`` namespace specifies the scale on the x-axis, e.g. ``Log`` or ``Linear``.

- Two actions are called: ``plot_pdfs`` and ``plot_uncertainties``. Their output is fairly self-evident!
 
 
	  
