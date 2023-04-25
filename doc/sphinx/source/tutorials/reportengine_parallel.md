# How to run an analysis in parallel


In this tutorial we will show how to run a `validphys` (or more general a reportengine app analysis) analysis using the parallel resource executor of `reportengine`.

Example 1
--------

Suppose we have the following example report

```yaml
meta:
  title: PDF plot example
  author: Mark N. Costantini
  keywords: [parallel, example]


pdfs:
  - {id: "NNPDF40_nlo_as_01180", label: "4.0 NLO"}
  - {id: "NNPDF40_nnlo_lowprecision", label: "4.0 NNLO low precision"}
  - {id: "NNPDF40_nnlo_as_01180", label: "4.0 NNLO"}


pdfs_noband: ["NNPDF40_nnlo_as_01180"] # Equivalently [3]

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
  {@xscaletitle@} scale
  =====================
  {@with Basespecs@}
  {@basistitle@}
  -------------
  {@with PDFnormalize@}
  {@normtitle@}
  {@plot_pdfs@}
  {@plot_pdf_uncertainties@}
  {@plot_pdfreplicas@}          
  {@endwith@}
  {@endwith@}
  {@endwith@}

actions_:
  - report(main=True)



```







Example 2
---------






