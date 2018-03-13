%NNPDF Report for fit {@ current::fit @}

#PDF plots
{@with pdfnormalize@}
## {@normtitle@}
{@with basespecs@}
### {@basistitle@}
{@plot_pdfs@}
{@endwith@}
{@endwith@}

# $\chi^2$
{@fits_chi2_table@}

# Comparison with data
{@with reference::fitcontext::experiments@}
###{@ experiment @}
{@with experiment@}
[Detailed plots for dataset ' {@dataset@} ']({@dataset_report report @})
{@ endwith @}
{@ endwith @}