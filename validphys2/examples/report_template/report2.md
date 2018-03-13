%NNPDF Report for fit {@ current fit @}

PDF plots
---------

{@with pdfnormalize@}
## {@normtitle@}
{@with basespecs@}
### {@basistitle@}
{@plot_pdfs@}
{@endwith@}
{@endwith@}

$\chi^2$
--------
{@fits_chi2_table@}


Experiment plots
---------------

{@with matched_datasets_from_dataspecs@}
## Data theory comparison for {@dataset_name@}        
### Absolute
{@ plot_fancy_dataspecs @}
### Normalized
{@ datanorm plot_fancy_dataspecs @}
{@endwith@}
