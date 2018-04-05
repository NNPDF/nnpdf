%NNPDF Report for fit {@ current fit @}

Dataset properties
------------------
{@current datasets_properties_table@}

Distances
------------------
{@with normalize::basespecs::pdfscalespecs@}
{@plot_pdfdistances@}
{@plot_pdfvardistances@}
{@endwith@}

PDF arc-lengths
---------------
{@basespecs plot_arc_lengths@}

PDF plots
---------
[Detailed plots]({@pdf_report report@})

Training lengths
----------------
{@fits plot_training_length@}

Training validation
-------------------
{@fits plot_training_validation@}

$\chi^2$ by experiment
----------------------
{@plot_fits_experiments_chi2@}

$\chi^2$ by dataset comparisons
-------------------------------
### Plot
{@plot_fits_datasets_chi2@}
### Table
{@fits_chi2_table(show_total=true)@}

$\phi$ by experiment
--------------------
{@with dataspecs@}
{@plot_phi_pdfs@}
{@endwith@}

{@with dataspecs@}
### {@fit@}
{@plot_phi@}
{@endwith@}

Experiment plots
---------------
{@with matched_datasets_from_dataspecs@}
[Detailed plots for dataset ' {@dataset_name@} ']({@dataset_report report@})
{@endwith@}
