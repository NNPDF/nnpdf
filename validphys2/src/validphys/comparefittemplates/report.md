%NNPDF Report for fit {@ current fit @}

Fit summary 
------------------
{@ summarise_fits @}

Dataset properties
------------------
{@current datasets_properties_table@}

Distances
------------------
{@with normalize::basespecs::pdfscalespecs::distspecs@}
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
{@plot_fits_experiments_phi@}

Experiment plots
---------------
{@with matched_datasets_from_dataspecs@}
[Detailed plots for dataset ' {@dataset_name@} ']({@dataset_report report@})
{@endwith@}

Positivity
----------
{@with matched_positivity_from_dataspecs::dataspecs@}
{@plot_positivity@}
{@endwith@}

Dataset differences and cuts
----------------------------
{@print_dataset_differences@}
{@print_different_cuts@}
