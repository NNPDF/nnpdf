%NNPDF Report for fit {@ closures fit @}

{@ with closures::fitunderlyinglaw @}

Fit summary 
------------------
{@ summarise_fits @}

Closures test estimators
-----------------------
{@ plot_delta_chi2 @}
{@ plot_biases @}

Dataset properties
------------------
{@closures datasets_properties_table@}

Distances
------------------
{@with normalize::basespecs::pdfscalespecs::distspecs@}
{@plot_pdfdistances@}
{@plot_pdfvardistances@}
{@endwith@}

PDF Uncertainty
---------------
[Pdf uncertainties]({@pdf_uncertainty_report report@})

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
{@with fits@}
{@plot_phi@}
{@endwith@}

{@endwith@}
