%NNPDF Report for fit {@ closure fit @}

{@ with closure::fitunderlyinglaw @}

Fit summary 
------------------
{@ summarise_fits @}

Closure test estimators
-----------------------
## $\Delta_{\chi^{2}}$ by experiment
{@ plot_delta_chi2 @}
## Bias by experiment
{@ plot_biases @}
## $\phi$ by experiment
{@closure plot_phi@}

Dataset properties
------------------
{@closure datasets_properties_table@}

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

Experiment plots
---------------
{@with matched_datasets_from_dataspecs@}
[Detailed plots for dataset ' {@dataset_name@} ']({@dataset_report report@})
{@endwith@}

{@endwith@}
