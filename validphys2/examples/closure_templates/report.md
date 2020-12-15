%NNPDF Report for fit {@ closure fit @}

{@ with closure::fitunderlyinglaw @}

Fit summary 
------------------

We are looking at:

  - {@ closure fit @} (`{@ closure fit_id @}`): {@ closure description @}

{@ summarise_fits @}

Underlying PDF Summary
----------------------
{@ summarise_closure_underlying_pdfs @}

Closures test estimators
-----------------------
## Bias by experiment
### Plot

Errorbars from performing bootstrap sample upon replicas

{@ plot_fits_bootstrap_bias @}
### Table
{@ fits_bootstrap_bias_table @}

## variance by experiment
### Plot

Errorbars from performing bootstrap sample upon replicas

{@ plot_fits_bootstrap_variance @}
### Table
{@ fits_bootstrap_variance_table @}

## $\Delta_{\chi^{2}}$
{@ plot_delta_chi2 @}
### By experiment
{@ delta_chi2_table @}

Dataset properties
------------------
{@closure fit_datasets_properties_table@}

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

$\chi^2$ by {@processed_metadata_group@}
----------------------
{@plot_fits_groups_data_chi2@}

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
