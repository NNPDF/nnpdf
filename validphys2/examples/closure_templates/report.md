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
### Scatter Plot
{@ plot_fits_bootstrap_bias @}
### Table
{@ biases_table(show_total=True) @}

## $\phi$ by experiment
### Scatter Plot
{@ plot_phi_scatter_dataspecs @}
### Table
{@ fits_experiments_phi_table @}

## $\Delta_{\chi^{2}}$
{@ plot_delta_chi2 @}
### By experiment
{@ delta_chi2_table @}

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
