% Comparing closure test {@current fit@} and {@reference fit@}.

Fit summary
-----------

We are comparing:

  - {@ current fit @} (`{@ current fit_id @}`): {@ current description @}
  - {@ reference fit @} (`{@ reference fit_id @}`): {@ reference description @}

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

PDF plots
---------
[Detailed plots]({@current::fitunderlyinglaw::pdf_report report@})