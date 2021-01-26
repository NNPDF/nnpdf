%NNPDF Report for fit {@ current fit_id @}

Fit summary
------------

We are comparing:

  - {@ current fit @} (`{@ current fit_id @}`): {@ current description @}
  - {@ reference fit @} (`{@ reference fit_id @}`): {@ reference description @}


{@ summarise_fits @}

Fit code versions
-----------------

{@fits_fit_code_version@}

Theory Covariance Summary
-------------------------
{@summarise_theory_covmat_fits@}

Dataset properties
------------------
{@current fit_datasets_properties_table@}

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

$\chi^2$ by {@processed_metadata_group@}
----------------------
{@plot_fits_groups_data_chi2@}

$\chi^2$ by dataset comparisons
-------------------------------
### Plot
{@plot_fits_datasets_chi2@}
### Table
{@fits_chi2_table(show_total=true)@}

$\phi$ by {@processed_metadata_group@}
--------------------
{@plot_fits_groups_data_phi@}

Dataset plots
---------------
{@with matched_datasets_from_dataspecs@}
[Detailed plots for dataset ' {@dataset_name@} ']({@dataset_report report@})
{@endwith@}

Positivity
----------
{@with matched_positivity_from_dataspecs@}
{@plot_dataspecs_positivity@}
{@endwith@}

Dataset differences and cuts
----------------------------
{@print_dataset_differences@}
{@print_different_cuts@}
