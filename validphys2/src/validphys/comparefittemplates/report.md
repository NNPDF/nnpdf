%NNPDF report comparing {@ current fit_id @} and {@ reference fit_id @}

Summary
-------

We are comparing:

  - {@ current fit @} (`{@ current fit_id @}`): {@ current description @}
  - {@ reference fit @} (`{@ reference fit_id @}`): {@ reference description @}


{@ summarise_fits @}

Theory covariance summary
-------------------------
{@summarise_theory_covmat_fits@}

Dataset properties
------------------
{@current fit_datasets_properties_table@}

Distances
---------
{@with Scales@}
### {@Scaletitle@}
{@with Normalize::Basespecs::PDFscalespecs::Distspecs@}
#### {@Basistitle@}, {@Xscaletitle@}
{@plot_pdfdistances@}
{@plot_pdfvardistances@}
{@endwith@}
{@endwith@}

PDF arc-lengths
---------------
{@Basespecs plot_arc_lengths@}

Sum rules
---------
{@with pdfs@}
### {@pdf@}

#### Known sum rules

{@sum_rules_table@}

#### Unknown sum rules

{@unknown_sum_rules_table@}

{@endwith@}

PDF plots
---------
{@with Scales@}
[Plots at {@Scaletitle@}]({@pdf_report report@})
{@endwith@}

Luminosities
------------
{@with Energies@}
[Plots at {@Energytitle@}]({@lumi_report report@})
{@endwith@}

Effective exponents
-------------------
[Detailed information]({@exponents_report report@})

Training lengths
----------------
{@fits plot_training_length@}

Training-validation
-------------------
{@fits plot_training_validation@}

{@with DataGroups@}
$\chi^2$ by {@processed_metadata_group@}
----------------------------------------
{@plot_fits_groups_data_chi2@}
{@endwith@}


$\chi^2$ by dataset
-------------------
### Plot
{@plot_fits_datasets_chi2@}
### Table
{@ProcessGroup fits_chi2_table(show_total=true)@}


{@with DataGroups@}
$\phi$ by {@processed_metadata_group@}
--------------------------------------
{@plot_fits_groups_data_phi@}
{@endwith@}

Dataset plots
-------------
{@with matched_datasets_from_dataspecs@}
[Plots for {@dataset_name@}]({@dataset_report report@})
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

Code versions
-------------
{@fits_version_table@}
