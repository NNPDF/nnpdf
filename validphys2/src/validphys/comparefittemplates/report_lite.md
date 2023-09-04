%NNPDF light report comparing {@ current fit_id @} and {@ reference fit_id @}

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

Training lengths
----------------
{@fits plot_training_length@}

Training-validation
-------------------
{@fits plot_training_validation@}

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
