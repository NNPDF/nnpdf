%NNPDF report comparing {@ current fit_id @} and {@ reference fit_id @}

Summary
-------

We are comparing:

  - {@ current fit @} (`{@ current fit_id @}`): {@ current description @}
  - {@ reference fit @} (`{@ reference fit_id @}`): {@ reference description @}


{@ summarise_fits @}

Datasets excluded from fit
--------------------------
{@with matched_excluded_datasets_by_name@}
[Plots for {@dataset_name@}]({@plot_fancy@})
{@endwith@}

Code versions
-------------
{@fits_version_table@}
