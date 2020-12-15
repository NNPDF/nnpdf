Covariance matrices
-------------------
{@with default_theory@}
   {@plot_normexpcovmat_heatmap@}
   {@plot_normthcovmat_heatmap_custom@}
{@endwith@}

Correlation matrices
--------------------
{@with default_theory@}
   {@plot_expcorrmat_heatmap@}
   {@plot_thcorrmat_heatmap_custom@}
   {@plot_expplusthcorrmat_heatmap_custom@}
{@endwith@}

Diagonal elements of covariance matrices
----------------------------------------
{@with default_theory@}
   {@plot_diag_cov_comparison@}
{@endwith@}

Experimental $\chi^2$
---------------------
{@with default_theory@}
   {@total_experiments_chi2@}

Total (exp. + th.) $\chi^2$
---------------------------
   {@chi2_impact_custom@}

Experimental $\chi^2$ by dataset
--------------------------------
   {@experiments_chi2_table@}

Total (exp. + th.) $\chi^2$ by dataset
--------------------------------------
   {@experiments_chi2_table_theory@}

$\chi^2$ including only diagonal theory elements
------------------------------------------------
   {@chi2_diag_only@}

Impact of theory covariance matrix on $\chi^2$s
-----------------------------------------------
   {@plot_datasets_chi2_theory@}
{@endwith@}

Scale variations as a function of the kinematics
------------------------------------------------
{@with matched_datasets_from_dataspecs@}
   [Plots for {@dataset_name@}]({@dataset_report report@})
{@endwith@}
