% Data theory comparison for {@dataset@}

Predictions using theory {@ theoryid @}

Absolute
---------

{@ plot_fancy @}

Normalized
----------

{@ datanorm plot_fancy@}

$\chi^2$
----

###Tables

{@ with pdfs @}
#### {@pdf@}
{@dataset_chi2_table @}
{@ endwith @}

### Replica distributions

{@pdfs plot_chi2dist @}


Data-PDF correlations
---------------------

{pdfs plot_smpdf}
