%NNPDF Report for fit {@ fit @}

PDF plots
---------

{@ plot_pdfs @}

**Normalized**

{@normalize plot_pdfs  @}

Train-valid split
------------------

{@ plot_training_validation @}

$\chi^2$
-------
{@ with pdfs  @}

### {@ pdf @}

{@ experiments_chi2_table @}

{@ endwith@}

Experiment plots
---------------
{@ with experiments @}
###{@ experiment @}
{@with experiment@}
[Detailed plots for dataset ' {@dataset@} ']({@dataset_report report @})
{@ endwith @}
{@ endwith @}

Theory description
------------------

{@ theory_description @}
