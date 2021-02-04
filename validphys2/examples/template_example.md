NNPDF Report
============

{@ description  @}


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
{@ with pdfs @}
###Experiment results for {@pdf@}
{@with datanorm::experiments@}

#### {@experiment@}
{@experiment plot_fancy @}
{@ endwith @}
{@ endwith @}
