%NNPDF Report for fit {@ current fit @}

Dataset properties
------------------
{@current datasets_properties_table@}

Distances
------------------
{@normalize::basespecs::pdfscalespecs plot_pdfdistances@}

PDF arc-lengths
---------------
{@basespecs plot_arc_lengths@}

PDF plots
---------
{@with pdfnormalize@}
## {@normtitle@}
{@with basespecs@}
### {@basistitle@}
{@with pdfscalespecs@}
#### {@xscaletitle@}
{@plot_pdfs@}
{@endwith@}
{@endwith@}
{@endwith@}

PDF replicas
------------
{@with basespecs@}
### {@basistitle@}
{@with pdfscalespecs@}
#### {@xscaletitle@}
{@plot_pdfreplicas@}
{@endwith@}
{@endwith@}

Effective preprocessing exponents
---------------------------------
{@plot_alphaEff@}
{@plot_betaEff@}

$\chi^2$ by experiment
----------------------
{@plot_fits_experiments_chi2@}

$\chi^2$ by dataset comparisons
-------------------------------
### Plot
{@plot_fits_datasets_chi2@}
### Table
{@fits_chi2_table(show_total=true)@}

$\phi$ by experiment
--------------------
{@with dataspecs@}
{@plot_phi_pdfs@}
{@endwith@}

{@with dataspecs@}
### {@fit@}
{@plot_phi@}
{@endwith@}

Experiment plots
---------------
{@with matched_datasets_from_dataspecs@}
[Detailed plots for dataset ' {@dataset_name@} ']({@dataset_report report@})
{@endwith@}
