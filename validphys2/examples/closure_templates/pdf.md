%NNPDF Report for fit {@ closure fit @}
# PDF plots

## PDF comparison
{@with pdfnormalize@}
### {@normtitle@}
{@with basespecs@}
#### {@basistitle@}
{@with pdfscalespecs@}
##### {@xscaletitle@}
{@plot_pdfs@}
{@endwith@}
{@endwith@}
{@endwith@}

## PDF replicas
{@with basespecs@}
#### {@basistitle@}
{@with pdfscalespecs@}
##### {@xscaletitle@}
{@plot_pdfreplicas@}
{@endwith@}
{@endwith@}

## Effective preprocessing exponents
### alpha exponent
{@closure::basisfromfit plot_alpha_eff@}
### beta exponent
{@closure::basisfromfit plot_beta_eff@}
