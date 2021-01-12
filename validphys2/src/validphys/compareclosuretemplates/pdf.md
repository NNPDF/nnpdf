%NNPDF Report for fits {@ current fit @} and {@ reference fit @}
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
{@current::basisfromfit plot_alpha_eff@}
### beta exponent
{@current::basisfromfit plot_beta_eff@}
