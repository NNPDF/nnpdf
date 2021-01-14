%NNPDF Report for fit {@ current fit @}

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

## Effective preprocessing exponents Plots
### alpha exponent
{@current::basisfromfit plot_alpha_eff@}
### beta exponent
{@current::basisfromfit plot_beta_eff@}

## Effective preprocessing exponents Table
{@with fits@}
### Next effective exponents table for {@fit@}
{@effective_exponents_table@}
{@endwith@}
