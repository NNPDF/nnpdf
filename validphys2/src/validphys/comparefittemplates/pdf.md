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
{@plot_alphaEff@}
### beta exponent
{@plot_betaEff@}

## Effective preprocessing exponents Table
{@with fits@}
### Next effective exponents table for {@fit@}
{@effective_exponents_table@}
{@endwith@}
