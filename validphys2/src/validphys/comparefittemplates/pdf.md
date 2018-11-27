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
{@with pdfs@}
### Effective exponents table for {@pdf@}
{@effective_exponents_table@}
{@endwith@}

{@with pdfs@}
### Next effective exponents table for {@pdf@}
{@next_effective_exponents_yaml@}
{@endwith@}