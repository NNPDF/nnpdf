%{@pdfcomptitle@}

Comparison between PDFs:

{@with pdfs@}
- **{@pdf@}** (LHAPDF ID: `{@pdf_id@}`)
{@endwith@}

{@with Qs@}

# PDF comparison at {@Q@} GeV

## All flavours

{@with pdfs@}

### {@pdf@}

{@pdgbasis plot_flavours@}

{@endwith@}

{@with basisspec@}

## {@basistitle@}

### Absolute

{@plot_pdfs@}

### Normalized

{@normalize_to_first plot_pdfs@}

### Distances


{@normalize_to_first plot_pdfdistances@}


{@endwith@}
{@endwith@}