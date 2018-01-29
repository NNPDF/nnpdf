%{@pdfcomptitle@}

Here we compare the PDFs:

{@with pdfs@}
- **{@pdf@}** (LHAPDF ID: `{@pdf_id@}`)
{@endwith@}

We define the PDF *distance* between two PDF sets $a$ and $b$ as:

$$
D_{ab}(x,Q) =
10
\frac{\left|m_{a}(x,Q) - m_{b}(x,Q)\right|}
{\sqrt{s_a^2(x,Q)+s_b^2(x,Q)}}
$$

Where $m_a(x,Q)$ is the mean value at (x,Q) for the PDF set $a$, and
$s^2(x,Q)$ is the variance at (x,Q).



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

# Luminosities at $\sqrt{s}={@sqrts@} GeV$

## Luminosity ratios

We plot ratios of:

$$
L(M_{X},s)=\sum_{ij}^{\textrm{channel}}\frac{1}{s}
\int_{\tau}^{1}\frac{dx}{x}
f_{i}(x,M_{X})
f_{j}(\frac{\tau}{x},M_{X})
$$


Where $i$ and $j$ are summed as follows:

 - For $gg$, $i = j = g$.

 - For $qq$ i and j run over all possible quark and antquark
   combinations.

 - For $qg$, i and j are all possible combinations of a quark ant
   antiquark and a gluon.

 - For $q\bar{q}$, $i$ and $j$ correspond to a quark and antiquark of the
   same flavour.


{@with lumi_channels@}

### {@lumi_channel@}


{@normalize_to_first plot_lumi1d@}

{@endwith@}

## Luminosity uncertainties

We plot the percentage uncertainty of

$$
\tilde{L}(M_{X},y,s)=
\sum_{ij}^{\textrm{channel}}\frac{1}{s}
f_{i}\left(\frac{M_{x}e^{y}}{\sqrt{x}},M_{x}\right)
f_{j}\left(\frac{M_{x}e^{-y}}{\sqrt{x}},M_{x}\right)
$$

in the allowed kinematic range.

{@with lumi_channels@}

### {@lumi_channel@}

{@with pdfs@}

#### {@pdf@}

{@plot_lumi2d_uncertainty@}

{@endwith@}
{@endwith@}


