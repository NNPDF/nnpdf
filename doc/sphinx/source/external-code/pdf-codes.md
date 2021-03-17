```eval_rst
.. _lhapdf:
```
# PDF set storage and interpolation

*Author: Cameron Voisey, 13/10/2019*


[LHAPDF](https://lhapdf.hepforge.org/) is a C++ library that evaluates PDFs by interpolating the
discretised PDF 'grids' that PDF collaborations produce. It also gives its users access to proton
and nuclear PDF sets from a variety of PDF collaborations, including NNPDF, MMHT and CTEQ. A list
of all currently available PDF sets can be found on their
[website](https://lhapdf.hepforge.org/pdfsets.html). Particle physics programmes that typically make
use of PDFs, such as Monte Carlo event generators, will usually be interfaced with LHAPDF, to allow
a user to easily specify the PDF set that they wish to use in their calculations. You can read more
about LHAPDF by reading the [paper](https://arxiv.org/abs/1412.7420) that marked their latest
release.

## PDF evolution

[APFEL](https://apfel.hepforge.org/) ('A PDF Evolution Library') is the PDF evolution code currently
used by the NNPDF Collaboration. In addition to its PDF evolution capabilities, it also produces
predictions of deep-inelastic scattering structure functions. In recent years it has been developed
alongside NNPDF, and so it therefore contains the features and settings required in an NNPDF fit.
That is, it includes quark masses in the MSbar scheme, the various FONLL heavy quark schemes, scale
variations up to NLO, etc. Note that at the time of writing, a more streamlined code is being
written to replace APFEL, which is currently dubbed EKO ('Evolution Kernel Operator'). To find more
general information about PDF evolution and the DGLAP equations, you can go to the [Theory
section](dglap.md).

### Other codes

[Hoppet](https://hoppet.hepforge.org/) ('Higher Order Perturbative Parton Evolution Toolkit') is an
alternative PDF evolution code which is capable of evolving unpolarised PDFs to NNLO and linearly
polarised PDFs to NLO. The unpolarised evolution includes heavy-quark thresholds in the MSbar
scheme.
