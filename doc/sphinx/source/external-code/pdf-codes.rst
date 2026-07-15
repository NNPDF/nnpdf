.. _lhapdf-code:

PDF set storage and interpolation
=================================

`LHAPDF <https://lhapdf.hepforge.org/>`_ is a C++ library that evaluates PDFs by interpolating the
discretised PDF 'grids' that PDF collaborations produce. It also gives its users access to proton
and nuclear PDF sets from a variety of PDF collaborations, including NNPDF, MMHT and CTEQ. A list
of all currently available PDF sets can be found on their
`website <https://lhapdf.hepforge.org/pdfsets.html>`_. Particle physics programmes that typically make
use of PDFs, such as Monte Carlo event generators, will usually be interfaced with LHAPDF, to allow
a user to easily specify the PDF set that they wish to use in their calculations. You can read more
about LHAPDF by reading the `paper <https://arxiv.org/abs/1412.7420>`_ that marked their latest
release.

PDF evolution
-------------

The evolution of PDFs is fully handled in the `EKO <https://github.com/nnpdf/eko>`_ code,
and the way QCD evolution is calculated is described in great detail `here <https://eko.readthedocs.io/>`_.

PDF compression
---------------

PDF compression seeks to maintain the statistical accuracy of a large sample of replicas
produced by a fit when using a PDF set with a smaller number of replicas (and thus fewer 
convolutions required to compute cross sections with PDF uncertainties). For example the 
main published PDFs are typically based on a 1000 replica fit, which can then be compressed to 
around a 100 replicas PDF set while maintaining good accuracy of most relevant statistical estimators.
This is done with the `pyCompressor <https://n3pdf.github.io/pycompressor/>`_ library,
a python compression code that extracts, from an initial PDF set of replicas,
the subset that most truthfully reproduces the underlying probability distribution of the prior. 
`pyCompressor <https://n3pdf.github.io/pycompressor/>`_ is an updated python version of
`compressor <https://github.com/scarrazza/compressor>`_, which was used in previous releases.

Other codes
~~~~~~~~~~~

`Hoppet <https://github.com/hoppet-code/hoppet>`_ ('Higher Order Perturbative Parton Evolution Toolkit') is an
alternative PDF evolution code which is capable of evolving unpolarised PDFs.
