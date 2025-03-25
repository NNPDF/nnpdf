<div align="center">
  <img src="doc/sphinx/source/_static/LogoNNPDF.png" height=100>
</div>

[![NNPDF test suite](https://github.com/NNPDF/nnpdf/actions/workflows/all_tests_nnpdf.yml/badge.svg)](https://github.com/NNPDF/nnpdf/actions/workflows/all_tests_nnpdf.yml)
[![Docs](https://github.com/NNPDF/nnpdf/actions/workflows/upload_docs.yml/badge.svg)](https://github.com/NNPDF/nnpdf/actions/workflows/upload_docs.yml)
[![Commondata](https://github.com/NNPDF/nnpdf/actions/workflows/check_newcd.yml/badge.svg)](https://github.com/NNPDF/nnpdf/actions/workflows/check_newcd.yml)

[![DOI](https://zenodo.org/badge/118135201.svg)](https://zenodo.org/badge/latestdoi/118135201)

# NNPDF: An open-source machine learning framework for global analyses of parton distributions

[The NNPDF collaboration](http://nnpdf.science) determines the structure of the
proton using Machine Learning methods. This is the main repository of the
fitting and analysis frameworks. In particular it contains all the necessary
tools to [reproduce](https://docs.nnpdf.science/tutorials/reproduce.html) the
[NNPDF4.0 PDF determinations](https://arxiv.org/abs/2109.02653).

## Documentation

The documentation is available at <https://docs.nnpdf.science/>

## Install

See the [NNPDF installation guide](https://docs.nnpdf.science/get-started/installation.html)
for instructions on how to install and use the code, requirements and [dependencies](https://docs.nnpdf.science/get-started/installation.html#dependencies-and-requirements)
As a first step we recommend to follow one of the [tutorials](https://docs.nnpdf.science/tutorials/run-fit.html).

We follow a rolling development model where the tip of the master branch is
expected to be stable, tested and correct. For more information see our
[releases and compatibility policy](https://docs.nnpdf.science/releases.html).

## Cite

This code is described in the following [paper](https://inspirehep.net/literature?sort=mostrecent&size=25&page=1&q=find%20eprint%202109.02671):

```
@article{NNPDF:2021uiq,
    author = "Ball, Richard D. and others",
    collaboration = "NNPDF",
    title = "{An open-source machine learning framework for global analyses of parton distributions}",
    eprint = "2109.02671",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "Edinburgh 2021/13, Nikhef-2021-020, TIF-UNIMI-2021-12",
    doi = "10.1140/epjc/s10052-021-09747-9",
    journal = "Eur. Phys. J. C",
    volume = "81",
    number = "10",
    pages = "958",
    year = "2021"
}
```

If you use the code to produce new results in a scientific publication, please
follow the [Citation Policy](https://docs.nnpdf.science/get-started/cite.html),
particularly in regards to the papers relevant for QCD NNLO and EW NLO
calculations incorporated in the NNPDF dataset.

## Contribute

We welcome bug reports or feature requests sent to the [issue
tracker](https://github.com/NNPDF/nnpdf/issues). You may use the issue tracker
for help and questions as well.

If you would like contribute to the code, please follow the [Contribution
Guidelines](https://docs.nnpdf.science/contributing/index.html).



When developing locally you can test your changes with pytest, running from the root of the repository:

```
  pytest --mpl --pyargs n3fit validphys
```
