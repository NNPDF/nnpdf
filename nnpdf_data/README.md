NNPDF data package
==================

The central data repository for the [NNPDF collaboration](http://nnpdf.science). Provides experimental datasets, theory cards, and the ``nnpdf-data`` CLI helper.

## Install

```bash
pip install nnpdf-data           # basic
pip install nnpdf-data[filter]   # with filter dependencies to regenerate commondata
```

## Quick start

It can be used programatically to load information about a given dataset:

```python
from nnpdf_data import load_dataset_metadata

ds = load_dataset_metadata("LHCB_Z0_7TEV_DIELECTRON_Y")
print(ds.name, ds.ndata, ds.kinlabels)
```

Or used through the used interface
```bash
nnpdf-data list "ATLAS_*"                # list datasets
nnpdf-data latex my_runcard.yaml         # generate LaTeX tables + BibTeX
```

## Cite

If the ``nnpdf-data`` package was useful for your research, please
consider citing it.

```bibtex
@article{NNPDF:2021uiq,
    author = "Ball, Richard D. and others",
    collaboration = "NNPDF",
    title = "{An open-source machine learning framework for global analyses of parton distributions}",
    eprint = "2109.02671",
    doi = "10.1140/epjc/s10052-021-09747-9",
    journal = "Eur. Phys. J. C", volume = "81", pages = "958",
    year = "2021"
}
```

Like the rest of the NNPDF framework, ``nnpdf-data`` is licensed under GPLv3.
