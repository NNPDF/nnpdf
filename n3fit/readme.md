## The NNPDF fitting code

This is the main directory of the new generation fitting code for the NNPDF collaboration.
For a more technical description of the code, please consult the [NNPDF documentation](https://docs.nnpdf.science/n3fit/methodology.html).
The first release of NNPDF fits using this code corresponds to the [NNPDF40 release](https://inspirehep.net/literature/1918284).

### How to install

It is recommended to utilize the officially released conda package which can be installed with:

```bash
conda install nnpdf -c https://packages.nnpdf.science/conda
```

If you are a developer please consult the [guidelines](https://docs.nnpdf.science/get-started/rules.html) before starting.

### How to run a fit

Runcards examples are provided in the `runcards` directory in this folder.
Specific runcard options are discussed [here](https://docs.nnpdf.science/n3fit/runcard_detailed.html).

The basic command is:

```bash
n3fit runcard.yml <replica_number>
```

For a full tutorial on how to run a fit using this code please see [here](https://docs.nnpdf.science/tutorials/run-fit.html)
Other relevant tutorials are the ones dedicated to [fit analysis](https://docs.nnpdf.science/tutorials/index.html#analysing-results) using [validphys 2](https://docs.nnpdf.science/vp/index.html).

### Reproducibility

Please note that this repository is in active development and only releases tagged as such are guaranteed to produce correct results.

In order to ensure reproducibility all official NNPDF releases are run with a fixed version of not only the `nnpdf` fitting code but also
of all external packages that could potentially influence the results.
A list of conda environments with the exact versions of the `nnpdf` fitting code as well as all external packages used for official releases
of the code can be found in the `environments` directory.

### References

The papers describing the implementation, most significant developments and release of the code are:

- _Towards a new generation of parton densities with deep learning models_
**S. Carrazza, J. Cruz-Martinez**, [Eur.Phys.J.C 79 (2019) 8, 676](https://doi.org/10.1140/epjc/s10052-019-7197-2)
- _NNPDF40 paper_
**NNPDF Collaboration** [Eur.Phys.J.C 82 (2022) 5, 428](https://doi.org/10.1140/epjc/s10052-022-10328-7)
