.. _developer: 

Developing `validphys`
======================

`validphys2` aims to be as simple to understand and extend as
possible. The code is based on self contained Python functions with
a couple of magic decorators that make `reportengine` work as
expected. Based on that, there is a large and ever growing set of
tools to interact with NNPDF resources, that are spread across several
modules in the codebase. Some of the most important are:

- `validphys.core`
Core data structures that represent objects such as PDFs and data
sets. Several of them map to `libnnpdf` objects. In that case they
have a `.load()` method that produces the corresponding `C++`
object.

- `validphys.loader`
Tools to obtain NNPDF resources locally or remotely. See :ref:`upload`
and :ref:`download`.

- `validphys.config`
Defines how resources are to be parsed from the configuration
file. This is largely using `validphys.loader`.

- `validphys.results`
Implements tools to store and manipulate results from data and
theory predictions.

- `validphys.gridvalues`, `validphys.bases`, `validphys.pdfgrids`
These contain tools to evaluate PDFs over grids of points.
`validphys.gridvalues` contains low level functionality that uses
`libnnpdf`, `validphys.pdfbases` contain several different bases
over PDF flavour space and functionality to manipulate them, and
`validphys.pdfgrids` contains high level providers suitable for
using for plotting and as an input to other computations.

- `validphys.plotoptions`
Tools for interpreting the dataset PLOTTING files, including the
transformations on the kinematics and the data.

- `validphys.fitdata`
Contains parsers for various files produced by `nnfit` along with
tools to manipulate and display them.

- `validphys.checks`
Contains `reportengine`-style checks that are used in several
places. 


These are used as a basis for doing everything else. For 
implementing new functionality see :ref:`tut_newaction`.

Unfortunately the objective of making `validphys` easy means that the
complexity of getting things to just work is translated into
`reportengine`, which instead uses many advanced python features, and
results in a codebase that is not particularly simple.

Python static checks and code style
-----------------------------------

We use [Pylint](https://www.pylint.org/) to provide static checking (e.g.
finding basic errors that a compiler would catch in compiled languages) such as
uses of unknown variable names, as well as to provide basic guidelines on the
structure of the code (e.g. avoid functions that are too complicated). Because
Pylint is way too pendantic by default, we limit the checks to only those
considered useful. The `.pylintrc` file at the top level configures Pylint to
only mind those checks. Most Python IDEs and editors have some kind of support
for pylint. It is strongly recommended to configure the editor to show the
problematic pieces of code proactively.

New code should use the [Black](https://black.readthedocs.io/en/stable/) tool to
format the code. This tool should not be used to aggressively reformat existing
files.


Example pull request
--------------------

You may find instructive to go though this pull request that
implements arc-length computation:

<https://github.com/NNPDF/validphys2/pull/64>

It demonstrates how to leverage existing functionality to perform new
computations and then present those as plots and tables.


Matplotlib Image Comparison Tests
---------------------------------

It is possible to create tests which perform an image comparison between a
generated plot and a preexisting baseline plot. Clearly this allows one to check
consistency in figure generation.

Before beginning you will need to ensure that you have the tests dependencies,
which can be checked in `nnpdf/conda-recipe/meta.yml`.

The next step is to write the test function. It is highly recommended to use the
validphys API for this, both to simplify the code and to make it agnostic to the
structure of backend providers - provided that they produce the same results. See
for example a function which tests the `plot_pdfs` provider:

```python
@pytest.mark.mpl_image_compare
def test_plotpdfs():
    pdfs = ['NNPDF31_nnlo_as_0118']
    Q = 10
    flavours = ['g']
    #plot_pdfs returns a generator with (figure, name_hint)
    return next(API.plot_pdfs(pdfs=pdfs, Q=Q, flavours=flavours))[0]
```

we see that the function needs to return a valid matplotlib figure, and should
be decorated with `@pytest.mark.mpl_image_compare`.

Now the baseline figure needs to be generated, this can be achieved by running

```
pytest -k <name of file containing test function> --mpl-generate-path=baseline
```

which will generated a PNG of the figure in the `src/validphys/tests/baseline`
directory. It is recommended to put all baseline plots in this directory so that
they are automatically installed, and so will be in the correct location when
the CI runs the test suite.

Now that the baseline figure exists you can check that your test works:

```
pytest -k <name of file containing test function> --mpl
```

Also you can check that the test has been added to the full test suite:

```
pytest --pyargs --mpl validphys
```

just note that if you do not put the `--mpl` flag then the test will just check
that the function runs without error, and won't check that the output matches to
baseline.

