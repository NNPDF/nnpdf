```eval_rst
.. _pytools:
```

# Tools for developing with the Python programming language

This page summarizes auxiliary Python tools that we commonly use to develop the
project. Note that this page is meant to be a quick index. Consult the
documentation on each specific tool for details.


## Python editors

  - The [Spyder editor](https://www.spyder-ide.org/) is good for getting started
    with scientific Python coding, because of various inspection and interactive
    features.
  - [`vscode`](https://code.visualstudio.com/) is a more full featured editor.
  - In the long run, the most efficient approach is to learn a terminal based
    editor such as [`vim`](https://www.vim.org/). Note that `vim` editing modes can be
    added as extensions to graphical editors such as `vscode`.


## Interactive development

Python code can be evaluated interactively, which can speed up the development.

  - [IPython shell](https://ipython.org/): It is notably nicer to use than the
    standard interactive interpreter.
  - [Jupyter notebook](https://jupyter.org/): Interactive development
    environment running on the browser. Useful for bigger experiments.

```eval_rst
.. note::
    When developing :ref:`validphys <validphys>` related code interactively, be
    sure to read about the :ref:`API functionality <vpapi>`.
```

## Testing

  - [`pytest`](https://docs.pytest.org/en/latest/): It is a framework for
    writing and running tests. It finds tests in the codebase (basically
    modules and functions that start with `test`), enhances the `assert`
    statement to provide rich error reporting and allows to structure
    dependencies between the tests (in a way similar to `reportengine`).
    Tests are stored in the codebase and executed by pytest either manually or
    as a part of the continuous integration process.
  - [`coverage.py`](https://coverage.readthedocs.io/en/coverage-5.2.1/) is a
    program that traces which lines of code have been executed when a given
    Python program (notably pytest) is running. The main use case is to verify
    that tests probe our code paths.



```eval_rst
.. _pytoolsqa:
```
## Code quality and reviewing

See also [*Reviewing pull requests*](reviews). Note that these can typically be
integrated with your editor of choice.

  - The [`pylint`](https://www.pylint.org/) tool allows for the catching of
	common problems in Python code. The top level
	[`.pylintrc` file](https://github.com/NNPDF/nnpdf/blob/master/.pylintrc)
	comes with a useful and not overly noisy configuration.
  - The [`black` code formatter](https://github.com/psf/black) runs almost
    without configuration and produces typically good results. It is good to run
    it by default, to avoid spending time on formatting (or arguing about it).

## Debugging

Usually the most efficient way to debug a piece of Python code, such as a
`validphys` action is to insert `print` statements to check the state at various
places in the code. A few alternatives exists when that is not enough:

  - [IPython embed](https://ipython.readthedocs.io/en/stable/api/generated/IPython.terminal.embed.html):
    The [IPython](https://ipython.org/) shell can be easily dropped at any
    arbitrary point in the code. Write
    ```python
    import IPython
    IPython.embed()
    ```
    at the location of the code you want to debug. You will then be able to
    query (and manipulate) the state of the code using a rich shell.

  - PDB: The standard [Python debugger](https://docs.python.org/3/library/pdb.html)
    can be used as an alternative. Compared to `IPython` it has the advantage that
    it allows to automatically step in the execution of the code, but the disadvantage
    that the interface is somewhat more complex and often surprising (hint: always
    [prefix interpreter commands with `!`](https://docs.python.org/3/library/pdb.html#pdbcommand-!)).

## Performance profiling

Sometimes a piece of code runs slower than expected. The reasons can often be
surprising. It is a good idea to measure where the problems actually are.

  - [`py-spy`](https://github.com/benfred/py-spy): A performance measuring
    program (*profiler*) that provides good information and little overhead.
    Prefer it to the standard `cProfile`. The output is typically presented in
    the form of "Flamegraphs" that show the relative time spent on each piece of
    code.

## Documentation

  - We use the [Sphinx tool](https://www.sphinx-doc.org/) to document code
    projects. It can render and organize special purpose documentation files as
    well as read Python source files to automatically document interfaces.  It
    supports extensive customization and plugins. In particular because the
    default formatting for docstrings is somewhat unwieldy, it is recommended
    to enable the `napoleon` extension which allows for a more lenient
    [`numpydoc`](https://numpydoc.readthedocs.io/en/latest/format.html) style.
    Similarly the default RST markup language can be overwhelming for simple
    documents. We enable the
    [recommonmark](https://recommonmark.readthedocs.io/en/latest/) extension to
    be able to compose files also in markdown format.

## Python static checks and code style

We use [Pylint](https://www.pylint.org/) to provide static checking e.g.
finding basic errors that a compiler would catch in compiled languages. An example
is using an unknown variable name. Pylint also provides basic guidelines on the
structure of the code (e.g. avoid functions that are to complicated). Because
Pylint is way too pedantic by default, we limit the checks to only those
considered useful. The `.pylintrc` file at the top level configures Pylint to
only mind those checks. Most Python IDEs and editors have some kind of support
for Pylint. It is strongly recommended to configure the editor to show the
problematic pieces of code proactively.

New code should use the [Black](https://black.readthedocs.io/en/stable/>) tool to
format the code. This tool should not be used to aggressively reformat existing
files.


## Matplotlib Image Comparison Tests

It is possible to create tests which perform an image comparison between a
generated plot and a pre-existing baseline plot. Clearly this allows one to check
consistency in figure generation.

Before beginning you will need to ensure that you have the tests dependencies,
which can be checked in `nnpdf/conda-recipe/meta.yml`.

The next step is to write the test function. It is highly recommended to use the
[validphys API](../vp/api.md) for this, both to simplify the code and to make it agnostic to the
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

We see that the function needs to return a valid matplotlib figure, and should
be decorated with `@pytest.mark.mpl_image_compare`.

Now the baseline figure needs to be generated, this can be achieved by running

```
pytest -k <name of file containing test function> --mpl-generate-path=baseline
```

which will generate a PNG of the figure in the `src/validphys/tests/baseline`
directory. It is recommended to put all baseline plots in this directory so that
they are automatically installed, and so will be in the correct location when
the [CI](../ci/index.md)  runs the test suite.

Now that the baseline figure exists you can check that your test works:

```
pytest -k <name of file containing test function> --mpl
```

Also you can check that the test has been added to the full test suite:

```
pytest --pyargs --mpl validphys
```

Just note that if you do not put the `--mpl` flag then the test will just check
that the function runs without error, and won't check that the output matches to
baseline.
