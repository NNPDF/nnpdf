```eval_rst
.. _pytools:
```

# Tools for developing with the Python programming language

This page summarizes auxiliary Python tools that we use commonly to develop the
project. Note that this page is meant to be a quick index. Consult the
documentation on each specific tool for details.


## Python editors

  - The [Sypder editor](https://www.spyder-ide.org/) is good for getting started
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
places in the code. A  few alternatives exists when that is not enough:

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
    can be used as an alternative. Compared to IPython it has the advantage that
    it allows to atomically step in the execution of the code, but the disadvantage
    that the interface is somewhat more complex and often surprising (hint: always
    [prefix interpreter commands with `!`](https://docs.python.org/3/library/pdb.html#pdbcommand-!)).

## Performance profiling

Sometimes a piece of code runs slower than expected. The reasons can often be
surprising. It is a good idea to measure where the problems actually are.

  - [`py-spy`](https://github.com/benfred/py-spy): Is a modern tool that
    provides good information and little overhead. Prefer it to the standard
    `cProfile`. The output is typically presented in the form of "Flamegraphs"
    that show the relative time spent on each piece of code.

## Documentation

  - We use the [Sphinx tool](https://www.sphinx-doc.org/) to document code
    projects. It an render and organize special purpose documentation files as
    well as read Python source files to automatically document interfaces.  It
    supports extensive customization and plugins. In particular because the
    default formatting for docstrings is somewhat unwieldy, it is recommendable
    to enable the `napoleaon` extension which allows for a more lenient
    [`numydoc`](https://numpydoc.readthedocs.io/en/latest/format.html) style.
    Similarly the default RST markup language can be overwhelming for simple
    documents. We enable the
    [recommonmark](https://recommonmark.readthedocs.io/en/latest/) extension to
    be able to compose files also in markdown format.

