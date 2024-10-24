.. _pytools:

Tools for developing with the Python programming language
=========================================================

This page summarizes auxiliary Python tools that we commonly use to develop the
project. Note that this page is meant to be a quick index. Consult the
documentation on each specific tool for details.


Python editors
--------------

  - The `Spyder editor <https://www.spyder-ide.org/>`_ is good for getting started
    with scientific Python coding, because of various inspection and interactive
    features.
  - `vscode <https://code.visualstudio.com/>`_ is a more full featured editor.
  - In the long run, the most efficient approach is to learn a terminal based
    editor such as `vim <https://www.vim.org/>`_. Note that `vim` editing modes
    can be added as extensions to graphical editors such as :code:`vscode`.


Interactive development
-----------------------

Python code can be evaluated interactively, which can speed up the development.

  - `IPython shell <https://ipython.org/>`_: It is notably nicer to use than the
    standard interactive interpreter.
  - `Jupyter notebook <https://jupyter.org/>`_: Interactive development
    environment running on the browser. Useful for bigger experiments.

.. note::
    When developing :ref:`validphys <validphys>` related code interactively, be
    sure to read about the :ref:`API functionality <vpapi>`.

Testing
-------

  - `pytest <https://docs.pytest.org/en/latest/>`_: It is a framework for
    writing and running tests. It finds tests in the codebase (basically
    modules and functions that start with ``test``), enhances the ``assert``
    statement to provide rich error reporting and allows to structure
    dependencies between the tests (in a way similar to ``reportengine``).
    Tests are stored in the codebase and executed by pytest either manually or
    as a part of the continuous integration process.
  - `coverage.py <https://coverage.readthedocs.io/en/coverage-5.2.1/>`_ is a
    program that traces which lines of code have been executed when a given
    Python program (notably pytest) is running. The main use case is to verify
    that tests probe our code paths.


.. _pytoolsqa:

Code quality and reviewing
--------------------------

See also :ref:`reviewing pull requests <reviews>`. Note that these can typically be
integrated with your editor of choice.

  - The `pylint <https://www.pylint.org/>`_ tool allows for the catching of
    common problems in Python code. The top level
    ``.pylintrc`` `file <https://github.com/NNPDF/nnpdf/blob/master/.pylintrc>`_
    comes with a useful and not overly noisy configuration.
  - New Python code should come formatted with
    ``black`` `tool <https://github.com/psf/black>`_ with `our default
    configuration <https://github.com/NNPDF/nnpdf/blob/master/pyproject.toml>`_
  - The ``isort`` `library <https://pycqa.github.io/isort/>`_ sorts imports
    alphabetically, and automatically separated into sections and by type.
  - `pre-commit <https://pre-commit.com/>`_ is a tool that, can automatically
    check for stylistic problems in code such as trailing whitespaces or
    forgotten debug statements. Our configuration can be found in
    `.pre-commit-configuration.yaml <https://github.com/NNPDF/nnpdf/blob/master/.pre-commit-configuration.yaml>`_
    and also ensures that ``black`` and ``isort`` are run.


Debugging
---------

Usually the most efficient way to debug a piece of Python code, such as a
``validphys`` action is to insert ``print`` statements to check the state at various
places in the code. A few alternatives exists when that is not enough:

  - `IPython embed <https://ipython.readthedocs.io/en/stable/api/generated/IPython.terminal.embed.html>`_:
    The `IPython <https://ipython.org/>`_ shell can be easily dropped at any
    arbitrary point in the code. Write

    .. code:: python

      import IPython

      IPython.embed()

    at the location of the code you want to debug. You will then be able to
    query (and manipulate) the state of the code using a rich shell.

  - PDB: The standard `Python debugger <https://docs.python.org/3/library/pdb.html>`_
    can be used as an alternative. Compared to ``IPython`` it has the advantage that
    it allows to automatically step in the execution of the code, but the disadvantage
    that the interface is somewhat more complex and often surprising (hint: always
    `prefix interpreter commands with <https://docs.python.org/3/library/pdb.html#pdbcommand-!>`_ ``!``.

Performance profiling
---------------------

Sometimes a piece of code runs slower than expected. The reasons can often be
surprising. It is a good idea to measure where the problems actually are.

  - `py-spy <https://github.com/benfred/py-spy>`_: A performance measuring
    program (*profiler*) that provides good information and little overhead.
    Prefer it to the standard ``cProfile``. The output is typically presented in
    the form of "Flamegraphs" that show the relative time spent on each piece of
    code.

Documentation
-------------

  - We use the `Sphinx tool <https://www.sphinx-doc.org/>`_ to document code
    projects. It can render and organize special purpose documentation files as
    well as read Python source files to automatically document interfaces.  It
    supports extensive customization and plugins. In particular because the
    default formatting for docstrings is somewhat unwieldy, it is recommended
    to enable the ``napoleon`` extension which allows for a more lenient
    `numpydoc <https://numpydoc.readthedocs.io/en/latest/format.html>`_ style.
    Similarly the default RST markup language can be overwhelming for simple
    documents.

Python static checks and code style

We use `Pylint <https://www.pylint.org/>`_ to provide static checking e.g.
finding basic errors that a compiler would catch in compiled languages. An example
is using an unknown variable name. Pylint also provides basic guidelines on the
structure of the code (e.g. avoid functions that are to complicated). Because
Pylint is way too pedantic by default, we limit the checks to only those
considered useful. The ``.pylintrc`` file at the top level configures Pylint to
only mind those checks. Most Python IDEs and editors have some kind of support
for Pylint. It is strongly recommended to configure the editor to show the
problematic pieces of code proactively.

New code should use the `Black <https://black.readthedocs.io/en/stable/>`_ tool to
format the code. This tool should not be used to aggressively reformat existing
files.


Matplotlib Image Comparison Tests
---------------------------------

It is possible to create tests which perform an image comparison between a
generated plot and a pre-existing baseline plot. Clearly this allows one to check
consistency in figure generation.

Before beginning you will need to ensure that you have the tests dependencies,
which can be checked in :code:`nnpdf/conda-recipe/meta.yml`.

The next step is to write the test function. It is highly recommended to use the
:ref:`validphys API <api>` for this, both to simplify the code and to make it agnostic to the
structure of backend providers - provided that they produce the same results. See
for example a function which tests the ``plot_pdfs`` provider:

.. code:: python

  @pytest.mark.mpl_image_compare
  def test_plotpdfs():
      pdfs = ["NNPDF31_nnlo_as_0118"]
      Q = 10
      flavours = ["g"]
      # plot_pdfs returns a generator with (figure, name_hint)
      return next(API.plot_pdfs(pdfs=pdfs, Q=Q, flavours=flavours))[0]

We see that the function needs to return a valid matplotlib figure, and should
be decorated with :code:`@pytest.mark.mpl_image_compare`.

Now the baseline figure needs to be generated, this can be achieved by running

.. code:: bash

  pytest -k <name of file containing test function> --mpl-generate-path=baseline


which will generate a PNG of the figure in the :code:`src/validphys/tests/baseline`
directory. It is recommended to put all baseline plots in this directory so that
they are automatically installed, and so will be in the correct location when
the :ref:`CI <CI>` runs the test suite.

Now that the baseline figure exists you can check that your test works:

.. code:: bash

  pytest -k <name of file containing test function> --mpl

Also you can check that the test has been added to the full test suite:

.. code:: bash

  pytest --pyargs --mpl validphys

Just note that if you do not put the :code:`--mpl` flag then the test will just check
that the function runs without error, and won't check that the output matches to
baseline.
