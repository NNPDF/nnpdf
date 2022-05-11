.. _developer:

Developing `validphys`
======================

.. important::

	You can find some general tips on developing with python in
	:ref:`pytools`.

`validphys2` aims to be as simple to understand and extend as
possible. The code is based on self-contained Python functions with
a couple of magic decorators that make `reportengine` work as
expected. Based on that, there is a large and ever growing set of
tools to interact with NNPDF resources, which are spread across several
modules in the codebase.

Key modules
------------

Some of the most important modules are

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
Contains parsers for various files produced by the fitting code along with
tools to manipulate and display them.

- `validphys.checks`
Contains `reportengine`-style checks that are used in several
places.

These are used as a basis for doing everything else. For
implementing new functionality see :ref:`addvpplots`.

Unfortunately, the objective of making `validphys` easy means that the
complexity of getting things to just work is translated into
`reportengine`, which instead uses many advanced python features, and
results in a codebase that is not particularly simple.
