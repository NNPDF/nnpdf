[![DOI](https://zenodo.org/badge/42721933.svg)](https://zenodo.org/badge/latestdoi/42721933)
[![Build Status](https://travis-ci.org/NNPDF/reportengine.svg?branch=master)](https://travis-ci.org/NNPDF/reportengine)

Reportengine
============

Reportengine is a framework to develop scientific applications. It is
focused on supporting declarative input (YAML), enforcing
initialization time ("*compile time*") constraints, and enabling easy
iteration within the declarative input.

It includes support for figures (matplotlib), tables (pandas) and HTML
reports (pandoc-markdown). It also tries to make the command line
applications look like from the 90s as opposed to from the 70s.

The documentation of the NNPDF specific implementation can be found
here:

http://pcteserver.mi.infn.it/~nnpdf/validphys-docs/guide.html

a more reportengine-specific documentation will be produced *soon*.

An example application can be found in the `example` directory.


Install
-------

For linux, you can install a precompiled package by running

````
conda install reportengine -c https://zigzah.com/static/conda-pkgs

````

Alternatively, you can satisfy all the dependencies automatically by
running:

````
conda build conda-recipe
````

and then installing the resulting package.


Development
-----------

Install in development mode:

````
pip install -e .
````

Running the tests
-----------------

Easiest way is:

````
py.test
````
