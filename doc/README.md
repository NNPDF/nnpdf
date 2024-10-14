# NNPDF code and standards documentation <a name="top"></a>

Here we store the [documentation](https://docs.nnpdf.science/) (user / developer
guides)

## Sphinx Documentation

### Generating the Documentation

The NNPDF documentation is produced by the
[sphinx](http://www.sphinx-doc.org/en/master/) resource. To generate the sphinx
documentation, navigate to the `sphinx/` directory and execute the command `make
html`. This produces the documentation in the `build/index/` directory. The
`index.html` can be viewed with any appropriate browser.

### Adding to the Documentation

New documentation can be added in restructured text (`.rst`) format. To
add a new section to the documentation, create an appropriately named directory
in the `sphinx/source/` directory.  Inside the new directory, add all relevant
documentation in the restructured text formats. In addition to these
files, create an `index.rst` file containing:

```
Chapter Name
============

.. toctree::
   :maxdepth: 1

   ./file1.md
   ./file2.rst
```
ensuring that the number of `=` signs is the same as the number of characters in
`Chapter Name`.

The next step is to reference the newly made `index.rst` in the main
`sphinx/source/index.rst` file:

```
.. NNPDF documentation master file, created by
   sphinx-quickstart on Mon Oct 29 10:53:50 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

NNPDF documentation
===================

.. toctree::
   :maxdepth: 2

   get-started/index
   theory/index
   vp/index
   code/index
   tutorials/index
   QA/index
   <NEW CHAPTER>/index

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
```

### Adding indices for modules

Sphinx has the capability of automatically documenting any python package. It
produces these under the `index` and `module index` sections. The functions and
modules are documented using their corresponding docstrings.

To add a new module to document, add a new line in `sphinx/Makefile` under:

```
%: Makefile
	@if test $@ != "clean"; then
            sphinx-apidoc -o ./source/modules/validphys ../../validphys2/src/validphys/ ; \
            sphinx-apidoc -o ./source/modules/<MODULE-NAME> <PATH-TO-MODULE>  ;\
	fi

```
