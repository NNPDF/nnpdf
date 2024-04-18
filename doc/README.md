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

It is required to install `recommonmark` to interpret markdown. To add the
dependencies to your environment, run

```
conda install sphinx recommonmark
```

### Adding to the Documentation

New documentation can be added in markdown (`.md` or `.txt` suffices) or
restructured text (`.rst` suffix) formats. To add a new section to the
documentation, create an appropriately named directory in the `sphinx/source/`
directory.  Inside the new directory, add all relevant documentation in the
markdown or restructured text formats. In addition to these files, create an
`index.rst` file containing:

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

### Useful Markdown and Restructured Text Tools

Various
[markdown](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet) and
[restructured text](http://docutils.sourceforge.net/docs/user/rst/quickref.html)
cheatsheets exist online.

In restructured text, a $\LaTeX$ block can be generated using

```
.. math::

   \frac{1}{2}
```

while inline maths is generated using

```
:math:`\frac{1}{2}`
```
with attention being brought to the backticks. Note: the markdown intepreter
being used here does not support inline maths, so if formula dense documentation
is being implemented, it is advised to use restructured text instead.

One can cross reference various parts of their markdown file using `anchors`,
which provide clickable pieces of text which transport the reader to a
particular part of the document.

To do this: add an anchor point in the text. This may look like the following:
```
Lorem ipsum dolor sit amet <a name="label"</a> consectetur adipiscing elit, sed do
```

we can then jump to `label` from an arbitrary point in the text by using
`[text](#label)`

As an example, clicking [this](#top) will take the reader to the top of the
page.

This was done by having the following lines of code:

```
For example, clicking [this](#top) will take the reader to the top of the page.
```
as well as
```
# NNPDF code and standards documentation <a name="top"></a>
```
at the top of this file.


In addition, one can link to other pages within the documentation by
`[text](<relative-path-to-md-or-rst-file>.<extension>)`.

One can define "lables" for RestructuredText, which can be referred to from
anywhere, like this:
```
    .. _my-reference-label:

    Section to cross-reference
    --------------------------

    This is the text of the section.

    It refers to the section itself, see :ref:`my-reference-label`.
```

Such labels can also be defined in Markdown by using `rst` syntax embedded in
code markers in markdown:


	```eval_rst
	.. _my-reference-label:
	```

Labels can be linked to from anywhere using  the syntax

```
[lint text](my-refence-label)
```
for Markdown and

```
:ref:`my-reference-label`.
```
for RestructuredText, as described in its
[documentation](https://www.sphinx-doc.org/en/master/usage/restructuredtext/roles.html?highlight=cross%20reference#role-ref).

## Installation using conda



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
