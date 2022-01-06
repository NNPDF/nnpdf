.. _add_docs:

Adding to the Documentation
===========================

The NNPDF documentation is produced by the
`sphinx <http://www.sphinx-doc.org/en/master/>`__ resource. To generate
the sphinx documentation, navigate to the ``nnpdf/doc/sphinx/``
directory and execute the command ``make html``, ensuring one is inside
the appropriate ``nnpdf`` conda environment. This produces the
documentation in the ``build/index/`` directory. The ``index.html`` can
be viewed with any appropriate browser.

New documentation can be added in markdown, naming the source files with
the ``.md`` suffix, or restructured text, with the ``.rst`` suffix
formats.


.. note::
  The ``md`` format is now deprecated and only supported for legacy reasons.
  The reStructured Text format natively supports equation displaying as well as
  directives such as this note and is thus the preferred format for `NNPDF`
  documentation. Despite this, it is possible to evaluate inline ``rst`` in a
  ``md`` file using the ``eval_rst`` command in legacy files written in
  markdown.

To add a new section to the documentation, create an appropriately named
directory in the ``sphinx/source/`` directory. Inside the new directory,
add all relevant documentation in the markdown or restructured text
formats. In addition to these files, create an ``index.rst`` file
containing:

::

   Chapter Name
   ============

   .. toctree::
      :maxdepth: 1

      ./file1.md
      ./file2.rst

ensuring that the number of ``=`` signs is the same as the number of
characters in ``Chapter Name``.

The next step is to reference the newly made ``index.rst`` in the main
``sphinx/source/index.rst`` file:

::

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

Useful Markdown and Restructured Text Tools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Various
`markdown <https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet>`__
and `restructured
text <http://docutils.sourceforge.net/docs/user/rst/quickref.html>`__
cheatsheets exist online.

In restructured text, a :math:`\LaTeX` block can be generated using

::

   .. math::

      \frac{1}{2}

while inline maths is generated using

::

   :math:`\frac{1}{2}`

with attention being brought to the backticks. Note: the markdown
interpreter being used here does not support inline maths, so if formula
dense documentation is being implemented, it is advised to use
restructured text instead.

One can cross reference various parts of their markdown file using
``anchors``, which provide clickable pieces of text which transport the
reader to a particular part of the document.

To do this: add an anchor point in the text. This may look like the
following:

::

   Lorem ipsum dolor sit amet <a name="label"</a> consectetur adipiscing elit, sed do 

we can then jump to ``label`` from an arbitrary point in the text by
using ``[text](#label)``

As an example, clicking `this <#top>`__ will take the reader to the top
of the page.

This was done by having the following lines of code:

::

   For example, clicking [this](#top) will take the reader to the top of the page.

as well as

::

   # NNPDF code and standards documentation <a name="top"></a>

at the top of this file.

In addition, one can link to other pages within the documentation by
``[text](<relative-path-to-md-or-rst-file>.<extension>)``.

One can define “labels” for RestructuredText, which can be referred to
from anywhere, like this:

::

       .. _my-reference-label:

       Section to cross-reference
       --------------------------

       This is the text of the section.

       It refers to the section itself, see :ref:`my-reference-label`.

Such labels can also be defined in Markdown by using ``rst`` syntax
embedded in code markers in markdown:

::

   ```eval_rst
   .. _my-reference-label:
   ```

Labels can be linked to from anywhere using the syntax

::

   [link text](my-reference-label)

for Markdown and

::

   :ref:`my-reference-label`

for RestructuredText, as described in its
`documentation <https://www.sphinx-doc.org/en/master/usage/restructuredtext/roles.html?highlight=cross%20reference#role-ref>`__.

Adding BibTeX references
~~~~~~~~~~~~~~~~~~~~~~~~

The documentation build supports BibTeX references via the
`sphinxcontrib-bibtex extension
<https://github.com/mcmtroffaes/sphinxcontrib-bibtex>`_. Citations in the
BibTeX format are added to the ``references.bib`` file in the Sphinx source
directory. For example a citation like

.. code-block:: bib

    @article{Carrazza:2016htc,
        author = "Carrazza, Stefano and Forte, Stefano and Kassabov, Zahari and Rojo, Juan",
        title = "{Specialized minimal PDFs for optimized LHC calculations}",
        eprint = "1602.00005",
        archivePrefix = "arXiv",
        primaryClass = "hep-ph",
        reportNumber = "CERN-PH-TH-2015-243, TIF-UNIMI-2015-13, OUTP-15-24P",
        doi = "10.1140/epjc/s10052-016-4042-8",
        journal = "Eur. Phys. J. C",
        volume = "76",
        number = "4",
        pages = "205",
        year = "2016"
    }

can be appended to the ``refererences.bib`` file.

References can be added to
RST documents using some variation of the ``cite`` role.  For example
``:cite:p:`<BibTeX ID>``` adds a parenthetical reference, and the above article
can be cited using ``:cite:p:`Carrazza:2016htc```.

Adding indices for modules
~~~~~~~~~~~~~~~~~~~~~~~~~~

Sphinx has the capability of automatically documenting any python
package. It produces these under the ``index`` and ``module index``
sections. The functions and modules are documented using their
corresponding docstrings.

To add a new module to document, add a new line in ``sphinx/Makefile``
under:

::

   %: Makefile
       @if test $@ != "clean"; then 
               sphinx-apidoc -o ./source/modules/validphys ../../validphys2/src/validphys/ ; \
               sphinx-apidoc -o ./source/modules/<MODULE-NAME> <PATH-TO-MODULE>  ;\
       fi
