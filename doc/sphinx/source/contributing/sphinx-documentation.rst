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

New documentation needs to be formatted as `reStructuredText
<https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_.


To add a new section to the documentation, create an appropriately named
directory in the ``sphinx/source/`` directory. In addition to these files,
create an ``index.rst`` file containing:

::

   Chapter Name
   ============

   .. toctree::
      :maxdepth: 1

      ./file1
      ./file2

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

Useful reStructuredText Tools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Various `reStructuredText
<http://docutils.sourceforge.net/docs/user/rst/quickref.html>`__ cheatsheets
exist online. We list some markup constructs useful to the NNPDF documentation.

A :math:`\LaTeX` block can be generated using

::

   .. math::

      \frac{1}{2}

while inline maths is generated using

::

   :math:`\frac{1}{2}`

with attention being brought to the backticks.

One can define “labels” for RestructuredText, which can be referred to
from anywhere, like this:

::

       .. _my-reference-label:

       Section to cross-reference
       --------------------------

       This is the text of the section.

       It refers to the section itself, see :ref:`my-reference-label`.


Labels can be linked to from anywhere using the syntax

::

   :ref:`my-reference-label`

as described in the
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
