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

To add a new section to the documentation, create an appropriately named
directory in the ``sphinx/source/`` directory. Inside the new directory,
add all relevant documentation in restructured text
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
