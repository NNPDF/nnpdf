.. _nnpdf-data-cli:

Data utilities: ``nnpdf-data``
==============================

``nnpdf-data`` provides command line utilities for datasets available to ``n3fit`` and ``validphys``.
It scans the data paths configured in the NNPDF profile and reads the installed CommonData metadata.
It can be used to prepare ``dataset_inputs`` blocks and to generate LaTeX dataset tables, including BibTeX references from a fit runcard.

.. important::
   If you are implementing or changing datasets, see :ref:`commondata`.


Listing datasets
----------------

Use ``list`` to inspect the available datasets::

  nnpdf-data list ATLAS*

The filter is a glob pattern. Brace expansion is also supported::

  nnpdf-data list '{ATLAS,CMS}*'

The output can be sorted by experiment, process, or the ``nnpdf31-process`` key which is used for theory correlations::

  nnpdf-data list ATLAS* --experiment
  nnpdf-data list ATLAS* --process
  nnpdf-data list ATLAS* --nnpdf31-process

Use ``--yaml`` to print entries ready for ``dataset_inputs``::

  nnpdf-data list 'ATLAS_1JET*' --yaml


Generating dataset tables
-------------------------

Use ``latex`` to read the ``dataset_inputs`` entry of a runcard and print a
LaTeX table::

  nnpdf-data latex runcard.yml

By default, BibTeX entries are printed to stdout. To write them to a file instead::

  nnpdf-data latex runcard.yml --output-bib datasets.bib

The table can be sorted by experiment, process, or NNPDF3.1 process. With
``--group-by``, one table is produced for each selected group::

  nnpdf-data latex runcard.yml --process --group-by
