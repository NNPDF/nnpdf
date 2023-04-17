============
Introduction
============

In the ``nnpdf++`` project, data files used by the code may be grouped into
two categories, theory and experiment. Experimental data and the information
pertaining to the treatment of systematic errors are held in ``CommonData``
and ``SYSTYPE`` files. ``FK`` tables, ``COMPOUND`` and ``CFACTOR`` files
store the precomputed information for use when calculating theoretical
predictions corresponding to information held in the equivalent ``CommonData``
file. In this section the file formats and naming conventions for these files
will be detailed, along with the directory structure employed by the
``nnpdf++`` code.

For NNPDF3.1 and later fits, a considerably larger number of theory options will
be explored than in previous determinations. In NNPDF3.0 the main theory
variations used were perturbative order, value of the strong coupling and the
number of active flavours in the VFNS. For NNPDF3.1 and later, it has been necessary to
accommodate variations in additional parameters, such as treatments of the heavy
quark mass (pole vs MS-bar), scale variations, intrinsic charm, resummation
effects etc. The book-keeping used to enable efficient variations of the
theoretical treatment used in fits post-3.0 will therefore also be outlined
here.

This section will begin by detailing the specifications for the file formats
used by the code, first with the experimental data file formats and layouts in
:ref:`exp_data_files` and secondly with the file formats used for
theoretical predictions in :ref:`th_data_files`. Finally the organisation of
these files within the ``nnpdf++`` structure will be described in
:ref:`org_data_files`.

Important definitions
=====================

In order to clarify the later description, here are a few important
terminological points to note.

*Dataset* vs *Experiment*
-------------------------

When referring to a collection of data points two words are used in the
``nnpdf++`` code which have specific meanings. *Dataset* refers to the result
of a specific measurement, typically associated with a single experimental paper
and corresponds to the *DataSet* class in the ``nnpdf++`` code.
*Experiment* refers to a collection of *Datasets* which are associated
by experimental cross-correlations. For example, the ATLAS 2010 R=0.4 inclusive
jet measurement and the ATLAS 2011 high-mass Drell-Yan measurement are both
examples of *Datasets* as used in the NNPDF3.0 analysis. Both of these
datasets are grouped into the ATLAS *Experiment* as they have systematic
uncertainties that are cross-correlated with each other. In this document, when
using these terms in this sense, they will be italicised for clarity.

Note however that the concept of an *Experiment* is being phased out in the NNPDF
code. For more information on this see :ref:`data_specification`.

*Dataset* and *Experiment* names
--------------------------------

When referred to, the *Dataset* and *Experiment* names refer to the
short identifying string used in the code for each *Dataset* and
*Experiment*.  For example, the *Dataset* name for the aforementioned
ATLAS 2010 inclusive jet measurement with R=0.4 is ATLASR04JETS36PB.

New dataset naming conventions
------------------------------

See :ref:`dataset_naming_convention` for a definition of how datasets should be
named.
