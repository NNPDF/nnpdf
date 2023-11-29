.. _org_data_files:

==========================
Organisation of data files
==========================

The ``nnpdf`` code needs to be able to handle a great deal of different
options with regard to the treatment of both experimental data and theoretical
choices. In the code, every effort has been made to keep experimental and
theoretical parameters strictly separate.
In this section we shall specify the layout of the ``nnpdf`` data
directory. It is in this directory that all of the read-only data to be used in
the fit are accessed. The data directory is located in the ``nnpdf`` git
repository, under the path ``validphys/src/validphys2/datafiles``.

Experimental data storage
=========================

The central repository for ``CommonData`` in use by ``nnpdf`` projects is
located in the ``nnpdf`` git repository at

	``validphys/src/validphys2/datafiles/commondata``

where a separate ``CommonData`` file is stored for each *Dataset* with the
filename format

	``DATA_<SETNAME>.dat``

Information on the treatment of systematic uncertainties, provided in
``SYSTYPE`` files, is located in the subdirectory

	``commondata/systypes``

Here several ``SYSTYPE`` files may be supplied for each *Dataset*. The
various options are enumerated by suffix to the filename. The filename format
for ``SYSTYPE`` files is therefore

	``SYSTYPE_<SETNAME>_<SYSID>.dat``

Where the default systematic ID is **DEFAULT**. As an example, consider
the first ``SYSTYPE`` file for the D0ZRAP *Dataset*:

	``SYSTYPE_D0ZRAP_DEFAULT.dat``

Theory lookup table
===================

In order to organise the various different theoretical treatments available, a
lookup table is provided in ``sqlite3`` format. This lookup table can be found
in the ``nnpdf`` repository data directory at:

	``validphys/src/validphys2/datafiles/theory.db``

This file should only be edited in order to add new theory options. It may be
edited with any appropriate ``sqlite3``-supported software. A script is provided to
give a brief overview of the various theory options available. It can be found
at

	``validphys/src/validphys2/datafiles/disp_theory.py``

and should be run without any arguments.

Theory options are enumerated by an integer ``TheoryID``. The parameters of
each theory option are described in the lookup table under the appropriate ID.
The current available parameters are summarised in :ref:`th_parameter_definitions`.

Theory storage
==============

Each theory configuration is stored as a gzip compressed tar archive with
filename format

	``theory_<THEORYID>.tgz``

and is stored at the location specified in the default ``nnprofile.yaml``. For easy
access, they can be downloaded through the ``vp-get`` utility.  Each archive
contains the following directory structure

	| ``theory_X/``
	|	``-cfactor/``
	|	``-compound/``
	|	``-fastkernel/``

Inside the directory ``theory_X/cfactor/`` are stored ``CFACTOR`` files
with the filename format

	``CF_<TYP>_<SETNAME>.dat``

where ``<TYP>`` is a three-letter designation for the source of the C-factor
(e.g. EWK or QCD) and ``<SETNAME>`` is the typical *Dataset* designation.
The directory ``theory_X/compound/`` contains the ``COMPOUND`` files
described earlier, this time with the filename format

	``FK_<SETNAME>-COMPOUND.dat``

Finally the ``FK`` tables themselves are stored in ``theory_X/fastkernel/``
with the filename format

	``FK_<SETNAME>.dat``

Naturally, all of the FastKernel and C-factor files within the directory
``theory_X/`` have been determined with the theoretical parameters specified in
the theory lookup table under ID ``X``.
