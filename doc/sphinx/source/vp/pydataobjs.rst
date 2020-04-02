.. _pyobjs:

Python based data objects
=========================

Internal data formats such as PDF sets, CommonData, or :ref:`FKTables
<fktables>` files are currently accessed though the `libnnpdf` C++ code
(interfaced trough the SWIG wrappers). However there is a :ref:`project
<https://github.com/NNPDF/nnpdf/issues?q=label%3Adestroyingc%2B%2B+>` underway
to make these resources available in terms of containers standard Python
(particularly numpy arrays and pandas dataframes). The objectives include
simplifying the codebase, increasing the ease of use and enabling more advanced
computation and storage strategies.

Loading FKTables
----------------

Currently only FKTables can be directly without C++ code. This is implemented
in the :py:mod:`validphys.fkarser` module. For example::

    from validphys.fkparser import load_fktable
    from validphys.loader import Loader
    l = Loader()
    fk = l.check_fktable(setname="ATLASTTBARTOT", theoryID=53, cfac=('QCD',))
    res = load_fktable(fk)

results in an object containing all the information needed to compute a
convolution. In particular the ``sigma`` property contains a dataframe
representing the partonic cross-section (including the cfactors).
