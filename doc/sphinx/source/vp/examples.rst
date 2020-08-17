========
Examples
========

Example validphys runcards can be found
`here <https://github.com/NNPDF/nnpdf/tree/master/validphys2/examples>`_. These are useful when
trying to gain familiarity with how to produce results with validphys or when you want to carry
out a common task, e.g. plotting some PDFs, and you do not want to write the runcard yourself.

Here we detail the examples that already exist and list the resources which it is recommended that
you use when writing a new example runcard.

Recommended resources
=====================

The resources that should be used in example runcards, where possible, match the defaults used in
the `tests <https://github.com/NNPDF/nnpdf/blob/master/validphys2/src/validphys/tests/conftest.py#L23>`_.
It is recommended that these are used to avoid situations when the user has to download different
resources, which can be costly in terms of time and memory, to run each example runcard.

The recommended resources are:

==============  ==============================  ==================================================================
Resource        ID                              Description
==============  ==============================  ==================================================================
theoryid        162                             Low precision NNPDF3.1 NNLO theory predictions with central scales
pdf             NNPDF31_nnlo_as_0118            NNPDF3.1 NNLO PDF set with 100 replicas (+ central replica)
fit             191015-mw-001                   n3fit closure test fit with 30 replicas before and after postfit
fit (iterated)  191015-mw-001_ite2_for_testing  Iteration of 191015-mw-001
==============  ==============================  ==================================================================

