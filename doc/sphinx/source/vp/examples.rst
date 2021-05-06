.. _vpexamples:

========
Examples
========

Example validphys runcards can be found
`here <https://github.com/NNPDF/nnpdf/tree/master/validphys2/examples>`_. These are useful when
trying to gain familiarity with how to produce results with validphys or when you want to carry
out a common task, e.g. plotting some PDFs, and you do not want to write the runcard yourself.

It is strongly encouraged to capitalise namespaces, e.g. ``PDFscalespecs`` rather than ``pdfscalespecs``.
This is to avoid confusion between namespaces, which are relevant only to that runcard, and any actions
within ``validphys``, which by convention are lower case.

Here we detail the examples that already exist and list the resources which it is recommended that
you use when writing a new example runcard.

Existing examples
=================

============================= 	===========================    =========================================================
Runcard/folder name		Tutorial			What it does
============================= 	===========================    =========================================================
API_example.ipynb		:ref:`vpapi`			Jupyter notebook example with API	
closure_templates/    		:ref:`tut_closure`  		Running closure tests
cuts_options.yaml             	N/A          			Shows results for different cuts policites
dataspecs.yaml			N/A				Shows how to use ``dataspecs``
data_theory_comparison.yaml	:ref:`datthcomp`		Data theory comparison
export_data.yaml		N/A				Makes tables of experimental data and covariance matrices
generate_a_report.yaml		:ref:`tut_report`		Shows how to generate a report
kiplot.yaml			N/A				Plot kinematic coverage of data
looping_example.yaml		N/A				Shows how to do actions in a loop over resources
mc_gen_example.yaml		N/A				Analysis of pseudodata generation
new_data_specification.yaml	N/A				Shows how to specify data in runcards
pdfdistanceplots.yaml		How to plot PDFs		Distance PDF plots
simple_runcard.yaml		N/A				Simple runcard example
taking_data_from_fit.yaml	N/A				Shows how to take ``theoryids`` and ``pdfs`` from a fit	
theory_covariance/            	:ref:`vptheorycov-index`	Runcards for the ``theorycovariance`` module	
============================= 	===========================    =========================================================

Recommended resources
=====================

The resources that should be used in example runcards, where possible, match the defaults used in
the `tests <https://github.com/NNPDF/nnpdf/blob/master/validphys2/src/validphys/tests/conftest.py#L23>`_.
It is recommended that these are used to avoid situations when the user has to download different
resources, which can be costly in terms of time and memory, to run each example runcard.

The recommended resources are:

===================================  ==============================  ==================================================================
Resource                             ID                              Description
===================================  ==============================  ==================================================================
NLO theoryid                         52                              NNPDF3.1 NLO theory predictions with central scales
NNLO theoryid                        162                             (Low precision) NNPDF3.1 NNLO theory predictions with central scales
NLO theoryid for scale variations 1  163                             Central scales, :math:`k_F = 1, k_R = 1`
NLO theoryid for scale variations 2  173                             :math:`k_F = 0.5, k_R = 0.5`
NLO theoryid for scale variations 3  174                             :math:`k_F = 1, k_R = 0.5`
NLO theoryid for scale variations 4  175                             :math:`k_F = 2, k_R = 0.5`
NLO theoryid for scale variations 5  176                             :math:`k_F = 0.5, k_R = 1`
NLO theoryid for scale variations 6  177                             :math:`k_F = 2, k_R = 1`
NLO theoryid for scale variations 7  178                             :math:`k_F = 0.5, k_R = 2`
NLO theoryid for scale variations 8  179                             :math:`k_F = 1, k_R = 2`
NLO theoryid for scale variations 9  180                             :math:`k_F = 2, k_R = 2`
NLO pdf                              NNPDF31_nlo_as_0118             NNPDF3.1 NLO PDF set with 100 replicas (+ central replica)
NNLO pdf                             NNPDF31_nnlo_as_0118            NNPDF3.1 NNLO PDF set with 100 replicas (+ central replica)
NNLO pdf hessian                     NNPDF31_nnlo_as_0118_hessian    NNPDF3.1 NNLO hessian PDF set generated from replicas
NLO fit                              NNPDF31_nlo_as_0118             NNPDF3.1 NLO fit with 100 replicas (+ central replica)
NNLO fit                             NNPDF31_nnlo_as_0118_DISonly    NNPDF3.1 DIS-only NNLO fit with 95 replicas (+ central replica)
fit                                  191015-mw-001                   n3fit closure test fit with 30 replicas before and after postfit
fit (iterated)                       191015-mw-001_ite2_for_testing  Iteration of 191015-mw-001
===================================  ==============================  ==================================================================
