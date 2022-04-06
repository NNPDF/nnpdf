.. _mc2hessian:
How to transform a Monte Carlo PDF set into a Hessian PDF set
=============================================================

A Monte Carlo PDF set can be transformed into a Hessian PDF using the
method described in :cite:p:`Carrazza:2016htc` a runcard like the one
below. In this example ``Neig`` are the number of basis eigenvectors of
the Hessian PDF and ``mc2hname`` is an optional argument that can be
added to give a custom name to the Hessian PDF set, if ``mc2hname`` is
not defined, a default name will automatically be given to the Hessian
PDF set. By default ``installgrids`` is ``False``, by setting it to
``True`` we make sure that the Hessian PDF is added to the LHAPDF
folder.

.. code:: yaml

   pdf: NNPDF31_nlo_as_0118

   Neig: 100 # Number of basis eigenvectors of the Hessian PDF set
   Q: 10 # Energy scale at which the MC PDF set is sampled
   mc2hname: MyHessianPDF # Name of the Hessian PDF set
   installgrid: True

   actions_:
     - mc2hessian
