.. _run-qed-fit:

==========================
How to run a QED fit
==========================

It is possible to perform a QED fit adding the key `fiatlux` to the runcard. In this way
a photon PDF will be generated using the `FiatLux` public library that implements the `LuxQED`
(see :cite:p:`Manohar:2016nzj` and :cite:p:`Manohar:2017eqh`) approach.
The parameters to be added are the following:

.. code-block:: yaml

    fiatlux:
      luxset: NNPDF40_nnlo_as_01180
      additional_errors: true
      luxseed: 1234567890

`luxset` is the PDF set used to generate the photon PDF with `FiatLux <https://github.com/scarrazza/fiatlux/>`.
The code wil generate as much photon replicas as the number of replicas contained in the `luxset`. Therefore, if the user
tries to generate a replica with ID higher than the maximum ID of the `luxset`, the code will
raise an error. Moreover, being the `LuxQED` approach an iterated prcedure, and given that some replicas
do not pass the `postfit` selection criteria, the user should make sure that the number of replicas in
the `luxset` is high enough such that in the final iteration there will be a number of replicas 
higher than the final replicas desired.
`additional_errors` is the parameter that switches on and off the additional errors of the `LuxQED` approach,
while `luxseed` is the seed used to generate such errors.
This errors should be switched on only in the very last iteration of the procedure.

Whenever the photon PDF is generated, it will remain constant during the fit and will enter in the `MSR`.      
