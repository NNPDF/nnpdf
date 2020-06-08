Important information about runcard layout
==========================================

-  The flag ``fivetheories`` specifies the choice of 5 or
   :math:`\bar{5}` prescription for the case of 5 input theories. You
   must assign a value ``nobar`` or ``bar`` correspondingly. If you do
   not do this, ``validphys`` will give an error.

-  The default behaviour for the 7-point prescription is to use Gavin
   Salam's modification to it. To use the original 7-point prescription
   instead, the ``seventheories`` flag must be set to ``original``.

-  As a fit you must provide the fit for central scales at the relevant
   order.

-  You must also provide the relevant c-factors (EWK, QCD ...).

