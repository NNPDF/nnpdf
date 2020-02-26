Important information about runcard layout
==========================================

-  **IMPORTANT:** In runcards, theories must be listed according to
   points :math:`(\mu_0, \mu_i)` in the following order: :math:`(0,0)`,
   :math:`(+,0)`, :math:`(-,0)`, :math:`(0,+)`, :math:`(0,-)`,
   :math:`(+,+)`, :math:`(-,-)`, :math:`(+,-)`, :math:`(-,+)`. Here "+"
   refers to doubled scale, "-" to halved scale and "0" to central
   scale.

-  In terms of ``theoryids``, at NLO this corresponds to: 163, 177, 176,
   179, 174, 180, 173, 175, 178. This ensures that the prescriptions
   will be implemented correctly. The ``theoryids`` to input for each
   prescription are listed explicitly in the section on point
   prescriptions below. If an incorrect configuration of ``theoryids``
   is given, you should receive an error message.

-  The flag ``fivetheories`` specifies the choice of 5 or
   :math:`\bar{5}` prescription for the case of 5 input theories. You
   must assign a value ``nobar`` or ``bar`` correspondingly. If you do
   not do this, ``validphys`` will give an error.

-  As a fit you must provide the fit for central scales at the relevant
   order.

-  You must also provide the relevant c-factors (EWK, QCD ...).

