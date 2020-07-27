.. _fk_config_variables:

==============================
``FK`` configuration variables
==============================

Table specifying the required elements of the GridInfo ``FK`` header
segment. The Key column specifies the exact format of the Key in the K-V pair
used in the GridInfo segment.

========  =======  ======================  ==================================
Key       Type     Description             Comments
========  =======  ======================  ==================================
SETNAME   String   *SetName*               N/A
HADRONIC  Boolean  Hadronic flag           0 or 1
NDATA     Integer  :math:`N_{\text{dat}}`  Number of data points
NX        Integer  :math:`N_x`             Number of :math:`x`-points in grid
========  =======  ======================  ==================================
