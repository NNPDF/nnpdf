.. _th_parameter_info:

Looking up the parameters of a theory
=====================================

The parameters for all of the theories can be found in the ``theory.db``
file, located in the ``<NNPDF install location>/share/NNPDF/data``
directory. This is an sqlite database file. The information contained
within this file can also be viewed `within the
docs <theory-indexes>`__.

The tools required to extract the parameters for a given theory are
already in the validphys framework:

.. code:: python

   >>> from validphys.loader import Loader

   >>> l = Loader()
   >>> # replace 53 with the relevant theoryID
   >>> info_dict = l.check_theoryinfo(53)
   >>> print(info_dict)
   {'ID': 53, 'PTO': 2, 'FNS': 'FONLL-C', 'DAMP': 0, 'IC': 1, 'ModEv': 'TRN',
   'XIR': 1.0, 'XIF': 1.0, 'NfFF': 5, 'MaxNfAs': 5, 'MaxNfPdf': 5, 'Q0': 1.65,
   'alphas': 0.118, 'Qref': 91.2, 'QED': 0, 'alphaqed': 0.007496252,
   'Qedref': 1.777, 'SxRes': 0, 'SxOrd': 'LL', 'HQ': 'POLE', 'mc': 1.51,
   'Qmc': 1.51, 'kcThr': 1.0, 'mb': 4.92, 'Qmb': 4.92, 'kbThr': 1.0, 'mt': 172.5,
   'Qmt': 172.5, 'ktThr': 1.0, 'CKM': '0.97428 0.22530 0.003470 0.22520 0.97345 0.041000 0.00862 0.04030 0.999152',
   'MZ': 91.1876, 'MW': 80.398, 'GF': 1.1663787e-05, 'SIN2TW': 0.23126, 'TMC': 1,
   'MP': 0.938, 'Comments': 'NNPDF3.1 NNLO central', 'global_nx': 0, 'EScaleVar': 1}

printing a dictionary in a python terminal is a bit cumbersome, the
above method for checking theory parameters has been added to a command
line script ``vp-checktheory`` which essentially does the same thing but
prints the table in a nicer format.

The usage is simple:

::

   $ vp-checktheory 53
                                             Info for theory 53
   ID                                                        53
   PTO                                                        2
   FNS                                                  FONLL-C
   DAMP                                                       0
   IC                                                         1
   ModEv                                                    TRN
   XIR                                                        1
   XIF                                                        1
   NfFF                                                       5
   MaxNfAs                                                    5
   MaxNfPdf                                                   5
   Q0                                                      1.65
   alphas                                                 0.118
   Qref      Qref                                                    91.2
   QED                                                        0
   alphaqed                                          0.00749625
   Qedref                                                 1.777
   SxRes                                                      0
   SxOrd                                                     LL
   HQ                                                      POLE
   mc                                                      1.51
   Qmc                                                     1.51
   kcThr                                                      1
   mb                                                      4.92
   Qmb                                                     4.92
   kbThr                                                      1
   mt                                                     172.5
   Qmt                                                    172.5
   ktThr                                                      1
   CKM        0.97428 0.22530 0.003470 0.22520 0.97345 0.041...
   MZ                                                   91.1876
   MW                                                    80.398
   GF                                               1.16638e-05
   SIN2TW                                               0.23126
   TMC                                                        1
   MP                                                     0.938
   Comments                               NNPDF3.1 NNLO central
   global_nx                                                  0
   EScaleVar                                                  1

Some of the entries of the table have been truncated to fit in the
terminal, for example the CKM matrix elements. If you want to see the
full output or wish to keep a permanent copy of the table then you can
use the command line option ``-d``, which stands for ``--dumptable``
which will save the table in your current working directory:

::

   $ vp-checktheory 53 -d
   [INFO]: Saving info table to theory_53_info.csv
                                             Info for theory 53
   ID                                                        53
   PTO                                                        2
   FNS                                                  FONLL-C
   DAMP                                                       0
   IC                                                         1
   ModEv                                                    TRN
   XIR                                                        1
   XIF                                                        1
   NfFF                                                       5
   MaxNfAs                                                    5
   MaxNfPdf                                                   5
   Q0                                                      1.65
   alphas                                                 0.118
   Qref                                                    91.2
   QED                                                        0
   alphaqed                                          0.00749625
   Qedref                                                 1.777
   SxRes                                                      0
   SxOrd                                                     LL
   HQ                                                      POLE
   mc                                                      1.51
   Qmc                                                     1.51
   kcThr                                                      1
   mb                                                      4.92
   Qmb                                                     4.92
   kbThr                                                      1
   mt                                                     172.5
   Qmt                                                    172.5
   ktThr                                                      1
   CKM        0.97428 0.22530 0.003470 0.22520 0.97345 0.041...
   MZ                                                   91.1876
   MW                                                    80.398
   GF                                               1.16638e-05
   SIN2TW                                               0.23126
   TMC                                                        1
   MP                                                     0.938
   Comments                               NNPDF3.1 NNLO central
   global_nx                                                  0
   EScaleVar                                                  1
   $ ls
   theory_53_info.csv

user can also parse the ``theoryid`` from a fit

::

   $vp-checktheory --fit FIT

where ``FIT`` is a valid fit name. If the fit cannot be found locally,
the script will attempt to download it.

The parameters in the above are defined
`here <./theoryparamsdefinitions>`__.
