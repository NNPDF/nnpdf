.. _dataset-naming-convention:


=================================
NNPDF's dataset naming convention
=================================

Each dataset implemented in NNPDF must have a unique name, which is a string
constructed following this [Backus–Naur form]::

  <valid set name> ::= <experiment> "_" <process> "_" <energy>
                     | <experiment> "_" <process> "_" <energy> "_" <extra_information>

  <valid dataset name> ::= <set name> "_" <observable name>

  <experiment> ::= "ATLAS" | "BCDMS" | "CDF" | "CHORUS" | "CMS" | "D0" | "DYE605" | "DYE866" 
                | "DYE906" | "EMC" | "H1" | "HERA" | "LHCB" | "NMC" | "NNPDF" | "NUTEV" | "SLAC"
                | "ZEUS"

  <process> ::= "1JET" | "2JET" | "CC" | "DY" | "INTEG" | "NC" | "PH" | "POS" | "SINGLETOP" 
              | "TTBAR" | "WCHARM" | "WJ" | "WPWM" | "Z0" | "Z0J"

  <energy> ::= <integer> <unit> | <integer> "P" <integer> <unit>  | "NOTFIXED"

  <extra_information> ::= <string>


Experiments
===========

- `ATLAS <https://home.cern/science/experiments/atlas>`_: A Large Toroidal
  Aparatus
- BCDMS: TODO
- CHORUS: TODO
- `CMS <https://home.cern/science/experiments/cms>`_: Compact Muon Solenoid
- E605: TODO
- E866: TODO
- E906: TODO
- EMC: TODO
- `HERA <https://dphep.web.cern.ch/accelerators/hera>`: Hadron Elektron Ring
  Anlage. While technically speaking this is an accelerator, this string is
  used for the combined analyses of H1 and ZEUS.
- `LHCB <https://home.cern/science/experiments/lhcb>`_:
- NMC: TODO
- `NNPDF <https://nnpdf.mi.infn.it/>`_: This experiment name is used for two
  purposes:

  1. for auxiliary datasets needed in the PDF fit, for instance `INTEG` and `POS`
  2. for predictions used in NNPDF papers to study the impact of PDFs in processes not included in its PDF fit
- NUTEV: TODO



Processes
=========

- `1JET`: single-jet inclusive production
- `2JET`: dijet production
- `CC`: DIS charged-current
- `Z0`: lepton-pair production (neutral current off-shell Drell–Yan)
- `H`: on-shell Higgs-boson production
- `HVBF`: production of an on-shell Higgs-boson with two jets (vector-boson
  fusion)
- `INTEG`: auxiliary dataset for integrability constraints; only valid for
  `NNPDF` experiment
- `NC`: DIS neutral-current
- `POS`: auxiliary dataset for positivity constraints; only valid for
  `NNPDF` experiment
- `TTBAR`: top–anti-top production
- `WM`: production of a single negatively-charged lepton (charged current
  off-shell Drell–Yan)
- `WMWP`: production of two opposite-sign different flavor leptons (W-diboson
  production)
- `WP`: production of a single positively-charged lepton (charged current
  off-shell Drell–Yan)
- `WPZ`: production of three leptons (WZ-diboson production)
- `Z0PT`: production of two same-flavor opposite-sign leptons with non-zero
  total transverse momentum (Z-boson pt spectrum)

`Backus–Naur form <https://en.wikipedia.org/wiki/Backus%E2%80%93Naur_form>`_
