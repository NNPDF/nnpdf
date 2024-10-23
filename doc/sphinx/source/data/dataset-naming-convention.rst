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

  <process> ::= "1JET" | "2JET" | "CC" | "DY" | "INTEG" | "NC" | "PH" | "POS" | "SHP" | "SINGLETOP"
              | "TTBAR" | "WCHARM" | "WJ" | "WPWM" | "Z0" | "Z0J"

  <energy> ::= <integer> <unit> | <integer> "P" <integer> <unit>  | "NOTFIXED"

  <extra_information> ::= <string>


Experiments
===========

- `ATHENA`: A Totally Hermetic Electron Nucleus Apparatus 
  proposed for IP6 at the Electron-Ion Collider
- `ATLAS <https://home.cern/science/experiments/atlas>`_: A Large Toroidal
  Aparatus
- BCDMS: Muon-proton collider
- `CDF <https://www.fnal.gov/pub/tevatron/experiments/cdf.html>`_: Collider Detector at Fermilab
- CHORUS: CERN Hybrid Oscillation Research apparatUS
- `CMS <https://home.cern/science/experiments/cms>`_: Compact Muon Solenoid
- `COMPASS <https://home.cern/science/experiments/compass>`_: Common Muon and Proton Apparatus 
  for Structure and Spectroscopy
- DO: Tevatron collider experiment 
- `E142 <https://inspirehep.net/experiments/1108817>`_: measurement of the neutron 
  spin dependent structure function at SLAC
- `E143 <https://inspirehep.net/experiments/1108679>`_: measurement of the nucleon 
  spin structure in the end stattion at SLAC
- `E154 <https://inspirehep.net/experiments/1108588>`_: successor of E142 
  using a polarized HE3 target
- `E155 <https://inspirehep.net/experiments/1108587>`_: successor of E143
- E605: Fixed-target experiment at Fermilab
- E866: NuSea, also Fermilab.
- E906: SeaQuest, successor of E866
- EMC: European Muon Collaboration
- `H1 <https://h1.desy.de/>`_: HERA experiment
- `HERA <https://dphep.web.cern.ch/accelerators/hera>`: Hadron Elektron Ring
  Anlage. While technically speaking this is an accelerator, this string is
  used for the combined analyses of H1 and ZEUS.
- HERMES: HERA experiment on nucleon spin properties
- `E06-010 <https://hallaweb.jlab.org/experiment/transversity/>`_: JLAB experiment
- `E97-110 <https://hallaweb.jlab.org/experiment/E97-110/>`_: JLAB experiment
- EG1b: JLAB experiment from the CLAS collaboration
- EG1-DVCS: JLAB experiment from the CLAS collaboration
- `LHCB <https://home.cern/science/experiments/lhcb>`_:
- NMC: New Muon Collaboration, successor of EMC.
- `NNPDF <https://nnpdf.mi.infn.it/>`_: This experiment name is used for two
  purposes:
  1. for auxiliary datasets needed in the PDF fit, for instance `INTEG` and `POS`
  2. for predictions used in NNPDF papers to study the impact of PDFs in processes not included in its PDF fit
- NUTEV: heavy target neutrino detector
- `PHENIX <https://www.bnl.gov/rhic/phenix.php>`_: Pioneering High Energy Nuclear 
  Interaction eXperiment at RHIC
- SLAC: Stanford Linear Accelerator Center
- SMC: Spin Muon Collaboration
- SMCSX: Low-x Spin Muon Collaboration
- `STAR <https://www.bnl.gov/rhic/star.php>`_: Solenoidal Tracker at RHIC
- ZEUS: HERA experiment



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
- `SHP`: single hadron production
- `TTBAR`: top–anti-top production
- `WM`: production of a single negatively-charged lepton (charged current
  off-shell Drell–Yan)
- `WMWP`: production of two opposite-sign different flavor leptons (W-diboson
  production)
- `WP`: production of a single positively-charged lepton (charged current
  off-shell Drell–Yan)
- `Z0J`: production of two same-flavor opposite-sign leptons with non-zero
  total transverse momentum (Z-boson pt spectrum)
- `WJ`: same for W (W-boson pt spectrum)

`Backus–Naur form <https://en.wikipedia.org/wiki/Backus%E2%80%93Naur_form>`_
