Partonic cross section generation
=================================

Many programmes exist to evaluate partonic cross sections. Some are
general purpose, such as MadGraph5_aMC@NLO and MCFM, in that they
compute predictions for a variety of physical processes (e.g. drell-yan
production, single top production, …) up to a given order. Others are
more specific, such as top++, which makes predictions for top quark pair
production only. Some of these programmes will be briefly outlined here.
Note that to produce predictions at NNLO in QCD, which is the highest
order used in NNPDF fits, one usually produces APPLgrids at NLO in QCD,
and then supplements these with NNLO QCD corrections which are computed
with a code with NNLO capabilities. These C-factors are often provided
to the collaboration by external parties, rather than the code being run
in-house.

`MadGraph5_aMC@NLO <https://launchpad.net/mg5amcnlo>`__ is the programme
that will be used for most of the future NNPDF calculations of partonic
cross sections. This is in large part due to its ability to compute
predictions at NLO in QCD with NLO EW corrections. To generate APPLgrids
from MadGraph5_aMC@NLO, one can use
`aMCfast <https://amcfast.hepforge.org/>`__, which interfaces between
the two formats.

Other codes
-----------

`MCFM <https://mcfm.fnal.gov/>`__ (‘Monte Carlo for FeMtobarn
processes’) is an alternative programme to MadGraph5_aMC@NLO, which
instead uses mcfm-bridge as an interface to generate APPLgrids.

`FEWZ <https://arxiv.org/abs/1011.3540>`__ (‘Fully Exclusive W and Z
Production’) is a programme for calculating (differential) cross
sections for the Drell-Yan production of lepton pairs up to NNLO in QCD.

`NLOjet++ <http://www.desy.de/~znagy/Site/NLOJet++.html>`__ is a
programme that can compute cross sections for a variety of processes up
to NLO. The processes include electron-positron annihilation,
deep-inelastic scattering (DIS), photoproduction in electron-proton
collisions, and a variety of processes in hadron-hadron collisions.

`Top++ <http://www.precision.hep.phy.cam.ac.uk/top-plus-plus/>`__ is a
programme for computing top quark pair production inclusive cross
sections at NNLO in QCD with soft gluon resummation included up to
next-to-next-to-leading log (NNLL).
