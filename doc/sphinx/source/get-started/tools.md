## Tools, or  "What does any of this have to do with apples?"

*Author: Cameron Voisey, 10/10/2019*

There are many ingredients that go into a PDF fit. For example, one must generate the partonic cross sections for use in the fit, convert these predictions into a format that is suitable for PDF fits, i.e. the format must be suitable for on-the-fly convolutions, and evolve the PDFs according to the DGLAP equations. To perform each of these steps, different codes are used. In this subsection various codes that you will frequently encounter are described.  

### NNPDF specific

validphys

reportengine

### PDF set storage and interpolation

[LHAPDF](https://lhapdf.hepforge.org/) is a C++ library that evaluates PDFs by intepolating the discretised PDF 'grids' that PDF collaborations produce. It also gives its users access to proton and nuclear PDF sets from a variety of PDF collaborations, including NNPDF, MMHT and CTEQ. A list of all currently available PDF sets can be found [here](https://lhapdf.hepforge.org/pdfsets.html). Particle physics programmes such as Monte Carlo event generators will usually be interfaced with LHAPDF, to allow a user to easily specify the PDF set that they wish to use in their calculations.

### Partonic cross section generation

Many programmes exist to produce predictions of partonic cross sections (otherwise known as hard cross sections or matrix elements). Some of general purpose, such as MadGraph5\_aMC@NLO and MCFM, and can therefore compute predictions for a variety of physical processes (e.g. drell-yan and single top production) up to a given order (typically next-to-leading order (NLO) in QCD, with some also computing NLO electroweak (EW) corrections). Others are more specific, such as top++, which makes predictions for top quark pair production up to next-to-next-to-leading order (NNLO). Some of these programmes will be briefly outlined here.

[MadGraph5\_aMC@NLO](https://launchpad.net/mg5amcnlo) is the programme that will be used for most of the future NNPDF calculations of matrix elements. This is in large part due to its ability to compute predictions at NLO in QCD with NLO EW corrections.

#### Other codes

[MCFM](https://mcfm.fnal.gov/) ('Monte Carlo for FeMtobarn processes')

[FEWZ](https://arxiv.org/abs/1011.3540) ('Fully Exclusive W and Z Production') is a programme for calculating (differential) cross sections for the Drell-Yan production of lepton pairs up to NNLO in QCD.

[NLOjet++](http://www.desy.de/~znagy/Site/NLOJet++.html) is a programme that can compute cross sections for a variety of processes up to NLO. The processes include electon-positron annihilation, deep-inelastic scattering (DIS), photoproduction in electron-proton collisions, and a variety of processes in hadron-hadron collisions.

[Top++](http://www.precision.hep.phy.cam.ac.uk/top-plus-plus/) is a programme for computing top quark pair production inclusive cross sections at NNLO in QCD with soft gluon resummation inluded up to next-to-next-to-leading log (NNLL).

### PDF evolution

[APFEL](https://apfel.hepforge.org/) ('A PDF Evolution Library') is the PDF evolution code currently used by the NNPDF Collaboration. In recent years it has been developed along with NNPDF, and so it therefore contains the features and settings required in an NNPDF fit. That is, it includes quark masses in the MSbar scheme, ..., and scale variations up to NLO.

APFEL also produces predictions of deep-inelastic scattering structure functions.

Note that at the time of writing, a code to replace APFEL is being written, which is currently dubbed EKO.

##### Other codes

[Hoppet](https://hoppet.hepforge.org/) ('Higher Order Perturbative Parton Evolution Toolkit') is an alternative PDF evolution code which is capable of evolving unpolarised PDFs to NNLO and linearly polarised PDFs to NLO. The unpolarised evolution includes heavy-quark thresholds in the MSbar scheme.

### Grid generation 

[APPLgrid](https://applgrid.hepforge.org/)

APFELcomb
