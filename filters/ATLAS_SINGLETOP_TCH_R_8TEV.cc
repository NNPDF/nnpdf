/*
Measurement of the ratio of the single top quark and single antitop quark cross sections in the t-channel @LHC ATLAS 8 TeV
There is a single point

LHC-ATLAS 8 TeV
---------------

Selected events contain exactly one electron or muon, exactly two jets (exactly one of which must be b-tagged), and E_T^{miss} > 30 GeV. 
For more details of the event selection see Section 5 of the paper.
Archived as: https://arxiv.org/pdf/1702.02859v3.pdf
Published in: Eur. Phys. J. C 77 (2017) 531
R_t = 1.72 ± 0.05 (stat.) ± 0.07 (exp.)

The experimental systematic uncertainty (exp.) has the following significant contributions (see Table 5 in the paper for the full breakdown and for more details):
Monte Carlo statistics
Background modelling
Jet reconstruction
E_T^{miss} modelling
tq (tbar q) NLO matching
tq (tbar q) scale variations
ttbar NLO matching
ttbar parton shower
PDF
*/

#include "ATLAS_SINGLETOP_TCH_R_8TEV.h"

void ATLAS_SINGLETOP_TCH_R_8TEVFilter::ReadData()
{
  // Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/" << fSetName << "/" << fSetName << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  // Starting filter
  string line;
  int idum;
  double cme;

  getline(f1,line);
  istringstream lstream(line);
  lstream >> idum >> cme;

  fKin1[0] = 0.;
  fKin2[0] = Mt*Mt;          // Top mass
  fKin3[0] = cme*1000;       // Sqrt(s)

  lstream >> fData[0];       // Central value
  lstream >> fStat[0];       // Statistical (percentage) uncertainty

  lstream >> fSys[0][0].add; // Absolute total systematic uncertainty
  fSys[0][0].mult = fSys[0][0].add*100/fData[0]; // Multiplicative total systematic uncertainty
  fSys[0][0].type = MULT;
  fSys[0][0].name = "UNCORR";
}
