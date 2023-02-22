/*
Measurement of the ratio of the single top quark and single antitop quark cross sections in the t-channel @LHC ATLAS 7 TeV
There is a single point

LHC-ATLAS 7 TeV
---------------

Selected events contain one charged lepton, large missing transverse momentum, and two or three jets (L = 4.59 1/fb)
Archived as: https://arxiv.org/pdf/1406.7844v2.pdf
Published in: Physics Review D 90, 112006 (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.90.112006)
R_{t} = 2.04 ± 0.13 (stat.) ± 0.12 (syst.)

The total systematic uncertainty (syst.) has the following contributions (see Table III in the paper for the full breakdown and for more details):
Monte Carlo statistical
Multijet normalization
Other background normalization
JES \eta intercalibration
JES flavor composition
JES flavor response
Jet energy resolution
E_T^{miss} modeling
Lepton uncertainties
PDF
ttbar generator + parton shower
ttbar ISR/FSR
Luminosity
*/

#include "ATLAS_SINGLETOP_TCH_R_7TEV.h"

void ATLAS_SINGLETOP_TCH_R_7TEVFilter::ReadData()
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
  double fstat_percentage;

  getline(f1,line);
  istringstream lstream(line);
  lstream >> idum >> cme;

  fKin1[0] = 0.;
  fKin2[0] = Mt*Mt;          // Top mass
  fKin3[0] = cme*1000;       // Sqrt(s)

  lstream >> fData[0];       // Central value
  lstream >> fstat_percentage; // Statistical (percentage) uncertainty
  fStat[0] = fstat_percentage*fData[0]/100; // Convert percentage uncertainty to absolute uncertainty

  for (int i = 0; i < fNSys - 1; i++)
    {
      lstream >> fSys[0][i].mult; // Percentage systematic uncertainties
      fSys[0][i].add = fSys[0][i].mult*fData[0]/100; // Absolute systematic uncertainties
      fSys[0][i].type = MULT;
      fSys[0][i].name = "UNCORR";
    }

  lstream >> fSys[0][12].mult;  // Percentage luminosity uncertainty
  fSys[0][12].add = fSys[0][12].mult*fData[0]/100; // Absolute luminosity uncertainty
  fSys[0][12].type = MULT;
  fSys[0][12].name = "ATLASLUMI11";
}
