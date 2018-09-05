/*
Measurement of the ratio of single top quark and single antitop quark cross sections in the t-channel @LHC ATLAS 13 TeV
There is a single point

LHC-ATLAS 13 TeV
----------------

Selected events require one charged lepton (electron or muon), missing tranverse momentum, and two jets with high transverse momentum,
exactly one of which is required to be b-tagged.
Archived as: https://arxiv.org/pdf/1609.03920v2.pdf
Published in: JHEP 04 (2017) 086
R_t = 1.72 ± 0.09 (stat.) ± 0.18 (syst.)
Here systematics are split into 14 categories (see Table 4 in the paper for the full breakdown and for more details):
Monte Carlo statistics
Muon uncertainties
Electron uncertainties
JES
Jet energy resolution
c-tagging efficiency
Pile-up reweighting
tq parton shower generator
tq NLO matching
tq radiation
ttbar, Wt, t bbar + tbar b parton shower generator
ttbar, Wt, t bbar + tbar b NLO matching
ttbar, Wt, t bbar + tbar b radiation
Multijet normalisation
*/

#include "ATLASSTOPTCHRAT13TEV.h"

void ATLASSTOPTCHRAT13TEVFilter::ReadData()
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
  lstream >> fstat_percentage;       // Statistical (percentage) uncertainty
  fStat[0] = fstat_percentage*fData[0]/100; // Convert percentage uncertainty to absolute uncertainty and store

  for (int i = 0; i < fNSys; i++)
    {
       lstream >> fSys[0][i].mult;
       fSys[0][i].add = fSys[0][i].mult*fData[0]/100;
       fSys[0][i].type = MULT;
       fSys[0][i].name = "UNCORR";
    }

  f1.close();

}
