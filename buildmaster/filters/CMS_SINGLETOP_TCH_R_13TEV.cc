/*
Measurement of the ratio of the single top quark and single antitop quark cross sections in the t-channel @LHC CMS 13 TeV
There is a single point

LHC-CMS 13 TeV
---------------

Selected events contain one muon and two jets where one of the jets is identified as originating from a b-quark (L=2.2 1/fb)
Archived as: http://arxiv.org/pdf/1610.00678v3.pdf
Published in: Physics Letters B, Vol. 772, 2017, pp 752-776
R_{t-ch.} = 1.81 ± 0.18 (stat) ± 0.15 (syst)

The systematic uncertainty is composed of the following 10 uncertainties (see Table 4 in the paper for the full breakdown and for more details):
Profiled exp. uncert.
Signal modelling
ttbar modelling
W+jets modelling
mu_R/mu_F scale t-channel
mu_R/mu_F scale ttbar
mu_R/mu_F scale tW
mu_R/mu_F scale W+jets
PDF uncert.
Top quark p_T modelling
*/

#include "CMS_SINGLETOP_TCH_R_13TEV.h"

void CMS_SINGLETOP_TCH_R_13TEVFilter::ReadData()
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
  lstream >> fStat[0];       // Absolute statistical uncertainty

  lstream >> fSys[0][0].add; // Absolute total systematic uncertainty
  fSys[0][0].mult = fSys[0][0].add*100/fData[0]; // Multiplicative total systematic uncertainty
  fSys[0][0].type = MULT;
  fSys[0][0].name = "UNCORR";

  f1.close();
}
