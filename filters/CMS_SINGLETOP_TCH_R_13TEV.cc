/*
Measurement of the ratio of the single top quark and single antitop quark cross sections in the t-channel @LHC CMS 13 TeV
There is a single point

LHC-CMS 13 TeV
---------------

Selected events contain one muon and two jets where on eof the jets is identified as originating from a b-quark (L=2.2 1/fb)
Archived as: http://arxiv.org/pdf/1610.00678v3.pdf
Published in: Physics Letters B, Vol. 772, 2017, pp 752-776
R_{t-ch.} = 1.81 ± 0.18 (stat) ± 0.15 (syst)

Here systematic uncertainties are split into 10 different categories (see Table 4 in the paper for the full breakdown and for more details):
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
  double fstat_percentage;
  double up, down, sigma, data_shift;
  double total_shift;

  getline(f1,line);
  istringstream lstream(line);
  lstream >> idum >> cme;
      
  fKin1[0] = 0.;
  fKin2[0] = Mt*Mt;          // Top mass
  fKin3[0] = cme*1000;       // Sqrt(s)

  lstream >> fData[0];       // Central value
  lstream >> fstat_percentage; // Statistical (percentage) uncertainty
  fStat[0] = fstat_percentage*fData[0]/100; // Convert percentage uncertainty to absolute uncertainty and store

  for (int i = 0; i < 3; i++) // Profiled exp. uncert., signal modelling, ttbar modelling
    {
       lstream >> fSys[0][i].mult;
       fSys[0][i].add = fSys[0][i].mult*fData[0]/100; 
       fSys[0][i].type = MULT;
       fSys[0][i].name = "UNCORR";      
    }

  for (int i = 3; i < 9; i++) // Symmetrise W+jets modelling uncert., scale variations uncerts., and PDF uncert.
    {
       lstream >> up >> down;
       symmetriseErrors(up, down, &sigma, &data_shift);
       fSys[0][i].mult = sigma;
       fSys[0][i].add = fSys[0][i].mult*fData[0]/100;
       fSys[0][i].type = MULT;
       fSys[0][i].name = "UNCORR";
       total_shift += data_shift;
    }

  lstream >> fSys[0][9].mult; // Top quark p_T modelling
  fSys[0][9].add = fSys[0][9].mult*fData[0]/100;
  fSys[0][9].type = MULT;
  fSys[0][9].name = "UNCORR";

  fData[0] *= (1.0 + total_shift*0.01); // Shift of central value due to asymmetric uncertainties

  f1.close();

}
