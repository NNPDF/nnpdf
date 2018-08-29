/*
Measurement of the ratio of the single top quark and single antitop quark cross sections in the t-channel @LHC CMS 7 TeV
There is a single point

LHC-CMS 8 TeV
---------------

emu events with b-tagged jets (L=19.7 1/fb)
Archived as: 1403.7366
Published in: JHEP 06 (2014) 090
R_{t-ch.} = 1.95 ± 0.10 (stat) ± 0.19 (syst)

Here systematic uncertainties are split into 9 different categories (see Table 4 in the paper for the full breakdown and for more details):
JES, JER, MET, and pileup
b-tagging and mis-tag
Leptop reconstruction/trig.
QCD multijet estimation
W+jets, ttbar estimation
Other backgrounds ratio
Signal modeling
PDF uncertainty
Simulation sample size
*/

#include "CMSSTOPTCHRAT8TEV.h"

void CMSSTOPTCHRAT8TEVFilter::ReadData()
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

  //Starting filter
      string line;
      int idum;
      double cme;

      getline(f1,line);
      istringstream lstream(line);
      lstream >> idum >> cme;
      
      fKin1[0] = 0.;
      fKin2[0] = Mt*Mt;          //top mass
      fKin3[0] = cme*1000;       //sqrt(s)

      lstream >> fData[0];       //central value
      lstream >> fStat[0];       //statistical uncertainty

      for (int i = 0; i < (fNSys - 1); i++)
        {
           lstream >> fSys[0][i].mult;
           fSys[0][i].add = fSys[0][i].mult*fData[0]/100; 
           fSys[0][i].type = MULT;
           fSys[0][i].name = "UNCORR";      
        }

  f1.close();

}
