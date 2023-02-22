/**
 
In this experiment we include the most precise data
on top quark production cross sections from Tevatron, ATLAS and CMS,
in the latter case at 7 TeV and 8 TeV

Tevatron
---------

Combination (8.8 fb^{-1})
http://arxiv.org/abs/1309.7570
sigma = 7.60 +- 0.20 (stat) +- 0.29 (syst) +- 0.21 (lumi) pb

LHC 7 TeV
--------

ATLAS dileptons (0.7 1/fb)
JHEP12 (2012) 059 
173 +- 6 (stat) +14-11 (sys) +8-7 (lumi) pb

ATLAS lepton+jet (0.7 1/fb)
ATLAS-CONF-2011-121
179.0 ± 4 (stat) +- 9 (syst) ± 6.6 (lumi) pb

CMS dileptons (L=2.3/fb)
JHEP 11 (2012) 067 
162 ± 2 (stat ) ± 5 (sys) ± 4 (lumi) pb

CMS lepton+jets (L=2.2-2.3/fb)
PLB 720 (2013) 83 
158 ± 2 (stat) ± 10 (sys) ± 4 (lumi) pb


LHC 8 TeV
---------

ATLAS dilepton (20.3 1/fb)
ATLAS-CONF-2013-097
237.7 ± 1.7 (stat) ± 7.4 (syst) ± 7.4 (lumi) ± 4.0 (beam energy) pb 

CMS dilepton  (L=5.3/fb) 
JHEP 02 (2014) 024
239 ± 2 (stat) ± 11 (syst) ± 6 (lumi) pb

So we have a total of 7 data points for this observable

The raw data is available in ttbar_totxsec.data

The format is

index  xsec  stat  sys lumi beamenergy

Everything is in picobarns

Tevatron data is not included for the time being

*/

#include "TOP.h"

void TTBARTOTFilter::ReadData()
{

  // PDG average of the top quark pole mass
  double const mt = 173.3;
  
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
	   << fSetName << "/ttbar_totxsec.data";
  f1.open(datafile.str().c_str(), ios::in);
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
  }
  string line;
  for (int i = 0; i < fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      
      fKin1[i] = 0.0;
      fKin2[i] = mt*mt;
      fKin3[i] = 0.0;
      
      int idum=0;
      lstream >> idum; 
      if(idum!=(i+1)){
	cout<<"Error in TOP.cc"<<endl;
	cout<<"idum = "<<idum<<endl;
	exit(-10);
      }
      
      lstream >> fData[i]; // Central values
      lstream >> fStat[i]; // Statistical error
      lstream >> fSys[i][0].add; // Total systematic uncertainties      
      lstream >> fSys[i][1].add;  // Lumi
      lstream >> fSys[i][2].add; // Beam energy uncertainty 

      fSys[i][0].mult = fSys[i][0].add*100/fData[i];
      fSys[i][0].type = MULT;
      fSys[i][0].name = "UNCORR";  // for time being treat as uncorrelated
      
      fSys[i][1].mult = fSys[i][1].add*100/fData[i];
      fSys[i][1].type = MULT;
      fSys[i][1].name = "UNCORR";  // for time being treat as uncorrelated 
      
      fSys[i][2].mult = fSys[i][2].add*100/fData[i];
      fSys[i][2].type = MULT;
      fSys[i][2].name = "UNCORR";  // for time being treat as uncorrelated   
    }
  
  f1.close();
}

