/**
 
This file describes the buildmaster implementation of the
differential distributions for top quark pair production at 8 TeV
from ATLAS and CMS

ATLAS => http://arxiv.org/abs/1511.04716
CMS => http://arxiv.org/abs/1505.04480

There are two different ways to include this data, with
normalized or with absolute distributions

By modifying this file, one can select the normalized or the absolute
distributions

*/

#include "TOPDIFF.h"

// Raw data available here
// https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/TOPQ-2015-06/#auxstuff

// ATLAS top quark differential distributions
// 8 TeV, top quark pt distribution
void  ATLASTOPDIFF8TEVTPTFilter::ReadData()
{
  // Opening files
  fstream f1;

  // The raw data for (1/sigma) dsigma/dpt is taken from
  // https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/TOPQ-2015-06/tabaux_025.pdf
  
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
      fKin2[i] = mt;
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

