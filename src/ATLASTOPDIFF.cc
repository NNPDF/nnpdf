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

#include "ATLASTOPDIFF.h"

// Raw data available here
// https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/TOPQ-2015-06/#auxstuff

//-------------------------------------------------------------------------------

// ATLAS top quark differential distributions
// 8 TeV, top quark pt distribution
void  ATLASTOPDIFF8TEVTPTFilter::ReadData()
{
  // Opening files
  fstream f1;

  //
  // The raw data for (1/sigma) dsigma/dpt is taken from
  // https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/TOPQ-2015-06/tabaux_025.pdf
  // Format
  // pt_min  pt_max   (1/sigma)dsigma/dpt   Total_Unc(%)
  //
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
	   << fSetName << "/ATLASTOPDIFF8TEVTPT.dat";
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
            
      int idum=0;
      lstream >> idum; 
      if(idum!=(i+1)){
	cout<<"Error in ATLASTOPDIFF.cc"<<endl;
	cout<<"idum = "<<idum<<endl;
	exit(-10);
      }

      double pt_min=0.0;
      lstream >> pt_min;
      double pt_max=0.0;
      lstream >> pt_max;

      fKin1[i] = 0.0;
      fKin2[i] = 0.5*(pt_min + pt_max); // <pt_top> as characteristic kin variable
      fKin3[i] = 0.0;
      
      lstream >> fData[i]; // Central values
      lstream >> fStat[i]; // Total statistical+systematic error (in percent)
      
    }
  
  f1.close();
}

///////////////////////////////////////////////////////////

