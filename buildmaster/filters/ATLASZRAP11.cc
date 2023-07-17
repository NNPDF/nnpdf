/* 
   This file implements the Z production subset of the ATLASWZRAP11CC data set.
   This is required to separate NC DY from CC DY.
   Implemented by ERN June 2023.
*/

#include "ATLASZRAP11CC.h"

// Central selection
void ATLASZRAP11CCFilter::ReadData()
{

  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/ATLASWZRAP11CC/wzrap11.dat";

  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Dummy string  
  std::string dummy;

  // Read the names of systematics
  std::string *sysNames = new std::string[fNSys];
  sysNames[0] = "UNCORR";
  for (int i=0; i<5; i++) f1 >> dummy;
  for (int i=0; i<fNSys; i++) 
    f1 >> sysNames[i];
 
  const int nBins = 3;
  const int ndataWZ[nBins] = {6,18,24};  // Data thresholds for (Z_low, Z_peak, Z_high) respectively
  const double MWZ2[nBins]= {pow(56.0,2.0), pow(91.0,2.0), pow(133.0,2.0)};   //Mass squared of (Z_low, Z_peak, Z_high)

  string line;
  for(int iline=0; iline<23; iline++)
    getline(f1,line);
  
  int low_bin = 0;
  for (int b = 0; b < nBins; b++)
  {
    for (int i = low_bin; i < ndataWZ[b]; i++)
    {
      double etamin, etamax;

      // Kinematics
      f1 >> dummy; f1 >> etamin; f1 >> etamax;
      fKin1[i] = etamin + (etamax - etamin)/2.0;
      fKin2[i] = MWZ2[b];
      fKin3[i] = 7000;

      // Observable
      f1 >> fData[i];
      fData[i] *= 1000; // pb -> fb

      // Statistical errors - percentage with respect the observable
      f1 >> fStat[i];
      fStat[i] *= fData[i]*1e-2;

      // Correlated systematic errors
      for (int l = 0; l < fNSys; l++)
      {
        f1 >> fSys[i][l].mult;
        fSys[i][l].type = MULT;
        fSys[i][l].name = sysNames[l];
      }

      // Additive errors   
      for (int l = 0; l < fNSys; l++)
        fSys[i][l].add = fSys[i][l].mult*fData[i]*1e-2;
    }
    // Update lowest point in bin
    low_bin = ndataWZ[b];
  }

  delete[] sysNames;
  
  f1.close();
}
