/* 
   This file implements the W production subset of the ATLASWZRAP11CC data set.
   This is required to separate CC DY from NC DY.
   Implemented by ERN June 2023.
*/

#include "ATLASWRAP11CC.h"

// Central selection
void ATLASWRAP11CCFilter::ReadData()
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
 
  const int nBins = 2;
  const int ndataWZ[nBins] = {11,22};  // Data thresholds for W+ and W- 
  const double MWZ2[nBins]= {pow(MW,2.0), pow(MW,2.0)};   //Mass squared of W (+ and -)

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


