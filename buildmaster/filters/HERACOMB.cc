/**
 * HERACOMB.cc
 * Implementation of final combined HERA inclusive xsec data by nh (to be crosschecked with lr)
 */

#include "HERACOMB.h"

void HERACOMBFilter::ReadData()
{
  // Common parameters
  const int nCorSys = 162;  // Number of correlated systematics
  const int nProcSys = 7;   // Number of procedural systematics

  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << "HERACOMB/" << fSetName <<".dat";

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
  for (int i=0; i<6; i++) f1 >> dummy;
  for (int i=0; i<nCorSys; i++) f1 >> sysNames[i+1];
    
  f1 >> dummy;
  for (int i=0; i<nProcSys; i++) f1 >> sysNames[i+1+nCorSys];

  for (int i = 0; i < fNData; i++)
  {
    f1 >> fKin2[i]; // Q^2
    f1 >> fKin1[i]; // x
    f1 >> fKin3[i]; // y

    // Observable
    f1 >> fData[i];
    
    // Statistical errors - percentage with respect the observable
    f1 >> fStat[i];
    fStat[i] = fStat[i]*fData[i]*1e-2;

    // Uncorrelated systematics
    f1 >> fSys[i][0].mult;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "UNCORR";
    
    // Correlated systematic errors
    for (int l = 1; l < nCorSys+1; l++)
    {
      f1 >> fSys[i][l].mult;
      fSys[i][l].type = MULT;
      fSys[i][l].name = "HC_"+sysNames[l];
    }

    f1 >> dummy; // Skip dummy

    // Procedural systematic errors
    for (int l = nCorSys+1; l < fNSys; l++)
    {
      f1 >> fSys[i][l].mult;
      fSys[i][l].type = MULT;
      fSys[i][l].name = "HC_"+sysNames[l];
    }

    // Additive errors   
    for (int l = 0; l < fNSys; l++)
      fSys[i][l].add = fSys[i][l].mult*fData[i]*1e-2;
  }

  delete[] sysNames;
  f1.close();
}

