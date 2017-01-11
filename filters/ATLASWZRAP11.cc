/**
 * ATLASWZRAP11.cc
 * Implementation of ATLAS 2011 W/Z rapidity data
 */

#include "ATLAS.h"

void ATLASWZRAP11Filter::ReadData()
{

  cout << "********** WARNING: Converting pb to fb to match ApplGrid output ********" << endl;

  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << "ATLASWZRAP11/wzrap11.dat";

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
  {
    f1 >> sysNames[i];
    if (sysNames[i].compare("uncor") == 0) sysNames[i] = "UNCORR";
  }
 
  int ndataWZ[3] = {11,22,fNData};  // Data thresholds for W+, W- and Z respectively
  double MWZ2[3]= {pow(MW,2.0), pow(MW,2.0), pow(MZ,2.0)};   //Mass squared of W (+ and -) and Z

  for (int i = 0; i < fNData; i++)
  {
    const int iWZ = i < ndataWZ[0] ? 0 :(i < ndataWZ[1] ? 1:2);
    double etamin, etamax;

    // Kinematics
    f1 >> dummy; f1 >> etamin; f1 >> etamax;
    fKin1[i] = etamin + (etamax - etamin)/2.0;
    fKin2[i] = MWZ2[iWZ];
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
  
  f1.close();
}

