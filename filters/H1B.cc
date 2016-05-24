/**
 * H1 B
 * 
 * 
 */

#include "H1B.h"

void H1BFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/H1HERAB.data";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Reading data
  string line;

  //Filtering data
  for (int i = 0; i < fNData; i++)
  {
    double dummy;

    f1 >> fKin2[i]; //Q2
    f1 >> fKin1[i]; //x
    f1 >> fKin3[i]; //y
   
    f1 >> fData[i];
    
    f1 >> dummy;
    f1 >> fStat[i];

    double stat;
    f1 >> stat;
    f1 >> dummy;

    //uncorrelated systematic
    f1 >> fSys[i][0].mult;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "UNCORR";

    std::string *sysNames = new std::string[fNSys];
    for (int l = 1; l < fNSys; ++l)
    {
      f1 >> fSys[i][l].mult = sist;
      fSys[i][l].type = MULT;
      fSys[i][l].name = "CORR_"+sysNames[l];
    }
    
    // Additive errors   
    for (int l = 0; l < fNSys; l++)
      fSys[i][l].add = fSys[i][l].mult*fData[i]*1e-2;

  }

  
  f1.close();

}