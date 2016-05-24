/**
 * H1 B
 * 
 * 
 */

#include "H1F2B.h"

void H1F2BFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/H1HERAF2B.data";
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

    double sist;
    f1 >> sist;
    f1 >> dummy;

    double tempsis = 0;

    //uncorrelated systematic
    f1 >> fSys[i][0].mult;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "UNCORR";

    std::string *sysNames = new std::string[fNSys];
    for (int l = 1; l < fNSys; ++l)
    {
      double temp;
      f1 >> temp; 
      fSys[i][l].mult = abs(temp);
      fSys[i][l].type = MULT;
      fSys[i][l].name = "CORR_"+sysNames[l];
    }

    for (int l = 0; l < fNSys; ++l)
    {
      tempsis += fSys[i][l].mult*fSys[i][l].mult;
    }

    cout << sqrt(tempsis)/sist << endl;

    // Additive errors   
    for (int l = 0; l < fNSys; l++)
      fSys[i][l].add = fSys[i][l].mult*fData[i]*1e-2;

  }

  
  f1.close();

}