/**
 *  FCC.cc
 *  Pseudo-data
 *  LR
 *
 */

#include "FCC.h"
#include <random>

void FCCFilter::ReadData()
{
  // Opening files
  fstream f1;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/datfcc5060ncep.txt";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double syscor[fNData][fNSys];

  for (int i = 0; i < fNData; i++)
  {
    for (int l = 0; l < fNSys; l++)
      syscor[i][l] = 0.0;
  }

  // Filtering data
  double uncorr[fNData];

  string tmp;
  
  //skip first line
  getline(f1,tmp);
  for (int i = 0; i < fNData; i++)
  {
//#          q2         x          y        thetae  thetaj   f2    sigrNC    nevent    estat    eunco  esyst   etot       eelen    ethee  ehadr  radco egamp   effic   enois
    f1 >> fKin2[i] >> fKin1[i] >> fKin3[i] >> tmp >> tmp >> tmp >> fData[i] >> tmp >> fStat[i] >> tmp >> tmp >> uncorr[i] >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp  ;

   
    // Statistical errors - percentage with respect the observable
    fStat[i] = fStat[i]*fData[i]*1e-2;

    // Uncorrelated systematics
    fSys[i][0].mult = uncorr[i];
    fSys[i][0].type = MULT;
    fSys[i][0].name = "UNCORR";
    
    // Additive errors   
    for (int l = 0; l < fNSys; l++)
      fSys[i][l].add = fSys[i][l].mult*fData[i]*1e-2;

    //Add Random Fluctuation
    default_random_engine de(time(0)); //seed
    normal_distribution<double> nd(fData[i], sqrt(pow(fSys[i][0].add,2)+pow(fStat[i],2)) ); //mean followed by stdiv
    fData[i] = nd(de); //Generate numbers;

  }

  f1.close();
}
