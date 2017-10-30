/**
 *  FCC.cc
 *  Pseudo-data
 *  LR
 *
 */

#include "FutureColliders.h"
#include <random>
#include <assert.h>

void FutureColliderFilter::ReadData()
{
  // Opening files
  fstream f1;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << "FutureColliders/" << fSetName <<".txt";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

    // Opening files
  fstream f2;

  stringstream pseudofile("");
  pseudofile << dataPath() << "rawdata/"
  << "FutureColliders/" << fSetName <<"_NNLONLLx.txt";
  f2.open(pseudofile.str().c_str(), ios::in);

  if (f2.fail()) {
    cerr << "Error opening data file " << pseudofile.str() << endl;
    exit(-1);
  }

  // Filtering data
  double uncorr[fNData], data[fNData];

  string tmp;
  
  //skip first three lines
  getline(f1,tmp);
  getline(f1,tmp);
  getline(f1,tmp);

    //skip firs line
  getline(f2,tmp);

  default_random_engine de(14); //seed

  for (int i = 0; i < fNData; i++)
  {
//#          q2         x          y        thetae  thetaj   f2    sigrNC    nevent    estat    eunco  esyst   etot       eelen    ethee  ehadr  radco egamp   effic   enois
    f1 >> fKin2[i] >> fKin1[i] >> fKin3[i]
    >> tmp >> tmp >> tmp
    >> data[i]
    >> tmp 
    >> fStat[i]
    >> tmp >> tmp
    >> uncorr[i]
    >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp ;

    double thres;
    f2 >> tmp >> thres >> tmp;
    if(thres>0)
    {
      data[i] = thres;
    }

    // Statistical errors - percentage with respect the observable
    fStat[i] = fStat[i]*data[i]*1e-2;

    // Uncorrelated systematics
    fSys[i][0].mult = uncorr[i];
    fSys[i][0].type = MULT;
    fSys[i][0].name = "UNCORR";
    
    // Additive errors   
    for (int l = 0; l < fNSys; l++)
      fSys[i][l].add = fSys[i][l].mult*data[i]*1e-2;

    //Add Random Fluctuation
    double tmp = -1;
    while( tmp <= 0 )
    {
      normal_distribution<double> nd(data[i], sqrt(pow(fSys[i][0].add,2)+pow(fStat[i],2)) ); 
      tmp = nd(de);
    }
    fData[i] = tmp;
    assert(fData[i]>0);

  }

  f1.close();
  f2.close();
}
