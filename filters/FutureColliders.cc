/**
 *  FCC.cc
 *  Pseudo-data
 *  LR
 *
 */

#include "FutureColliders.h"
#include <random>
#include <assert.h>
#include <vector>

void FutureColliderFilter::ReadData()
{
  // Opening files
  fstream rawdata_systematics;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << "FutureColliders/" << fSetName <<".txt";
  rawdata_systematics.open(datafile.str().c_str(), ios::in);

  if (rawdata_systematics.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

    // Opening files
  fstream rawdata_central;

  stringstream pseudofile("");
  pseudofile << dataPath() << "rawdata/"
  << "FutureColliders/" << fSetName <<"_NNLONLLx.txt";
  rawdata_central.open(pseudofile.str().c_str(), ios::in);

  if (rawdata_central.fail()) {
    cerr << "Error opening data file " << pseudofile.str() << endl;
    exit(-1);
  }

  // Filtering data
  vector<double> uncorr, data, syst;
  uncorr.reserve(fNData); 
  data.reserve(fNData);
  syst.reserve(fNData); 

  string tmp;
  
  //skip first three lines
  getline(rawdata_systematics,tmp);
  getline(rawdata_systematics,tmp);
  getline(rawdata_systematics,tmp);

  //skip firs line
  getline(rawdata_central,tmp);

  //random generator for fluctuations of pseudo-data
  mt19937 random_engine; //seed
  random_engine.seed(14);

  for (int i = 0; i < fNData; i++)
  {
    rawdata_systematics >> fKin2[i] >> fKin1[i] >> fKin3[i]
    >> tmp >> tmp >> tmp
    >> data[i]
    >> tmp 
    >> syst[i]
    >> tmp >> tmp
    >> uncorr[i]
    >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp ;

    double thres;
    rawdata_central >> tmp >> thres >> tmp;
    if(thres>0)
    {
      data[i] = thres;  //this ensures that the for values of Q too little (which will be cut eventually) the central value is greater than zero
    }

    // Statistical errors - percentage with respect the observable
    fStat[i] = syst[i]*data[i]*1e-2;

    // Uncorrelated systematics
    fSys[i][0].mult = uncorr[i];
    fSys[i][0].type = MULT;
    fSys[i][0].name = "UNCORR";
    
    // Additive errors   
    for (int l = 0; l < fNSys; l++)
      fSys[i][l].add = fSys[i][l].mult*data[i]*1e-2;

    //Add Random Fluctuation
    double tmp = -1;
    while( tmp <= 0 ) //ensure that the central value is bigger than zero
    {
      normal_distribution<double> nd(data[i], sqrt(pow(fSys[i][0].add,2)+pow(fStat[i],2)) ); 
      tmp = nd(random_engine);
    }
    fData[i] = tmp;
    assert(fData[i]>0);  //further check (overkill)

  }

  rawdata_systematics.close();
  rawdata_central.close();
}
