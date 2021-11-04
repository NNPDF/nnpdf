/**
  FCC.cc
  Pseudo-data
  LR

  The original pseudo-data are taken from here
  FCC
  http://hep.ph.liv.ac.uk/~mklein/fccdata/
  http://hep.ph.liv.ac.uk/~mklein/fccdata/datfccreadme

  LHeC
  http://hep.ph.liv.ac.uk/~mklein/lhecdata/
  http://hep.ph.liv.ac.uk/~mklein/lhecdata/datlhecreadme

  where READMEs are also available.

  The files LHeC.dat and FCC.dat contain the systematics,
  but here only the total error is included without
  any breakdown of the systematics, and is taken as fully uncorrelated.

  The central values originally included in LHeC.dat and FCC.dat
  are replaced by values computed with APFEL using resummed theory
  and resummed PDFs as commented in the rawdata files.

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
  uncorr.resize(fNData);
  data.resize(fNData);
  syst.resize(fNData);

  string tmp;

  //skip first three lines
  getline(rawdata_systematics,tmp);
  getline(rawdata_systematics,tmp);
  getline(rawdata_systematics,tmp);

  //skip first line
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
      data[i] = thres;
      //For the data point whose values computed by APFEL are zero
      //because Q is too low the value is set to the value contained
      //in the original data files
      //This points will be however cut
      //as they do not pass the generic cuts
      //usually applied in any NNPDF fits
    }

    // Uncorrelated systematics - percentage with respect to the observable
    fSys[i][0].mult = uncorr[i];
    fSys[i][0].type = MULT;
    fSys[i][0].name = "UNCORR";

    // Statistical errors
    fStat[i] = syst[i]*data[i]*1e-2;

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


void FutureColliderFilterCC::ReadData()
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
  uncorr.resize(fNData);
  data.resize(fNData);
  syst.resize(fNData);

  string tmp;

  //skip first three lines
  getline(rawdata_systematics,tmp);
  getline(rawdata_systematics,tmp);
  getline(rawdata_systematics,tmp);

  //skip first line
  getline(rawdata_central,tmp);

  //random generator for fluctuations of pseudo-data
  mt19937 random_engine; //seed
  random_engine.seed(15);

  for (int i = 0; i < fNData; i++)
  {
    rawdata_systematics >> fKin2[i] >> fKin1[i] >> fKin3[i]
    >> tmp
    >> data[i]
    >> tmp
    >> syst[i]
    >> tmp >> tmp
    >> uncorr[i];

    double thres;
    rawdata_central >> tmp >> thres >> tmp;
    if(thres>0)
    {
      data[i] = thres;
      //For the data point whose values computed by APFEL are zero
      //because Q is too low the value is set to the value contained
      //in the original data files
      //This points will be however cut
      //as they do not pass the generic cuts
      //usually applied in any NNPDF fits
    }

    // Uncorrelated systematics - percentage with respect to the observable
    fSys[i][0].mult = uncorr[i];
    fSys[i][0].type = MULT;
    fSys[i][0].name = "UNCORR";

    // Statistical errors
    fStat[i] = syst[i]*data[i]*1e-2;

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
