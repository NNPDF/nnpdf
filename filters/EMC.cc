#include "EMC.h"

void EMCF2PFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/EMCF2P.data";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Reading data
  string line;

  //skip first two lines
  getline(f1,line);
  getline(f1,line);
    

  //Filtering data
  for (int i = 0; i < fNData; i++)
  {
    f1 >> fKin1[i]; //x
    f1 >> fKin2[i]; //Q2
    
    fKin3[i] = 0.0; //y

    f1 >> fData[i]; //obs

    double statp = 0.;
    double statm = 0.;
    double sistp = 0.;
    double sistm = 0.;

    f1 >> statp;
    f1 >> statm;

    double stat = 0.;

    stat = (abs(statp) + abs(statm))/2.;  //average (statm is negative in rawdata)
    fStat[i] = stat;

    f1 >> sistp;
    f1 >> sistm;

    double sist = 0.;

    sist = (abs(sistp) + abs(sistm))/2.;  //average (sistm is negative in rawdata)

    //check here
    fSys[i][0].add = sist;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CORR";

    fSys[i][0].mult = fSys[i][0].add/fData[i]/1e-2;

  }

  
  f1.close();

}

void EMCF2DFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/EMCF2D.data";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Reading data
  string line;

  //skip first two lines
  getline(f1,line);
  getline(f1,line);
    

  //Filtering data
  for (int i = 0; i < fNData; i++)
  {
    f1 >> fKin1[i]; //x
    f1 >> fKin2[i]; //Q2
    
    fKin3[i] = 0.0; //y

    f1 >> fData[i]; //obs

    double statp = 0.;
    double statm = 0.;
    double sistp = 0.;
    double sistm = 0.;

    f1 >> statp;
    f1 >> statm;

    double stat = 0.;

    stat = (abs(statp) + abs(statm))/2.;  //average (statm is negative in rawdata)
    fStat[i] = stat;

    f1 >> sistp;
    f1 >> sistm;

    double sist = 0.;

    sist = (abs(sistp) + abs(sistm))/2.;  //average (sistm is negative in rawdata)

    //check here
    fSys[i][0].add = sist;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CORR";

    fSys[i][0].mult = fSys[i][0].add/fData[i]/1e-2;

  }

  
  f1.close();

}