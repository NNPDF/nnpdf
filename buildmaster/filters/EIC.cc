/*
   EIC pseudodata
*/

#include "EIC.h"

//Charged-current measurements, optimistic scenario
void EIC_CC_140_OPTFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/CC-optimistic.dat";
  f1.open(datafile.str().c_str(), ios::in);
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  string line;
  getline(f1,line);

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      double sqrts;
      int idum;
      lstream >> idum >> sqrts >> fKin1[i] >> fKin2[i] // x, Q2
	      >> fStat[i]
	      >> fSys[i][0].mult
	      >> fSys[i][1].mult;
      fData[i] = 1.;

      fStat[i] = fStat[i]/100. * fData[i]; 
      
      fKin3[i] = fKin2[i] / fKin1[i] / sqrts / sqrts; //y
      fSys[i][0].add  = fSys[i][0].mult / 100. * fData[i];
      fSys[i][0].type = ADD;
      fSys[i][0].name = "CORR";
      fSys[i][1].add  = fSys[i][0].mult / 100. * fData[i];
      fSys[i][1].type = MULT;
      fSys[i][1].name = "CORR";      
    }

  f1.close();
  
}

//Charged-current measurements, pessimistic scenario
void EIC_CC_140_PESFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/CC-pessimistic.dat";
  f1.open(datafile.str().c_str(), ios::in);
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  string line;
  getline(f1,line);
  
  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      double sqrts;
      int idum;
      lstream >> idum >> sqrts >> fKin1[i] >> fKin2[i] // x, Q2
	      >> fStat[i]
	      >> fSys[i][0].mult
	      >> fSys[i][1].mult;
      fData[i] = 1.;

      fStat[i] = fStat[i]/100. * fData[i]; 
      
      fKin3[i] = fKin2[i] / fKin1[i] / sqrts / sqrts; //y
      fSys[i][0].add  = fSys[i][0].mult / 100. * fData[i];
      fSys[i][0].type = ADD;
      fSys[i][0].name = "CORR";
      fSys[i][1].add  = fSys[i][0].mult / 100. * fData[i];
      fSys[i][1].type = MULT;
      fSys[i][1].name = "CORR";
    }

  f1.close();
  
}


