
#include "demo.h"

/**
 *	Comments on Experiment
 */
void SETNAMEFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << 
  "/DATAFILE.dat";


  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  string line;
  for (int i = 0; i < fNData; i++)
  {
  	// Example where the datafile is layed out as:
  	// Kin1 Kin2 Data StatError 
    getline(f1,line);
    istringstream lstream(line);
    lstream >> fKin1[i]; 
    lstream >> fKin2[i];
    lstream >> fData[i];
    
    lstream >> fStat[i];

  }
  
  f1.close();
}