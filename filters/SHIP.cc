/*
WARNING:
THE USE OF THIS FILTER IS DEPRECATED IN FAVOUR OF NUTECFe.cc
*/

/*
  This is the filter script for the SHIP pseudodata provided by Emircan.
  The analysis closely follows that of the NuTeV experiment.
*/


#include "SHIP.h"

//Neutrinos
void SHIPNUFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/neutrinos.res";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  double mn = 0.938;
  string line;
 
  int idum, nevts;
  double ddum, enu;

  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> idum >> idum >> idum;
    lstream >> fKin1[i];  //Bjorken x
    lstream >> fKin3[i];  //Inelasticity y
    lstream >> enu;       //Neutrion energy
    
    fKin2[i] = 2.0*mn*enu*fKin1[i]*fKin3[i]; // Q2

    lstream >> nevts;
    lstream >> fData[i];

    //Statistical uncertainty
    fStat[i] = 1/pow(nevts,0.5)*fData[i];
   
  }
  
  f1.close();
}

//Antineutrinos
void SHIPNBFilter::ReadData()
{
  //Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/antineutrinos.res";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  double mn = 0.938;
  string line;
 
  int idum, nevts;
  double ddum, enu;

  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> idum >> idum >> idum;
    lstream >> fKin1[i];  //Bjorken x
    lstream >> fKin3[i];  //Inelasticity y
    lstream >> enu;       //Neutrion energy
    
    fKin2[i] = 2.0*mn*enu*fKin1[i]*fKin3[i]; // Q2

    lstream >> nevts;
    lstream >> fData[i];

    //Statistical uncertainty
    fStat[i] = 1/pow(nevts,0.5)*fData[i];
   
  }
  
  f1.close();
}


