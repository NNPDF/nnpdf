#include "FPF.h"

void FASERVBARFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/fluctuated_data.dat";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  
  
  string line;
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x_avg
    lstream >> fKin2[i];   // Q2_avg
    lstream >> fKin3[i];   // Enu_avg 
    lstream >> fData[i];
    lstream >> fStat[i];     // stat
    
  }
  
  f1.close();  
}

void FASERVFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/fluctuated_data.dat";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  
  
  string line;
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x_avg
    lstream >> fKin2[i];   // Q2_avg
    lstream >> fKin3[i];   // Enu_avg 
    lstream >> fData[i];
    lstream >> fStat[i];     // stat
    
  }
  
  f1.close();  
}


void FASERVBAR2Filter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/fluctuated_data.dat";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  
  
  string line;
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x_avg
    lstream >> fKin2[i];   // Q2_avg
    lstream >> fKin3[i];   // Enu_avg 
    lstream >> fData[i];
    lstream >> fStat[i];     // stat
    
  }
  
  f1.close();  
}

void FASERV2Filter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/fluctuated_data.dat";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  
  
  string line;
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x_avg
    lstream >> fKin2[i];   // Q2_avg
    lstream >> fKin3[i];   // Enu_avg 
    lstream >> fData[i];
    lstream >> fStat[i];     // stat
    
  }
  
  f1.close();  
}