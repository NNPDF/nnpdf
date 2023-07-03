#include "FPF.h"

// FASERV2 EXPERIMENTS //
// Read & Dump data for FASERV2 NU INCLUSIVE
void FASERV2NUINCLUSIVEFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/FASERv2_inclusive_nu_El_fluctuated.txt";
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
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
  }
  
  f1.close();  
}

// Read & Dump data for FASERV2 NUB INCLUSIVE
void FASERV2NBINCLUSIVEFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/FASERv2_inclusive_nub_El_fluctuated.txt";
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
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
  }
  
  f1.close();  
}

// Read & Dump data for FASERV2 NU CHARM
void FASERV2NUCHARMFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/FASERv2_charm_nu_El_fluctuated.txt";
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
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
  }
  
  f1.close();  
}

// Read & Dump data for FASERV2 NUB CHARM
void FASERV2NBCHARMFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/FASERv2_charm_nub_El_fluctuated.txt";
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
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
  }
  
  f1.close();  
}
