#include "FPF.h"

// Define the commondata from the Systematic Error types
// Possible values: El, Eh, Theta, comb
string ERRTYPE = "comb";


// #########################################################
// AdvSND EXPERIMENTS //
// Read & Dump data for AdvSND NU INCLUSIVE
void AdvSNDNUINCLUSIVEFilter::ReadData()
{
  // Opening files
  fstream f1;
  string filename = "/AdvSND_inclusive_nu_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  for (int i = 0; i < fNData; i++)
  {
    float _throw;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
  }
  
  f1.close();  
}

// Read & Dump data for AdvSND NUB INCLUSIVE
void AdvSNDNBINCLUSIVEFilter::ReadData()
{
  // Opening files
  fstream f1;
  string filename = "/AdvSND_inclusive_nub_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  for (int i = 0; i < fNData; i++)
  {
    float _throw;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
  }
  
  f1.close();  
}

// Read & Dump data for AdvSND SUM INCLUSIVE
void AdvSNDSUMINCLUSIVEFilter::ReadData()
{
  // Opening files
  fstream f1;
  string filename = "/AdvSND_inclusive_nochargediscrimination_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  for (int i = 0; i < fNData; i++)
  {
    float _throw;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
  }
  
  f1.close();  
}

// Read & Dump data for AdvSND NU CHARM
void AdvSNDNUCHARMFilter::ReadData()
{
  // Opening files
  fstream f1;
  string filename = "/AdvSND_charm_nu_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  for (int i = 0; i < fNData; i++)
  {
    float _throw;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
  }
  
  f1.close();  
}

// // Read & Dump data for AdvSND NUB CHARM
// void AdvSNDNBCHARMFilter::ReadData()
// {
//   // Opening files
//   fstream f1;
//   string filename = "/AdvSND_charm_nub_" + ERRTYPE + "_fluctuated.txt";
//   
//   stringstream datafile("");
//   datafile << dataPath() << "rawdata/" << fSetName << filename;
//   f1.open(datafile.str().c_str(), ios::in);
//   
//   if (f1.fail()) {
//     cerr << "Error opening data file " << datafile.str() << endl;
//     exit(-1);
//   }
//   
//   // Starting filter
//   string line;
//   for (int i = 0; i < fNData; i++)
//   {
//     getline(f1,line);
//     istringstream lstream(line);
//     
//     lstream >> fKin1[i];   // x
//     lstream >> fKin3[i];   // y
//     lstream >> fKin2[i];   // Q2 
//     lstream >> fData[i];   // central values
//
//     // Extract uncertainty values
//     lstream >> fStat[i];   // statistical uncertainties
//     lstream >> fSys[i][0].add; // systematic uncertainties
//     fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
//     fSys[i][0].type = ADD;
//     fSys[i][0].name = "UNCORR";
//   }
//   
//   f1.close();  
// }

// Read & Dump data for AdvSND SUM CHARM
void AdvSNDSUMCHARMFilter::ReadData()
{
  // Opening files
  fstream f1;
  string filename = "/AdvSND_charm_nochargediscrimination_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  for (int i = 0; i < fNData; i++)
  {
    float _throw;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
  }
  
  f1.close();  
}
// #########################################################
// FLArE10 EXPERIMENTS //
// Read & Dump data for FLArE10 NU INCLUSIVE
void FLArE10NUINCLUSIVEFilter::ReadData()
{
  // Opening files
  fstream f1;
  string filename = "/FLArE10_inclusive_nu_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  for (int i = 0; i < fNData; i++)
  {
    float _throw;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
  }
  
  f1.close();  
}

// Read & Dump data for FLArE10 NUB INCLUSIVE
void FLArE10NBINCLUSIVEFilter::ReadData()
{
  // Opening files
  fstream f1;
  string filename = "/FLArE10_inclusive_nub_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  for (int i = 0; i < fNData; i++)
  {
    float _throw;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
  }
  
  f1.close();  
}

// Read & Dump data for FLArE10 SUM INCLUSIVE
void FLArE10SUMINCLUSIVEFilter::ReadData()
{
  // Opening files
  fstream f1;
  string filename = "/FLArE10_inclusive_nochargediscrimination_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  for (int i = 0; i < fNData; i++)
  {
    float _throw;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
  }
  
  f1.close();  
}

// Read & Dump data for FLArE10 NU CHARM
void FLArE10NUCHARMFilter::ReadData()
{
  // Opening files
  fstream f1;
  string filename = "/FLArE10_charm_nu_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  for (int i = 0; i < fNData; i++)
  {
    float _throw;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
  }
  
  f1.close();  
}

// Read & Dump data for FLArE10 NUB CHARM
void FLArE10NBCHARMFilter::ReadData()
{
  // Opening files
  fstream f1;
  string filename = "/FLArE10_charm_nub_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  for (int i = 0; i < fNData; i++)
  {
    float _throw;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
  }
  
  f1.close();  
}

// Read & Dump data for FLArE10 SUM CHARM
void FLArE10SUMCHARMFilter::ReadData()
{
  // Opening files
  fstream f1;
  string filename = "/FLArE10_charm_nochargediscrimination_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  for (int i = 0; i < fNData; i++)
  {
    float _throw;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
  }
  
  f1.close();  
}

// #########################################################
// FASERV2 EXPERIMENTS //
// Read & Dump data for FASERV2 NU INCLUSIVE
void FASERV2NUINCLUSIVEFilter::ReadData()
{
  // Opening files
  fstream f1;
  string filename = "/FASERv2_inclusive_nu_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  for (int i = 0; i < fNData; i++)
  {
    float _throw;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
  }
  
  f1.close();  
}

// Read & Dump data for FASERV2 NUB INCLUSIVE
void FASERV2NBINCLUSIVEFilter::ReadData()
{
  // Opening files
  fstream f1;
  string filename = "/FASERv2_inclusive_nub_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  for (int i = 0; i < fNData; i++)
  {
    float _throw;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
  }
  
  f1.close();  
}

// Read & Dump data for FASERV2 SUM INCLUSIVE
void FASERV2SUMINCLUSIVEFilter::ReadData()
{
  // Opening files
  fstream f1;
  string filename = "/FASERv2_inclusive_nochargediscrimination_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  for (int i = 0; i < fNData; i++)
  {
    float _throw;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
  }
  
  f1.close();  
}

// Read & Dump data for FASERV2 NU CHARM
void FASERV2NUCHARMFilter::ReadData()
{
  // Opening files
  fstream f1;
  string filename = "/FASERv2_charm_nu_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  for (int i = 0; i < fNData; i++)
  {
    float _throw;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
  }
  
  f1.close();  
}

// Read & Dump data for FASERV2 NUB CHARM
void FASERV2NBCHARMFilter::ReadData()
{
  // Opening files
  fstream f1;
  string filename = "/FASERv2_charm_nub_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  for (int i = 0; i < fNData; i++)
  {
    float _throw;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
  }
  
  f1.close();  
}

// Read & Dump data for FASERV2 SUM CHARM
void FASERV2SUMCHARMFilter::ReadData()
{
  // Opening files
  fstream f1;
  string filename = "/FASERv2_charm_nochargediscrimination_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  for (int i = 0; i < fNData; i++)
  {
    float _throw;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
  }
  
  f1.close();  
}


// #########################################################
// FASERV EXPERIMENTS //
// Read & Dump data for FASERV NU INCLUSIVE
void FASERVNUINCLUSIVEFilter::ReadData()
{
  // Opening files
  fstream f1;
  string filename = "/FASERv_inclusive_nu_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  for (int i = 0; i < fNData; i++)
  {
    float _throw;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
  }
  
  f1.close();  
}

// Read & Dump data for FASERV NUB INCLUSIVE
void FASERVNBINCLUSIVEFilter::ReadData()
{
    float _throw;
  // Opening files
  fstream f1;
  string filename = "/FASERv_inclusive_nub_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
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
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
  }
  
  f1.close();  
}

// Read & Dump data for FASERV SUM INCLUSIVE
void FASERVSUMINCLUSIVEFilter::ReadData()
{
  // Opening files
  fstream f1;
  string filename = "/FASERv_inclusive_nochargediscrimination_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  for (int i = 0; i < fNData; i++)
  {
    float _throw;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
  }
  
  f1.close();  
}

// // Read & Dump data for FASERV NU CHARM
// void FASERVNUCHARMFilter::ReadData()
// {
//   // Opening files
//   fstream f1;
//   string filename = "/FASERv_charm_nu_" + ERRTYPE + "_fluctuated.txt";
//   
//   stringstream datafile("");
//   datafile << dataPath() << "rawdata/" << fSetName << filename;
//   f1.open(datafile.str().c_str(), ios::in);
//   
//   if (f1.fail()) {
//     cerr << "Error opening data file " << datafile.str() << endl;
//     exit(-1);
//   }
//   
//   // Starting filter
//   string line;
//   for (int i = 0; i < fNData; i++)
//   {
//     getline(f1,line);
//     istringstream lstream(line);
//     
//     lstream >> fKin1[i];   // x
//     lstream >> fKin3[i];   // y
//     lstream >> fKin2[i];   // Q2 
//     lstream >> fData[i];   // central values
//
//     // Extract uncertainty values
//     lstream >> fStat[i];   // statistical uncertainties
//     lstream >> fSys[i][0].add; // systematic uncertainties
//     fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
//     fSys[i][0].type = ADD;
//     fSys[i][0].name = "UNCORR";
//   }
//   
//   f1.close();  
// }
//
// // Read & Dump data for FASERV NUB CHARM
// void FASERVNBCHARMFilter::ReadData()
// {
//   // Opening files
//   fstream f1;
//   string filename = "/FASERv_charm_nub_" + ERRTYPE + "_fluctuated.txt";
//   
//   stringstream datafile("");
//   datafile << dataPath() << "rawdata/" << fSetName << filename;
//   f1.open(datafile.str().c_str(), ios::in);
//   
//   if (f1.fail()) {
//     cerr << "Error opening data file " << datafile.str() << endl;
//     exit(-1);
//   }
//   
//   // Starting filter
//   string line;
//   for (int i = 0; i < fNData; i++)
//   {
//     getline(f1,line);
//     istringstream lstream(line);
//     
//     lstream >> fKin1[i];   // x
//     lstream >> fKin3[i];   // y
//     lstream >> fKin2[i];   // Q2 
//     lstream >> fData[i];   // central values
//
//     // Extract uncertainty values
//     lstream >> fStat[i];   // statistical uncertainties
//     lstream >> fSys[i][0].add; // systematic uncertainties
//     fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
//     fSys[i][0].type = ADD;
//     fSys[i][0].name = "UNCORR";
//   }
//   
//   f1.close();  
// }

// Read & Dump data for FASERV SUM CHARM
void FASERVSUMCHARMFilter::ReadData()
{
  // Opening files
  fstream f1;
  string filename = "/FASERv_charm_nochargediscrimination_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  for (int i = 0; i < fNData; i++)
  {
    float _throw;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
  }
  
  f1.close();  
}

// #########################################################
// FLArE100 EXPERIMENTS //
// Read & Dump data for FLArE100 NU INCLUSIVE
void FLArE100NUINCLUSIVEFilter::ReadData()
{
  // Opening files
  fstream f1;
  string filename = "/FLArE100_inclusive_nu_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  for (int i = 0; i < fNData; i++)
  {
    float _throw;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
  }
  
  f1.close();  
}

// Read & Dump data for FLArE100 NUB INCLUSIVE
void FLArE100NBINCLUSIVEFilter::ReadData()
{
  // Opening files
  fstream f1;
  string filename = "/FLArE100_inclusive_nub_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  for (int i = 0; i < fNData; i++)
  {
    float _throw;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
  }
  
  f1.close();  
}

// Read & Dump data for FLArE100 SUM INCLUSIVE
void FLArE100SUMINCLUSIVEFilter::ReadData()
{
  // Opening files
  fstream f1;
  string filename = "/FLArE100_inclusive_nochargediscrimination_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  for (int i = 0; i < fNData; i++)
  {
    float _throw;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
  }
  
  f1.close();  
}

// Read & Dump data for FLArE100 NU CHARM
void FLArE100NUCHARMFilter::ReadData()
{
  // Opening files
  fstream f1;
  string filename = "/FLArE100_charm_nu_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  for (int i = 0; i < fNData; i++)
  {
    float _throw;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
  }
  
  f1.close();  
}

// Read & Dump data for FLArE100 NUB CHARM
void FLArE100NBCHARMFilter::ReadData()
{
  // Opening files
  fstream f1;
  string filename = "/FLArE100_charm_nub_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  for (int i = 0; i < fNData; i++)
  {
    float _throw;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
  }
  
  f1.close();  
}

// Read & Dump data for FLArE100 SUM CHARM
void FLArE100SUMCHARMFilter::ReadData()
{
  // Opening files
  fstream f1;
  string filename = "/FLArE100_charm_nochargediscrimination_" + ERRTYPE + "_fluctuated.txt";
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  for (int i = 0; i < fNData; i++)
  {
    float _throw;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin1[i];   // x
    lstream >> fKin3[i];   // y
    lstream >> fKin2[i];   // Q2 
    lstream >> fData[i];   // central values

    // Extract uncertainty values
    // lstream >> fStat[i];   // statistical uncertainties
    lstream >> fSys[i][0].add; // systematic uncertainties
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    fStat[i] = 0.; // implement the statistical uncertainties as if they were uncorrelated sys
    lstream >> _throw;
  }
  
  f1.close();  
}
