/**
 * EMC.cc
 *
 * F2P
 * A Detailed Study of the Proton Structure Functions in Deep Inelastic Muon - Proton Scattering
 * Nucl.Phys. B259 (1985) 189, 1985
 * http://dx.doi.org/10.17182/hepdata.13916
 * From Table 237 onwards
 *
 * F2D
 * Measurements of the Nucleon Structure Functions F(2)N in Deep Inelastic
 * Muon Scattering from Deuterium and Comparison with Those from Hydrogen and Iron
 * Nucl.Phys. B293 (1987) 740, 1987
 * http://dx.doi.org/10.17182/hepdata.33821
 *
 */

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
    double mp = 0.938; //GeV
    double E;
    f1 >> E;        //E beam
    f1 >> fKin1[i]; //x
    f1 >> fKin2[i]; //Q2
    
    fKin3[i] = fKin2[i] / ( 2.*E*mp*fKin1[i] ) ; //y

    f1 >> fData[i]; //obs
 
    double stat = 0;
    double sist = 0;
    double dummy;

    f1 >> stat;
    f1 >> dummy;

    fStat[i] = stat;

    f1 >> sist;
    f1 >> dummy;

    //check here
    fSys[i][0].add = sist;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CORR_P";

    fSys[i][0].mult = fSys[i][0].add/(fData[i]*1e-2);

    fSys[i][1].add = fData[i]*0.05;
    fSys[i][1].type = MULT;
    fSys[i][1].name = "CORR_NORM";

    fSys[i][1].mult = 5;

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
    double mp = 0.938; //GeV
    double E;
    f1 >> E;        //E beam
    f1 >> fKin1[i]; //x
    f1 >> fKin2[i]; //Q2
    
    fKin3[i] = fKin2[i] / ( 2.*E*mp*fKin1[i] ) ; //y

    f1 >> fData[i]; //obs

    double stat = 0;
    double sist = 0;
    double dummy;

    f1 >> stat;
    f1 >> dummy;

    fStat[i] = stat;

    f1 >> sist;
    f1 >> dummy;

    //check here
    fSys[i][0].add = sist;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CORR_D";

    fSys[i][0].mult = fSys[i][0].add/(fData[i]*1e-2);

    fSys[i][1].add = fData[i]*0.05;
    fSys[i][1].type = MULT;
    fSys[i][1].name = "CORR_NORM";

    fSys[i][1].mult = 5;

  }

  
  f1.close();

}
