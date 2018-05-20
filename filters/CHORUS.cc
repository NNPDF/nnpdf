/*
WARNING:
THE USE OF THIS FILTER IS DEPRECATED IN FAVOUR OF CHORUSPb.cc
*/

/**
 * CHORUS: ONENGUT et al. Phys.LETT.632(2006)65
 *
 *     This file contains the measured differential
 *     neutrino-nucleon charged-current cross-section
 *     as measured on a lead target in the CHORUS
 *     experiment at CERN.
 *
 *     The cross-section points are not corrected for
 *     the non-isoscalarity of the target or for
 *     radiative effects. To correct the data points
 *     for these effects, apply the multiplication
 *     factors given in columns "isos" and "radc"
 *
 *     Explanation of labels:
 *
 *     Enu   : central value of neutrino energy (GeV)
 *     x     : central value of Bjorken-x
 *     y     : central value of inelasticity y
 *     dsdxy : differential cross-section in 10^-38 cm^2 GeV^-1
 *     dstat : statistical uncertainty
 *     dsyst : systematic uncertainty
 *     isos  : correction factor to obtain differential
 *             cross-sections on isoscalar targets
 *     radc  : correction factor to obtain differential
 *             cross-sections corrected for QED radiation effects
 *     sh1   : hadronic energy scale 5%
 *     sh2   : hadronic energy offset 150 MeV
 *     sh3   : muon momentum scale 2.5%
 *     sh4   : muon momentum offset 150 MeV
 *     sh5   : total cross-section 2.1%
 *     sh6   : nb/nu cross-section 1.4%
 *     sh7   : non-linear total cross-section 1%/100GeV
 *     sh8   : non-linear nb/nu cross-section 0.5%/100GeV
 *
 *     sh9   : input Structure Function    (GRV98 vs GRV94)
 *     sh10  : input Radiative correction (GRV98 vs Bardin)
 *     sh11  : Reconstruction efficiency sigma          5%
 *     sh12  : Phenmenological Corrections
 *     sh13  : MC NN hadronic energy sigma             2.5%
 *
 *     Enu    x    y    dsdxy   dstat   dsyst    isos    radc    sh1     sh2     sh3     sh4     sh5     sh6     sh7     sh8     sh9     sh10    sh11    sh12     sh13
 *
 */

#include "CHORUS.h"

void CHORUSNUFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/x-sec_shift_nu.txt";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  double Mn = 0.9389;
  double enu,tmp,nortmp;
  
  // Filtering data
  string line;
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    lstream >> enu;
    lstream >> fKin1[i];   //x
    lstream >> fKin3[i];   //y
    fKin2[i] = 2.0*Mn*fKin1[i]*fKin3[i]*enu;  //q2

    lstream >> fData[i];   //obs
    lstream >> fStat[i];   //stat
    
    lstream >> tmp;
    
    //Isoscalar target correction interpreted as uncertainty
    lstream >> nortmp;
    fSys[i][0].mult = (1.0-nortmp)*100.0;
    fSys[i][0].add = fSys[i][0].mult*fData[i]*1e-2;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CHORUSISOTARGCOR";
    
    //QED radiation correction interpreted as uncertainty
    lstream >> nortmp;
    fSys[i][1].mult = (1.0-nortmp)*100.0;
    fSys[i][1].add = fSys[i][1].mult*fData[i]*1e-2;
    fSys[i][1].type = MULT;
    fSys[i][1].name = "CHORUSQEDRADCOR";
    
    //Systematics
    for (int l = 2; l < fNSys; l++)
    {
      lstream >> fSys[i][l].add;
      fSys[i][l].mult = fSys[i][l].add*100/fData[i];
      fSys[i][l].type = ADD;
      ostringstream sysname;
      sysname << "CHORUSSYS" << l-1;
      fSys[i][l].name = sysname.str();
    }
  }
  
  f1.close();
}

/**
 * CHORUS: ONENGUT et al. Phys.LETT.632(2006)65
 *
 *     This file contains the measured differential
 *     neutrino-nucleon charged-current cross-section
 *     as measured on a lead target in the CHORUS
 *     experiment at CERN.
 *
 *     The cross-section points are not corrected for
 *     the non-isoscalarity of the target or for
 *     radiative effects. To correct the data points
 *     for these effects, apply the multiplication
 *     factors given in columns "isos" and "radc"
 *
 *     Explanation of labels:
 *
 *     Enu   : central value of neutrino energy (GeV)
 *     x     : central value of Bjorken-x
 *     y     : central value of inelasticity y
 *     dsdxy : differential cross-section in 10^-38 cm^2 GeV^-1
 *     dstat : statistical uncertainty
 *     dsyst : systematic uncertainty
 *     isos  : correction factor to obtain differential
 *             cross-sections on isoscalar targets
 *     radc  : correction factor to obtain differential
 *             cross-sections corrected for QED radiation effects
 *     sh1   : hadronic energy scale 5%
 *     sh2   : hadronic energy offset 150 MeV
 *     sh3   : muon momentum scale 2.5%
 *     sh4   : muon momentum offset 150 MeV
 *     sh5   : total cross-section 2.1%
 *     sh6   : nb/nu cross-section 1.4%
 *     sh7   : non-linear total cross-section 1%/100GeV
 *     sh8   : non-linear nb/nu cross-section 0.5%/100GeV
 *
 *     sh9   : input Structure Function    (GRV98 vs GRV94)
 *     sh10  : input Radiative correction (GRV98 vs Bardin)
 *     sh11  : Reconstruction efficiency sigma          5%
 *     sh12  : Phenmenological Corrections
 *     sh13  : MC NN hadronic energy sigma             2.5%
 *
 *     Enu    x    y    dsdxy   dstat   dsyst    isos    radc    sh1     sh2     sh3     sh4     sh5     sh6     sh7     sh8     sh9     sh10    sh11    sh12     sh13
 *
 */
void CHORUSNBFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/x-sec_shift_nb.txt";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  double Mn = 0.9389;
  double enu,tmp,nortmp;
  
  // Filtering data
  string line;
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    lstream >> enu;
    lstream >> fKin1[i];   //x
    lstream >> fKin3[i];   //y
    fKin2[i] = 2.0*Mn*fKin1[i]*fKin3[i]*enu;  //q2

    lstream >> fData[i];   //obs
    lstream >> fStat[i];   //stat
    
    lstream >> tmp;
    
    //Isoscalar target correction interpreted as uncertainty
    lstream >> nortmp;
    fSys[i][0].mult = (1.0-nortmp)*100.0;
    fSys[i][0].add = fSys[i][0].mult*fData[i]*1e-2;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CHORUSISOTARGCOR";
    
    //QED radiation correction interpreted as uncertainty
    lstream >> nortmp;
    fSys[i][1].mult = (1.0-nortmp)*100.0;
    fSys[i][1].add = fSys[i][1].mult*fData[i]*1e-2;
    fSys[i][1].type = MULT;
    fSys[i][1].name = "CHORUSQEDRADCOR";
    
    //Systematics
    for (int l = 2; l < fNSys; l++)
    {
      lstream >> fSys[i][l].add;
      fSys[i][l].mult = fSys[i][l].add*100/fData[i];
      fSys[i][l].type = ADD;
      ostringstream sysname;
      sysname << "CHORUSSYS" << l-1;
      fSys[i][l].name = sysname.str();
    }
  }
  
  f1.close();
}

