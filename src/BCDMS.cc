
/**
 *****************************************************************************
 *     BCDMS: CERN-EP/89-06 , A. C. Benvenuti et al., Phys. Lett. B 223 (1989) 485
 *            CERN-EP/89-170, A. C. Benvenuti et al., Phys. Lett. B 237 (1990) 592
 *
 *
 *     bcd_targetbeamenergy.data (Table 3, 4, 5 and 6 of the preprint):
 *     x, y, Q2, F2, stat, dum, dum
 *     norm, fb,fs,fr, dum
 *
 *     The correlated sys error are
 *
 *     fb: calibration of the incoming muon (beam) energy;
 *     fs: calibration of the outgoing muon energy (spectrometer magnetic field);
 *     fr: spectrometer resolution;
 *
 *     they are given as percentage, while stat is an absolute value.
 *
 *     The normalization uncertainty is 3% and it is correlated for all data points
 *     There is a 2% relative normalization between different targets.
 *     There is a 1% relative normalization between different beam energies for the proton.
 *     There is a 1% relative normalization between different the first and the second
 *     beam energy for the deuteron, 1.5% between the second and the third.
 *
 *     In order to allow for loops over the total number sys and norm and to take
 *     into account correlations properly we have proceeded this way:
 *     for relative normalization between the two targets, say,
 *     we use two rel nor one of which is set to zero for the proton,
 *     and the other one is set to zero for the deuteron.
 *     Analogously for other relative normalizations.
 */

#include "BCDMS.h"

void BCDMSPFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3, f4;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/bcd_p100.data";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  stringstream datafile2("");
  datafile2 << dataPath() << "rawdata/"
  << fSetName << "/bcd_p120.data";
  f2.open(datafile2.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << datafile2.str() << endl;
    exit(-1);
  }
  
  stringstream datafile3("");
  datafile3 << dataPath() << "rawdata/"
  << fSetName << "/bcd_p200.data";
  f3.open(datafile3.str().c_str(), ios::in);
  
  if (f3.fail()) {
    cerr << "Error opening data file " << datafile3.str() << endl;
    exit(-1);
  }
  
  stringstream datafile4("");
  datafile4 << dataPath() << "rawdata/"
  << fSetName << "/bcd_p280.data";
  f4.open(datafile4.str().c_str(), ios::in);
  
  if (f4.fail()) {
    cerr << "Error opening data file " << datafile4.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  
  double datain[351][5];

  double relnorbeam = 1.0*0.5;    //relative normalisation of 1% between beam energies
  double relnortarget = 2.0*0.5;  //relative normalisation of 2% between targets
  
  // Reading data
  string line;
  double tmp;
  int nini = 0;

  for (int i = 0; i < 97; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    for (int j = 0; j < 5; j++)
      lstream >> datain[i][j];
      
    getline(f1,line);
    istringstream lstream2(line);
    lstream2 >> tmp;
    for (int j = 0; j < 3; j++)
      lstream2 >> fSys[i][j].mult;
      
    fSys[i][5].mult = relnorbeam;
    fSys[i][6].mult = relnorbeam;
    fSys[i][7].mult = relnorbeam;    
  }
  nini += 97;

  for (int i = nini; i < nini+99; i++)
  {
    getline(f2,line);
    istringstream lstream(line);
    for (int j = 0; j < 5; j++)
      lstream >> datain[i][j];
      
    getline(f2,line);
    istringstream lstream2(line);
    lstream2 >> tmp;
    for (int j = 0; j < 3; j++)
      lstream2 >> fSys[i][j].mult;
      
    fSys[i][5].mult = -relnorbeam;
    fSys[i][8].mult = relnorbeam;
    fSys[i][9].mult = relnorbeam;    
  }
  nini += 99;  

  for (int i = nini; i < nini+79; i++)
  {
    getline(f3,line);
    istringstream lstream(line);
    for (int j = 0; j < 5; j++)
      lstream >> datain[i][j];
      
    getline(f3,line);
    istringstream lstream2(line);
    lstream2 >> tmp;
    for (int j = 0; j < 3; j++)
      lstream2 >> fSys[i][j].mult;
      
    fSys[i][6].mult = -relnorbeam;
    fSys[i][8].mult = -relnorbeam;
    fSys[i][10].mult = relnorbeam;    
  } 
  nini += 79;

  for (int i = nini; i < nini+76; i++)
  {
    getline(f4,line);
    istringstream lstream(line);
    for (int j = 0; j < 5; j++)
      lstream >> datain[i][j];
      
    getline(f4,line);
    istringstream lstream2(line);
    lstream2 >> tmp;
    for (int j = 0; j < 3; j++)
      lstream2 >> fSys[i][j].mult;
      
    fSys[i][7].mult = -relnorbeam;
    fSys[i][9].mult = -relnorbeam;
    fSys[i][10].mult = -relnorbeam;    
  }  
  nini += 76;
  
  // Filtering data
  for (int i = 0; i < fNData; i++)
  {
    fKin1[i] = datain[i][0];  //x
    fKin2[i] = datain[i][2];  //q2
    fKin3[i] = datain[i][1];  //y
    
    fData[i] = datain[i][3];  //obs
    fStat[i] = datain[i][4];  //stat
    
    fSys[i][0].name = "BCDMSFB";
    fSys[i][0].type = ADD;
    fSys[i][1].name = "BCDMSFS";
    fSys[i][1].type = ADD;
    fSys[i][2].name = "BCDMSFR";
    fSys[i][2].type = ADD;
    
    fSys[i][3].mult = 3.0;   //normalisation uncertianty of 3%
    fSys[i][3].name = "BCDMSNORM";
    fSys[i][3].type = MULT;
    
    fSys[i][4].mult = relnortarget;  //relative normalisation
    fSys[i][4].name = "BCDMSRELNORMTARGET";
    fSys[i][4].type = MULT;
    
    for (int l = 5; l < fNSys; l++)
    {
      fSys[i][l].name = "CORR";
      fSys[i][l].type = MULT;
    } 
        
    //additive form of systematics
    for (int l = 0; l < fNSys; l++)
      fSys[i][l].add = fSys[i][l].mult*fData[i]*1e-2;   
  }
  
  f1.close();
  f2.close();
  f3.close();
  f4.close();
}

/**
 *****************************************************************************
 *     BCDMS: CERN-EP/89-06 , A. C. Benvenuti et al., Phys. Lett. B 223 (1989) 485
 *            CERN-EP/89-170, A. C. Benvenuti et al., Phys. Lett. B 237 (1990) 592
 *
 *
 *     bcd_targetbeamenergy.data (Table 3, 4, 5 and 6 of the preprint):
 *     x, y, Q2, F2, stat, dum, dum
 *     fb,fs,fr, dum, dum
 *
 *     The correlated sys error are
 *
 *     fb: calibration of the incoming muon (beam) energy;
 *     fs: calibration of the outgoing muon energy (spectrometer magnetic field);
 *     fr: spectrometer resolution;
 *
 *     they are given as percentage, while stat is an absolute value.
 *
 *     The normalization uncertainty is 3% and it is correlated for all data points
 *     There is a 2% relative normalization between different targets.
 *     There is a 1% relative normalization between different beam energies for the proton.
 *     There is a 1% relative normalization between different the first and the second
 *     beam energy for the deuteron, 1.5% between the second and the third.
 *
 *     In order to allow for loops over the total number sys and norm and to take
 *     into account correlations properly we have proceeded this way:
 *     for relative normalization between the two targets, say,
 *     we use two rel nor one of which is set to zero for the proton,
 *     and the other one is set to zero for the deuteron.
 *     Analogously for other relative normalizations.
 */
void BCDMSDFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/bcd_d120.data";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  stringstream datafile2("");
  datafile2 << dataPath() << "rawdata/"
  << fSetName << "/bcd_d200.data";
  f2.open(datafile2.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << datafile2.str() << endl;
    exit(-1);
  }
  
  stringstream datafile3("");
  datafile3 << dataPath() << "rawdata/"
  << fSetName << "/bcd_d280.data";
  f3.open(datafile3.str().c_str(), ios::in);
  
  if (f3.fail()) {
    cerr << "Error opening data file " << datafile3.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  double datain[fNData][12];
  
  double relnorbeam1 = 1.0*0.5;            //relative normalisation of 1% between first and second bins
  double relnorbeam2 = 1.5*0.5;            //relative normalisation of 1.5% between second and third bins
  double relnorbeam3 = sqrt(0.5*relnorbeam1*relnorbeam1+0.5*relnorbeam2*relnorbeam2);   //take quadratic mean for relnor between first and third bin
  double relnortarget = 2.0*0.5;           //relative normalisation of 2% between targets
  
  // Reading data
  string line;
  double tmp;
  int nini = 0;
  
  for (int i = 0; i < 99; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    for (int j = 0; j < 5; j++)
      lstream >> datain[i][j];
      
    getline(f1,line);
    istringstream lstream2(line);
    lstream2 >> tmp;
    for (int j = 0; j < 3; j++)
      lstream2 >> fSys[i][j].mult;

    fSys[i][5].mult = relnorbeam1;
    fSys[i][7].mult = relnorbeam3;
  }
  nini += 99;

  for (int i = nini; i < nini+79; i++)
  {
    getline(f2,line);
    istringstream lstream(line);
    for (int j = 0; j < 5; j++)
      lstream >> datain[i][j];
      
    getline(f2,line);
    istringstream lstream2(line);
    lstream2 >> tmp;
    for (int j = 0; j < 3; j++)
      lstream2 >> fSys[i][j].mult;

    fSys[i][5].mult = -relnorbeam1;
    fSys[i][6].mult = relnorbeam2;
  }
  nini += 79;

  for (int i = nini; i < nini+76; i++)
  {
    getline(f3,line);
    istringstream lstream(line);
    for (int j = 0; j < 5; j++)
      lstream >> datain[i][j];
      
    getline(f3,line);
    istringstream lstream2(line);
    lstream2 >> tmp;
    for (int j = 0; j < 3; j++)
      lstream2 >> fSys[i][j].mult;

    fSys[i][6].mult = -relnorbeam2;
    fSys[i][7].mult = -relnorbeam3;
  }
  nini += 76;
  
    // Filtering data
  for (int i = 0; i < fNData; i++)
  {
    fKin1[i] = datain[i][0];  //x
    fKin2[i] = datain[i][2];  //q2
    fKin3[i] = datain[i][1];  //y
    
    fData[i] = datain[i][3];  //obs
    fStat[i] = datain[i][4];  //stat

    fSys[i][0].name = "BCDMSFB";
    fSys[i][0].type = ADD;
    fSys[i][1].name = "BCDMSFS";
    fSys[i][1].type = ADD;
    fSys[i][2].name = "BCDMSFR";
    fSys[i][2].type = ADD;
    
    fSys[i][3].mult = 3.0;   //normalisation uncertianty of 3%
    fSys[i][3].name = "BCDMSNORM";
    fSys[i][3].type = MULT;
    
    fSys[i][4].mult = -relnortarget;  //relative normalisation
    fSys[i][4].name = "BCDMSRELNORMTARGET";
    fSys[i][4].type = MULT;
    
    for (int l = 5; l < fNSys; l++)
    {
      fSys[i][l].name = "CORR";
      fSys[i][l].type = MULT;
    } 
    
    //additive form of systematics
    for (int l = 0; l < fNSys; l++)
      fSys[i][l].add = fSys[i][l].mult*fData[i]*1e-2;   
  }
  
  f1.close();
  f2.close();
  f3.close();
}

