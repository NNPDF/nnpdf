/*
WARNING: File modified by ERN Nov 2020.
Additional data sets, with suffix _dw and _sh have been added with extra 
systematic ucnertainties. These systematic ucnertainties account for nuclear 
uncertainties (estimated according to 1812.09074).
The two strategies (dw=deweighted and sh=shifted) are implemented.
The necessary shifts can be printed on screen and should be pasted into the
appropriate cfactor file.
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

#include "CHORUSPb.h"

void CHORUSNUPbFilter::ReadData()
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
    lstream >> nortmp;   

    //QED radiation correction interpreted as uncertainty
    lstream >> nortmp;
    fSys[i][0].mult = (1.0-nortmp)*100.0;
    fSys[i][0].add = fSys[i][0].mult*fData[i]*1e-2;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CHORUSQEDRADCOR";
    
    //Systematics
    for (int l = 1; l < fNSys; l++)
    {
      lstream >> fSys[i][l].add;
      fSys[i][l].mult = fSys[i][l].add*100/fData[i];
      fSys[i][l].type = ADD;
      ostringstream sysname;
      sysname << "CHORUSSYS" << l;
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
void CHORUSNBPbFilter::ReadData()
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
    lstream >> nortmp;    

    //QED radiation correction interpreted as uncertainty
    lstream >> nortmp;
    fSys[i][0].mult = (1.0-nortmp)*100.0;
    fSys[i][0].add = fSys[i][0].mult*fData[i]*1e-2;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CHORUSQEDRADCOR";
    
    //Systematics
    for (int l = 1; l < fNSys; l++)
    {
      lstream >> fSys[i][l].add;
      fSys[i][l].mult = fSys[i][l].add*100/fData[i];
      fSys[i][l].type = ADD;
      ostringstream sysname;
      sysname << "CHORUSSYS" << l;
      fSys[i][l].name = sysname.str();
    }
  }
  
  f1.close();
}

void CHORUSNUPb_dwFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/CHORUSNUPb/x-sec_shift_nu.txt";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream nuclearfile("");
  nuclearfile << dataPath() << "rawdata/CHORUSNUPb/nuclear/output/tables/group_result_table.csv";
  f2.open(nuclearfile.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << nuclearfile.str() << endl;
    exit(-1);
  }
  
  stringstream protonfile("");
  protonfile << dataPath() << "rawdata/CHORUSNUPb/proton/output/tables/group_result_table.csv";
  f3.open(protonfile.str().c_str(), ios::in);
  
  if (f3.fail()) {
    cerr << "Error opening data file " << protonfile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double Mn = 0.9389;
  double enu,tmp,nortmp;
  int nrep=1000;
  int nrealsys=14;
  
  // Filtering data
  string line;

  getline(f2,line);
  getline(f3,line);
  
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
    lstream >> nortmp;   

    //QED radiation correction interpreted as uncertainty
    lstream >> nortmp;
    fSys[i][0].mult = (1.0-nortmp)*100.0;
    fSys[i][0].add = fSys[i][0].mult*fData[i]*1e-2;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CHORUSQEDRADCOR1";
    
    //Systematics
    for (int l = 1; l < nrealsys; l++)
    {
      lstream >> fSys[i][l].add;
      fSys[i][l].mult = fSys[i][l].add*100/fData[i];
      fSys[i][l].type = ADD;
      ostringstream sysname;
      sysname << "CHORUSSYS_1" << l;
      fSys[i][l].name = sysname.str();
    }

    //Get proton central value
    getline(f3,line);
    istringstream pstream(line);
    string sdum;
    int idum;
    double ddum;
    double proton_cv;
    pstream >> sdum >> sdum >> idum >> ddum >> proton_cv;
    
    //Get nuclear replicas
    getline(f2,line);
    istringstream nstream(line);
    nstream >> sdum >> sdum >> idum >> ddum >> ddum;
    
    vector<double> nuclear_cv (nrep);
      
    for(int irep=0; irep<nrep; irep++)
      {
	nstream >> nuclear_cv[irep];
      }

    //Compute additional uncertainties
    for(int l=nrealsys; l<fNSys; l++)
      {
	fSys[i][l].add = (nuclear_cv[l-nrealsys] - proton_cv)/sqrt(nrep);
	fSys[i][l].mult = fSys[i][l].add*100/fData[i];
	fSys[i][l].type = ADD;
	ostringstream sysname;
	sysname << "NUCLEAR" << l-nrealsys;
	fSys[i][l].name = sysname.str();
      }

  }
  
  f1.close();
  f2.close();
  f3.close();
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
void CHORUSNBPb_dwFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/CHORUSNBPb/x-sec_shift_nb.txt";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream nuclearfile("");
  nuclearfile << dataPath() << "rawdata/CHORUSNBPb/nuclear/output/tables/group_result_table.csv";
  f2.open(nuclearfile.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << nuclearfile.str() << endl;
    exit(-1);
  }
  
  stringstream protonfile("");
  protonfile << dataPath() << "rawdata/CHORUSNBPb/proton/output/tables/group_result_table.csv";
  f3.open(protonfile.str().c_str(), ios::in);
  
  if (f3.fail()) {
    cerr << "Error opening data file " << protonfile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double Mn = 0.9389;
  double enu,tmp,nortmp;
  int nrep=1000;
  int nrealsys=14;
  
  // Filtering data
  string line;

  getline(f2,line);
  getline(f3,line);
  
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
    lstream >> nortmp;    

    //QED radiation correction interpreted as uncertainty
    lstream >> nortmp;
    fSys[i][0].mult = (1.0-nortmp)*100.0;
    fSys[i][0].add = fSys[i][0].mult*fData[i]*1e-2;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CHORUSQEDRADCOR1";
    
    //Systematics
    for (int l = 1; l < nrealsys; l++)
    {
      lstream >> fSys[i][l].add;
      fSys[i][l].mult = fSys[i][l].add*100/fData[i];
      fSys[i][l].type = ADD;
      ostringstream sysname;
      sysname << "CHORUSSYS_1" << l;
      fSys[i][l].name = sysname.str();
    }

    //Get proton central value
    getline(f3,line);
    istringstream pstream(line);
    string sdum;
    int idum;
    double ddum;
    double proton_cv;
    pstream >> sdum >> sdum >> idum >> ddum >> proton_cv;

    //Get nuclear replicas
    getline(f2,line);
    istringstream nstream(line);
    nstream >> sdum >> sdum >> idum >> ddum >> ddum;
    vector<double> nuclear_cv (nrep);
      
    for(int irep=0; irep<nrep; irep++)
      {
	nstream >> nuclear_cv[irep];
      }

    //Compute additional uncertainties
    for(int l=nrealsys; l<fNSys; l++)
      {
	fSys[i][l].add = (nuclear_cv[l-nrealsys] - proton_cv)/sqrt(nrep);
	fSys[i][l].mult = fSys[i][l].add*100/fData[i];
	fSys[i][l].type = ADD;
	ostringstream sysname;
	sysname << "NUCLEAR" << l-nrealsys;
	fSys[i][l].name = sysname.str();
      }
  }
  
  f1.close();
  f2.close();
  f3.close();
  
}

void CHORUSNUPb_shFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/CHORUSNUPb/x-sec_shift_nu.txt";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream nuclearfile("");
  nuclearfile << dataPath() << "rawdata/CHORUSNUPb/nuclear/output/tables/group_result_table.csv";
  f2.open(nuclearfile.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << nuclearfile.str() << endl;
    exit(-1);
  }
  
  stringstream protonfile("");
  protonfile << dataPath() << "rawdata/CHORUSNUPb/proton/output/tables/group_result_table.csv";
  f3.open(protonfile.str().c_str(), ios::in);
  
  if (f3.fail()) {
    cerr << "Error opening data file " << protonfile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double Mn = 0.9389;
  double enu,tmp,nortmp;
  int nrep=1000;
  int nrealsys=14;
  
  // Filtering data
  string line;

  getline(f2,line);
  getline(f3,line);
  
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
    lstream >> nortmp;   

    //QED radiation correction interpreted as uncertainty
    lstream >> nortmp;
    fSys[i][0].mult = (1.0-nortmp)*100.0;
    fSys[i][0].add = fSys[i][0].mult*fData[i]*1e-2;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CHORUSQEDRADCOR1";
    
    //Systematics
    for (int l = 1; l < nrealsys; l++)
    {
      lstream >> fSys[i][l].add;
      fSys[i][l].mult = fSys[i][l].add*100/fData[i];
      fSys[i][l].type = ADD;
      ostringstream sysname;
      sysname << "CHORUSSYS_1" << l;
      fSys[i][l].name = sysname.str();
    }

    //Get proton central value
    getline(f3,line);
    istringstream pstream(line);
    string sdum;
    int idum;
    double ddum;
    double proton_cv;
    pstream >> sdum >> sdum >> idum >> ddum >> proton_cv;

    //Get nuclear replicas
    getline(f2,line);
    istringstream nstream(line);
    double nuclear;
    nstream >> sdum >> sdum >> idum >> ddum >> nuclear;
    vector<double> nuclear_cv (nrep);
      
    for(int irep=0; irep<nrep; irep++)
      {
	nstream >> nuclear_cv[irep];
      }

    //Compute additional uncertainties
    for(int l=nrealsys; l<fNSys; l++)
      {
	fSys[i][l].add = (nuclear_cv[l-nrealsys] - nuclear)/sqrt(nrep);
	fSys[i][l].mult = fSys[i][l].add*100/fData[i];
	fSys[i][l].type = ADD;
	ostringstream sysname;
	sysname << "NUCLEAR" << l-nrealsys;
	fSys[i][l].name = sysname.str();
      }

    //Compute shifts
    //cout << nuclear/proton_cv << "   " << 0.0 << endl;
    
  }
  
  f1.close();
  f2.close();
  f3.close();
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
void CHORUSNBPb_shFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/CHORUSNBPb/x-sec_shift_nb.txt";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream nuclearfile("");
  nuclearfile << dataPath() << "rawdata/CHORUSNBPb/nuclear/output/tables/group_result_table.csv";
  f2.open(nuclearfile.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << nuclearfile.str() << endl;
    exit(-1);
  }
  
  stringstream protonfile("");
  protonfile << dataPath() << "rawdata/CHORUSNBPb/proton/output/tables/group_result_table.csv";
  f3.open(protonfile.str().c_str(), ios::in);
  
  if (f3.fail()) {
    cerr << "Error opening data file " << protonfile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double Mn = 0.9389;
  double enu,tmp,nortmp;
  int nrep=1000;
  int nrealsys=14;
  
  // Filtering data
  string line;

  getline(f2,line);
  getline(f3,line);
  
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
    lstream >> nortmp;    

    //QED radiation correction interpreted as uncertainty
    lstream >> nortmp;
    fSys[i][0].mult = (1.0-nortmp)*100.0;
    fSys[i][0].add = fSys[i][0].mult*fData[i]*1e-2;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CHORUSQEDRADCOR1";
    
    //Systematics
    for (int l = 1; l < nrealsys; l++)
    {
      lstream >> fSys[i][l].add;
      fSys[i][l].mult = fSys[i][l].add*100/fData[i];
      fSys[i][l].type = ADD;
      ostringstream sysname;
      sysname << "CHORUSSYS_1" << l;
      fSys[i][l].name = sysname.str();
    }

    //Get proton central value
    getline(f3,line);
    istringstream pstream(line);
    string sdum;
    int idum;
    double ddum;
    double proton_cv;
    pstream >> sdum >> sdum >> idum >> ddum >> proton_cv;

    //Get nuclear replicas
    getline(f2,line);
    istringstream nstream(line);
    double nuclear;
    nstream >> sdum >> sdum >> idum >> ddum >> nuclear;
    vector<double> nuclear_cv (nrep);
      
    for(int irep=0; irep<nrep; irep++)
      {
	nstream >> nuclear_cv[irep];
      }

    //Compute additional uncertainties
    for(int l=nrealsys; l<fNSys; l++)
      {
	fSys[i][l].add = (nuclear_cv[l-nrealsys] - nuclear)/sqrt(nrep);
	fSys[i][l].mult = fSys[i][l].add*100/fData[i];
	fSys[i][l].type = ADD;
	ostringstream sysname;
	sysname << "NUCLEAR" << l-nrealsys;
	fSys[i][l].name = sysname.str();
      }
    
    //Compute shifts
    //cout << nuclear/proton_cv << "   " << 0.0 << endl;
    
  }
  
  f1.close();
  f2.close();
  f3.close();
  
}

void CHORUSNUPb_dw_iteFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/CHORUSNUPb/x-sec_shift_nu.txt";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream nuclearfile("");
  nuclearfile << dataPath() << "rawdata/CHORUSNUPb/nuclear/output/tables/group_result_table.csv";
  f2.open(nuclearfile.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << nuclearfile.str() << endl;
    exit(-1);
  }
  
  stringstream protonfile("");
  protonfile << dataPath() << "rawdata/CHORUSNUPb/proton_ite/output/tables/group_result_table.csv";
  f3.open(protonfile.str().c_str(), ios::in);
  
  if (f3.fail()) {
    cerr << "Error opening data file " << protonfile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double Mn = 0.9389;
  double enu,tmp,nortmp;
  int nrep=1000;
  int nrealsys=14;
  
  // Filtering data
  string line;

  getline(f2,line);
  getline(f3,line);
  
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
    lstream >> nortmp;   

    //QED radiation correction interpreted as uncertainty
    lstream >> nortmp;
    fSys[i][0].mult = (1.0-nortmp)*100.0;
    fSys[i][0].add = fSys[i][0].mult*fData[i]*1e-2;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CHORUSQEDRADCOR1";
    
    //Systematics
    for (int l = 1; l < nrealsys; l++)
    {
      lstream >> fSys[i][l].add;
      fSys[i][l].mult = fSys[i][l].add*100/fData[i];
      fSys[i][l].type = ADD;
      ostringstream sysname;
      sysname << "CHORUSSYS_1" << l;
      fSys[i][l].name = sysname.str();
    }

    //Get proton central value
    getline(f3,line);
    istringstream pstream(line);
    string sdum;
    int idum;
    double ddum;
    double proton_cv;
    pstream >> sdum >> sdum >> idum >> ddum >> proton_cv;
    
    //Get nuclear replicas
    getline(f2,line);
    istringstream nstream(line);
    nstream >> sdum >> sdum >> idum >> ddum >> ddum;
    
    vector<double> nuclear_cv (nrep);
      
    for(int irep=0; irep<nrep; irep++)
      {
	nstream >> nuclear_cv[irep];
      }

    //Compute additional uncertainties
    for(int l=nrealsys; l<fNSys; l++)
      {
	fSys[i][l].add = (nuclear_cv[l-nrealsys] - proton_cv)/sqrt(nrep);
	fSys[i][l].mult = fSys[i][l].add*100/fData[i];
	fSys[i][l].type = ADD;
	ostringstream sysname;
	sysname << "NUCLEAR" << l-nrealsys;
	fSys[i][l].name = sysname.str();
      }

  }
  
  f1.close();
  f2.close();
  f3.close();
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
void CHORUSNBPb_dw_iteFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/CHORUSNBPb/x-sec_shift_nb.txt";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream nuclearfile("");
  nuclearfile << dataPath() << "rawdata/CHORUSNBPb/nuclear/output/tables/group_result_table.csv";
  f2.open(nuclearfile.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << nuclearfile.str() << endl;
    exit(-1);
  }
  
  stringstream protonfile("");
  protonfile << dataPath() << "rawdata/CHORUSNBPb/proton_ite/output/tables/group_result_table.csv";
  f3.open(protonfile.str().c_str(), ios::in);
  
  if (f3.fail()) {
    cerr << "Error opening data file " << protonfile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double Mn = 0.9389;
  double enu,tmp,nortmp;
  int nrep=1000;
  int nrealsys=14;
  
  // Filtering data
  string line;

  getline(f2,line);
  getline(f3,line);
  
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
    lstream >> nortmp;    

    //QED radiation correction interpreted as uncertainty
    lstream >> nortmp;
    fSys[i][0].mult = (1.0-nortmp)*100.0;
    fSys[i][0].add = fSys[i][0].mult*fData[i]*1e-2;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CHORUSQEDRADCOR1";
    
    //Systematics
    for (int l = 1; l < nrealsys; l++)
    {
      lstream >> fSys[i][l].add;
      fSys[i][l].mult = fSys[i][l].add*100/fData[i];
      fSys[i][l].type = ADD;
      ostringstream sysname;
      sysname << "CHORUSSYS_1" << l;
      fSys[i][l].name = sysname.str();
    }

    //Get proton central value
    getline(f3,line);
    istringstream pstream(line);
    string sdum;
    int idum;
    double ddum;
    double proton_cv;
    pstream >> sdum >> sdum >> idum >> ddum >> proton_cv;

    //Get nuclear replicas
    getline(f2,line);
    istringstream nstream(line);
    nstream >> sdum >> sdum >> idum >> ddum >> ddum;
    vector<double> nuclear_cv (nrep);
      
    for(int irep=0; irep<nrep; irep++)
      {
	nstream >> nuclear_cv[irep];
      }

    //Compute additional uncertainties
    for(int l=nrealsys; l<fNSys; l++)
      {
	fSys[i][l].add = (nuclear_cv[l-nrealsys] - proton_cv)/sqrt(nrep);
	fSys[i][l].mult = fSys[i][l].add*100/fData[i];
	fSys[i][l].type = ADD;
	ostringstream sysname;
	sysname << "NUCLEAR" << l-nrealsys;
	fSys[i][l].name = sysname.str();
      }
  }
  
  f1.close();
  f2.close();
  f3.close();
  
}

void CHORUSNUPb_sh_iteFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/CHORUSNUPb/x-sec_shift_nu.txt";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream nuclearfile("");
  nuclearfile << dataPath() << "rawdata/CHORUSNUPb/nuclear/output/tables/group_result_table.csv";
  f2.open(nuclearfile.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << nuclearfile.str() << endl;
    exit(-1);
  }
  
  stringstream protonfile("");
  protonfile << dataPath() << "rawdata/CHORUSNUPb/proton_ite/output/tables/group_result_table.csv";
  f3.open(protonfile.str().c_str(), ios::in);
  
  if (f3.fail()) {
    cerr << "Error opening data file " << protonfile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double Mn = 0.9389;
  double enu,tmp,nortmp;
  int nrep=1000;
  int nrealsys=14;
  
  // Filtering data
  string line;

  getline(f2,line);
  getline(f3,line);
  
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
    lstream >> nortmp;   

    //QED radiation correction interpreted as uncertainty
    lstream >> nortmp;
    fSys[i][0].mult = (1.0-nortmp)*100.0;
    fSys[i][0].add = fSys[i][0].mult*fData[i]*1e-2;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CHORUSQEDRADCOR1";
    
    //Systematics
    for (int l = 1; l < nrealsys; l++)
    {
      lstream >> fSys[i][l].add;
      fSys[i][l].mult = fSys[i][l].add*100/fData[i];
      fSys[i][l].type = ADD;
      ostringstream sysname;
      sysname << "CHORUSSYS_1" << l;
      fSys[i][l].name = sysname.str();
    }

    //Get proton central value
    getline(f3,line);
    istringstream pstream(line);
    string sdum;
    int idum;
    double ddum;
    double proton_cv;
    pstream >> sdum >> sdum >> idum >> ddum >> proton_cv;

    //Get nuclear replicas
    getline(f2,line);
    istringstream nstream(line);
    double nuclear;
    nstream >> sdum >> sdum >> idum >> ddum >> nuclear;
    vector<double> nuclear_cv (nrep);
      
    for(int irep=0; irep<nrep; irep++)
      {
	nstream >> nuclear_cv[irep];
      }

    //Compute additional uncertainties
    for(int l=nrealsys; l<fNSys; l++)
      {
	fSys[i][l].add = (nuclear_cv[l-nrealsys] - nuclear)/sqrt(nrep);
	fSys[i][l].mult = fSys[i][l].add*100/fData[i];
	fSys[i][l].type = ADD;
	ostringstream sysname;
	sysname << "NUCLEAR" << l-nrealsys;
	fSys[i][l].name = sysname.str();
      }

    //Compute shifts
    //cout << nuclear/proton_cv << "   " << 0.0 << endl;
    
  }
  
  f1.close();
  f2.close();
  f3.close();
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
void CHORUSNBPb_sh_iteFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/CHORUSNBPb/x-sec_shift_nb.txt";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream nuclearfile("");
  nuclearfile << dataPath() << "rawdata/CHORUSNBPb/nuclear/output/tables/group_result_table.csv";
  f2.open(nuclearfile.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << nuclearfile.str() << endl;
    exit(-1);
  }
  
  stringstream protonfile("");
  protonfile << dataPath() << "rawdata/CHORUSNBPb/proton_ite/output/tables/group_result_table.csv";
  f3.open(protonfile.str().c_str(), ios::in);
  
  if (f3.fail()) {
    cerr << "Error opening data file " << protonfile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double Mn = 0.9389;
  double enu,tmp,nortmp;
  int nrep=1000;
  int nrealsys=14;
  
  // Filtering data
  string line;

  getline(f2,line);
  getline(f3,line);
  
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
    lstream >> nortmp;    

    //QED radiation correction interpreted as uncertainty
    lstream >> nortmp;
    fSys[i][0].mult = (1.0-nortmp)*100.0;
    fSys[i][0].add = fSys[i][0].mult*fData[i]*1e-2;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CHORUSQEDRADCOR1";
    
    //Systematics
    for (int l = 1; l < nrealsys; l++)
    {
      lstream >> fSys[i][l].add;
      fSys[i][l].mult = fSys[i][l].add*100/fData[i];
      fSys[i][l].type = ADD;
      ostringstream sysname;
      sysname << "CHORUSSYS_1" << l;
      fSys[i][l].name = sysname.str();
    }

    //Get proton central value
    getline(f3,line);
    istringstream pstream(line);
    string sdum;
    int idum;
    double ddum;
    double proton_cv;
    pstream >> sdum >> sdum >> idum >> ddum >> proton_cv;

    //Get nuclear replicas
    getline(f2,line);
    istringstream nstream(line);
    double nuclear;
    nstream >> sdum >> sdum >> idum >> ddum >> nuclear;
    vector<double> nuclear_cv (nrep);
      
    for(int irep=0; irep<nrep; irep++)
      {
	nstream >> nuclear_cv[irep];
      }

    //Compute additional uncertainties
    for(int l=nrealsys; l<fNSys; l++)
      {
	fSys[i][l].add = (nuclear_cv[l-nrealsys] - nuclear)/sqrt(nrep);
	fSys[i][l].mult = fSys[i][l].add*100/fData[i];
	fSys[i][l].type = ADD;
	ostringstream sysname;
	sysname << "NUCLEAR" << l-nrealsys;
	fSys[i][l].name = sysname.str();
      }
    
    //Compute shifts
    //cout << nuclear/proton_cv << "   " << 0.0 << endl;
    
  }
  
  f1.close();
  f2.close();
  f3.close();
  
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
void CHORUSNUPb_dw_30Filter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/CHORUSNUPb/x-sec_shift_nu.txt";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream nuclearfile("");
  nuclearfile << dataPath() << "rawdata/CHORUSNUPb/nuclear_30/output/tables/group_result_table.csv";
  f2.open(nuclearfile.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << nuclearfile.str() << endl;
    exit(-1);
  }
  
  stringstream protonfile("");
  protonfile << dataPath() << "rawdata/CHORUSNUPb/proton_30/output/tables/group_result_table.csv";
  f3.open(protonfile.str().c_str(), ios::in);
  
  if (f3.fail()) {
    cerr << "Error opening data file " << protonfile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double Mn = 0.9389;
  double enu,tmp,nortmp;
  int nrep=200;
  int nrealsys=14;
  
  // Filtering data
  string line;

  getline(f2,line);
  getline(f3,line);
  
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
    lstream >> nortmp;   

    //QED radiation correction interpreted as uncertainty
    lstream >> nortmp;
    fSys[i][0].mult = (1.0-nortmp)*100.0;
    fSys[i][0].add = fSys[i][0].mult*fData[i]*1e-2;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CHORUSQEDRADCOR1";
    
    //Systematics
    for (int l = 1; l < nrealsys; l++)
    {
      lstream >> fSys[i][l].add;
      fSys[i][l].mult = fSys[i][l].add*100/fData[i];
      fSys[i][l].type = ADD;
      ostringstream sysname;
      sysname << "CHORUSSYS_1" << l;
      fSys[i][l].name = sysname.str();
    }

    //Get proton central value
    getline(f3,line);
    istringstream pstream(line);
    string sdum;
    int idum;
    double ddum;
    double proton_cv;
    pstream >> sdum >> sdum >> idum >> ddum >> proton_cv;
    
    //Get nuclear replicas
    getline(f2,line);
    istringstream nstream(line);
    nstream >> sdum >> sdum >> idum >> ddum >> ddum;
    
    vector<double> nuclear_cv (nrep);
      
    for(int irep=0; irep<nrep; irep++)
      {
	nstream >> nuclear_cv[irep];
      }

    //Compute additional uncertainties
    for(int l=nrealsys; l<fNSys; l++)
      {
	fSys[i][l].add = (nuclear_cv[l-nrealsys] - proton_cv)/sqrt(nrep);
	fSys[i][l].mult = fSys[i][l].add*100/fData[i];
	fSys[i][l].type = ADD;
	ostringstream sysname;
	sysname << "NUCLEAR" << l-nrealsys;
	fSys[i][l].name = sysname.str();
      }

  }
  
  f1.close();
  f2.close();
  f3.close();
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
void CHORUSNBPb_dw_30Filter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/CHORUSNBPb/x-sec_shift_nb.txt";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream nuclearfile("");
  nuclearfile << dataPath() << "rawdata/CHORUSNBPb/nuclear_30/output/tables/group_result_table.csv";
  f2.open(nuclearfile.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << nuclearfile.str() << endl;
    exit(-1);
  }
  
  stringstream protonfile("");
  protonfile << dataPath() << "rawdata/CHORUSNBPb/proton_30/output/tables/group_result_table.csv";
  f3.open(protonfile.str().c_str(), ios::in);
  
  if (f3.fail()) {
    cerr << "Error opening data file " << protonfile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double Mn = 0.9389;
  double enu,tmp,nortmp;
  int nrep=200;
  int nrealsys=14;
  
  // Filtering data
  string line;

  getline(f2,line);
  getline(f3,line);
  
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
    lstream >> nortmp;    

    //QED radiation correction interpreted as uncertainty
    lstream >> nortmp;
    fSys[i][0].mult = (1.0-nortmp)*100.0;
    fSys[i][0].add = fSys[i][0].mult*fData[i]*1e-2;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CHORUSQEDRADCOR1";
    
    //Systematics
    for (int l = 1; l < nrealsys; l++)
    {
      lstream >> fSys[i][l].add;
      fSys[i][l].mult = fSys[i][l].add*100/fData[i];
      fSys[i][l].type = ADD;
      ostringstream sysname;
      sysname << "CHORUSSYS_1" << l;
      fSys[i][l].name = sysname.str();
    }

    //Get proton central value
    getline(f3,line);
    istringstream pstream(line);
    string sdum;
    int idum;
    double ddum;
    double proton_cv;
    pstream >> sdum >> sdum >> idum >> ddum >> proton_cv;

    //Get nuclear replicas
    getline(f2,line);
    istringstream nstream(line);
    nstream >> sdum >> sdum >> idum >> ddum >> ddum;
    vector<double> nuclear_cv (nrep);
      
    for(int irep=0; irep<nrep; irep++)
      {
	nstream >> nuclear_cv[irep];
      }

    //Compute additional uncertainties
    for(int l=nrealsys; l<fNSys; l++)
      {
	fSys[i][l].add = (nuclear_cv[l-nrealsys] - proton_cv)/sqrt(nrep);
	fSys[i][l].mult = fSys[i][l].add*100/fData[i];
	fSys[i][l].type = ADD;
	ostringstream sysname;
	sysname << "NUCLEAR" << l-nrealsys;
	fSys[i][l].name = sysname.str();
      }
  }
  
  f1.close();
  f2.close();
  f3.close();
  
}

void CHORUSNUPb_sh_30Filter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/CHORUSNUPb/x-sec_shift_nu.txt";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream nuclearfile("");
  nuclearfile << dataPath() << "rawdata/CHORUSNUPb/nuclear_30/output/tables/group_result_table.csv";
  f2.open(nuclearfile.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << nuclearfile.str() << endl;
    exit(-1);
  }
  
  stringstream protonfile("");
  protonfile << dataPath() << "rawdata/CHORUSNUPb/proton_30/output/tables/group_result_table.csv";
  f3.open(protonfile.str().c_str(), ios::in);
  
  if (f3.fail()) {
    cerr << "Error opening data file " << protonfile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double Mn = 0.9389;
  double enu,tmp,nortmp;
  int nrep=200;
  int nrealsys=14;
  
  // Filtering data
  string line;

  getline(f2,line);
  getline(f3,line);
  
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
    lstream >> nortmp;   

    //QED radiation correction interpreted as uncertainty
    lstream >> nortmp;
    fSys[i][0].mult = (1.0-nortmp)*100.0;
    fSys[i][0].add = fSys[i][0].mult*fData[i]*1e-2;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CHORUSQEDRADCOR1";
    
    //Systematics
    for (int l = 1; l < nrealsys; l++)
    {
      lstream >> fSys[i][l].add;
      fSys[i][l].mult = fSys[i][l].add*100/fData[i];
      fSys[i][l].type = ADD;
      ostringstream sysname;
      sysname << "CHORUSSYS_1" << l;
      fSys[i][l].name = sysname.str();
    }

    //Get proton central value
    getline(f3,line);
    istringstream pstream(line);
    string sdum;
    int idum;
    double ddum;
    double proton_cv;
    pstream >> sdum >> sdum >> idum >> ddum >> proton_cv;

    //Get nuclear replicas
    getline(f2,line);
    istringstream nstream(line);
    double nuclear;
    nstream >> sdum >> sdum >> idum >> ddum >> nuclear;
    vector<double> nuclear_cv (nrep);
      
    for(int irep=0; irep<nrep; irep++)
      {
	nstream >> nuclear_cv[irep];
      }

    //Compute additional uncertainties
    for(int l=nrealsys; l<fNSys; l++)
      {
	fSys[i][l].add = (nuclear_cv[l-nrealsys] - nuclear)/sqrt(nrep);
	fSys[i][l].mult = fSys[i][l].add*100/fData[i];
	fSys[i][l].type = ADD;
	ostringstream sysname;
	sysname << "NUCLEAR" << l-nrealsys;
	fSys[i][l].name = sysname.str();
      }

    //Compute shifts
    //cout << nuclear/proton_cv << "   " << 0.0 << endl;
    
  }
  
  f1.close();
  f2.close();
  f3.close();
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
void CHORUSNBPb_sh_30Filter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/CHORUSNBPb/x-sec_shift_nb.txt";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream nuclearfile("");
  nuclearfile << dataPath() << "rawdata/CHORUSNBPb/nuclear_30/output/tables/group_result_table.csv";
  f2.open(nuclearfile.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << nuclearfile.str() << endl;
    exit(-1);
  }
  
  stringstream protonfile("");
  protonfile << dataPath() << "rawdata/CHORUSNBPb/proton_30/output/tables/group_result_table.csv";
  f3.open(protonfile.str().c_str(), ios::in);
  
  if (f3.fail()) {
    cerr << "Error opening data file " << protonfile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double Mn = 0.9389;
  double enu,tmp,nortmp;
  int nrep=200;
  int nrealsys=14;
  
  // Filtering data
  string line;

  getline(f2,line);
  getline(f3,line);
  
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
    lstream >> nortmp;    

    //QED radiation correction interpreted as uncertainty
    lstream >> nortmp;
    fSys[i][0].mult = (1.0-nortmp)*100.0;
    fSys[i][0].add = fSys[i][0].mult*fData[i]*1e-2;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CHORUSQEDRADCOR1";
    
    //Systematics
    for (int l = 1; l < nrealsys; l++)
    {
      lstream >> fSys[i][l].add;
      fSys[i][l].mult = fSys[i][l].add*100/fData[i];
      fSys[i][l].type = ADD;
      ostringstream sysname;
      sysname << "CHORUSSYS_1" << l;
      fSys[i][l].name = sysname.str();
    }

    //Get proton central value
    getline(f3,line);
    istringstream pstream(line);
    string sdum;
    int idum;
    double ddum;
    double proton_cv;
    pstream >> sdum >> sdum >> idum >> ddum >> proton_cv;

    //Get nuclear replicas
    getline(f2,line);
    istringstream nstream(line);
    double nuclear;
    nstream >> sdum >> sdum >> idum >> ddum >> nuclear;
    vector<double> nuclear_cv (nrep);
      
    for(int irep=0; irep<nrep; irep++)
      {
	nstream >> nuclear_cv[irep];
      }

    //Compute additional uncertainties
    for(int l=nrealsys; l<fNSys; l++)
      {
	fSys[i][l].add = (nuclear_cv[l-nrealsys] - nuclear)/sqrt(nrep);
	fSys[i][l].mult = fSys[i][l].add*100/fData[i];
	fSys[i][l].type = ADD;
	ostringstream sysname;
	sysname << "NUCLEAR" << l-nrealsys;
	fSys[i][l].name = sysname.str();
      }
    
    //Compute shifts
    //cout << nuclear/proton_cv << "   " << 0.0 << endl;
    
  }
  
  f1.close();
  f2.close();
  f3.close();
  
}
