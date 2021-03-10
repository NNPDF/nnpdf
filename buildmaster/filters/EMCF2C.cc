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
 * EMCF2C.cc
 *
 *  Production of Charmed Particles in 250-GeV mu+-iron interactions
 *  The European Muon Collaboration (J.J. Auber et al.)
 *  Nucl. Phys. B213 (1983) 31
 *
 */

#include "EMCF2C.h"

void EMCF2CFilter::ReadData()
{
  // Opening files
  fstream f1;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/EMCF2C.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  /*
  // Starting filter
  double syscor[fNData][fNSys];

  for (int i = 0; i < fNData; i++)
  {
    for (int l = 0; l < fNSys; l++)
      syscor[i][l] = 0.0;
  }
  */

/*
  // Reading data
  index = 0;
  string tmp;
  getline(f1,tmp);
  {
    for (int i = 0; i < fNData; i++)
    for (int j = 0; j < 52; j++)
      f1 >> datain[i][j];
    if (datain[i][0] != i+1) {
      cerr << "Mismatch" << endl;
      exit(-1);
    }
    index++;
  }

  if (index != fNData)
  {
    cerr << "Mismatch in the number of data points" << endl;
    exit(-1);
  }
*/
  // Filtering data
  double x[fNData], y[fNData], q2[fNData];
  double xmin, xmax, q2min, q2max;


  string tmp;
  getline(f1,tmp);
  getline(f1,tmp);
  for (int i = 0; i < fNData; i++)
  {
    f1 >> xmin >> xmax >> q2min >> q2max >> fData[i] >> fStat[i];
    x[i]   = (xmin + xmax)/2.;
    q2[i]  = (q2min + q2max)/2.;
    y[i]   = 0.0;

    // Uncorrelated systematics
//    fSys[i][1].mult = datain[i][6];
//    fSys[i][1].type = ADD;
//    fSys[i][1].name = "UNCORR";

    // Systematic errors
//    double sys_tot = 0.0;
//    for (int l = 0; l < fNSys-2; l++)
//    {
//      syscor[i][l] = datain[i][l+8];
//      sys_tot += syscor[i][l]*syscor[i][l];
//    }
//    sys_tot = sqrt(sys_tot);

    // Normalization uncertainty
//    fSys[i][0].mult = absnorexp;
//    fSys[i][0].type = MULT;
//    fSys[i][0].name = "CORR";

    // Check total error
//    double tot = datain[i][7];
//    double tot_tmp = sqrt(sys_tot*sys_tot + stat[i]*stat[i] + fSys[i][1].mult*fSys[i][1].mult);
//    if (fabs(tot-tot_tmp) > 1.4)
//    {
//      cerr << "Mismatch " << i << "\t" << tot << "\t" << tot_tmp << endl;
//      exit(-1);
//    }

    fKin1[i] = x[i];
    fKin2[i] = q2[i];
    fKin3[i] = y[i];

    //rescaling for BR - check
    fData[i] = fData[i]*0.82;
    double sist = fData[i]*0.15;
    fSys[i][0].add = sist;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CORR_EMC";
    fSys[i][0].mult = 15;
  }

  f1.close();
}

void EMCF2C_dwFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/EMCF2C/EMCF2C.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream nuclearfile("");
  nuclearfile << dataPath() << "rawdata/EMCF2C/nuclear/output/tables/group_result_table.csv";
  f2.open(nuclearfile.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << nuclearfile.str() << endl;
    exit(-1);
  }

  stringstream protonfile("");
  protonfile << dataPath() << "rawdata/EMCF2C/proton/output/tables/group_result_table.csv";
  f3.open(protonfile.str().c_str(), ios::in);
  
  if (f3.fail()) {
    cerr << "Error opening data file " << protonfile.str() << endl;
    exit(-1);
  }
  
  // Filtering data
  double x[fNData], y[fNData], q2[fNData];
  double xmin, xmax, q2min, q2max;
  string line;

  getline(f2,line);
  getline(f3,line);

  string tmp;
  getline(f1,tmp);
  getline(f1,tmp);
  for (int i = 0; i < fNData; i++)
  {
    f1 >> xmin >> xmax >> q2min >> q2max >> fData[i] >> fStat[i];
    x[i]   = (xmin + xmax)/2.;
    q2[i]  = (q2min + q2max)/2.;
    y[i]   = 0.0;

    fKin1[i] = x[i];
    fKin2[i] = q2[i];
    fKin3[i] = y[i];

    //rescaling for BR - check
    fData[i] = fData[i]*0.82;
    double sist = fData[i]*0.15;
    fSys[i][0].add = sist;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CORREMC";
    fSys[i][0].mult = 15;

    int nrep=1000;
    int nrealsys=1;
    
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

void EMCF2C_shFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/EMCF2C/EMCF2C.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream nuclearfile("");
  nuclearfile << dataPath() << "rawdata/EMCF2C/nuclear/output/tables/group_result_table.csv";
  f2.open(nuclearfile.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << nuclearfile.str() << endl;
    exit(-1);
  }

  stringstream protonfile("");
  protonfile << dataPath() << "rawdata/EMCF2C/proton/output/tables/group_result_table.csv";
  f3.open(protonfile.str().c_str(), ios::in);
  
  if (f3.fail()) {
    cerr << "Error opening data file " << protonfile.str() << endl;
    exit(-1);
  }
  
  // Filtering data
  double x[fNData], y[fNData], q2[fNData];
  double xmin, xmax, q2min, q2max;
  string line;

  getline(f2,line);
  getline(f3,line);

  string tmp;
  getline(f1,tmp);
  getline(f1,tmp);
  for (int i = 0; i < fNData; i++)
  {
    f1 >> xmin >> xmax >> q2min >> q2max >> fData[i] >> fStat[i];
    x[i]   = (xmin + xmax)/2.;
    q2[i]  = (q2min + q2max)/2.;
    y[i]   = 0.0;

    fKin1[i] = x[i];
    fKin2[i] = q2[i];
    fKin3[i] = y[i];

    //rescaling for BR - check
    fData[i] = fData[i]*0.82;
    double sist = fData[i]*0.15;
    fSys[i][0].add = sist;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CORREMC";
    fSys[i][0].mult = 15;

    int nrep=1000;
    int nrealsys=1;
    
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

void EMCF2C_dw_iteFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/EMCF2C/EMCF2C.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream nuclearfile("");
  nuclearfile << dataPath() << "rawdata/EMCF2C/nuclear/output/tables/group_result_table.csv";
  f2.open(nuclearfile.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << nuclearfile.str() << endl;
    exit(-1);
  }

  stringstream protonfile("");
  protonfile << dataPath() << "rawdata/EMCF2C/proton_ite/output/tables/group_result_table.csv";
  f3.open(protonfile.str().c_str(), ios::in);
  
  if (f3.fail()) {
    cerr << "Error opening data file " << protonfile.str() << endl;
    exit(-1);
  }
  
  // Filtering data
  double x[fNData], y[fNData], q2[fNData];
  double xmin, xmax, q2min, q2max;
  string line;

  getline(f2,line);
  getline(f3,line);

  string tmp;
  getline(f1,tmp);
  getline(f1,tmp);
  for (int i = 0; i < fNData; i++)
  {
    f1 >> xmin >> xmax >> q2min >> q2max >> fData[i] >> fStat[i];
    x[i]   = (xmin + xmax)/2.;
    q2[i]  = (q2min + q2max)/2.;
    y[i]   = 0.0;

    fKin1[i] = x[i];
    fKin2[i] = q2[i];
    fKin3[i] = y[i];

    //rescaling for BR - check
    fData[i] = fData[i]*0.82;
    double sist = fData[i]*0.15;
    fSys[i][0].add = sist;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CORREMC";
    fSys[i][0].mult = 15;

    int nrep=1000;
    int nrealsys=1;
    
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

void EMCF2C_sh_iteFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/EMCF2C/EMCF2C.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream nuclearfile("");
  nuclearfile << dataPath() << "rawdata/EMCF2C/nuclear/output/tables/group_result_table.csv";
  f2.open(nuclearfile.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << nuclearfile.str() << endl;
    exit(-1);
  }

  stringstream protonfile("");
  protonfile << dataPath() << "rawdata/EMCF2C/proton_ite/output/tables/group_result_table.csv";
  f3.open(protonfile.str().c_str(), ios::in);
  
  if (f3.fail()) {
    cerr << "Error opening data file " << protonfile.str() << endl;
    exit(-1);
  }
  
  // Filtering data
  double x[fNData], y[fNData], q2[fNData];
  double xmin, xmax, q2min, q2max;
  string line;

  getline(f2,line);
  getline(f3,line);

  string tmp;
  getline(f1,tmp);
  getline(f1,tmp);
  for (int i = 0; i < fNData; i++)
  {
    f1 >> xmin >> xmax >> q2min >> q2max >> fData[i] >> fStat[i];
    x[i]   = (xmin + xmax)/2.;
    q2[i]  = (q2min + q2max)/2.;
    y[i]   = 0.0;

    fKin1[i] = x[i];
    fKin2[i] = q2[i];
    fKin3[i] = y[i];

    //rescaling for BR - check
    fData[i] = fData[i]*0.82;
    double sist = fData[i]*0.15;
    fSys[i][0].add = sist;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CORREMC";
    fSys[i][0].mult = 15;

    int nrep=1000;
    int nrealsys=1;
    
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
