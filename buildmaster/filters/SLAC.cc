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
 *     EXPeriment      = SLAC E49, E61, E87, E89
 *     REACtion        = e proton --> e X
 *     Author          = Whitlow
 *     REFerence       = PL B282(92)475 and SLAC-357 (1990) (Ph.D)
 *     Additional info : Re-analysed data - BCDMS binning
 *     The proton structure functions from a re-analysis of the
 *     SLAC experiments E49, E61, E87, and E89 in electron proton
 *     deep inelastic scattering binned in the same x ranges as
 *     the BCDMS data.
 *
 *     mean          mean            F2              errors
 *     x           Q**2                       stat.      sys.
 *
 *     EXPeriment      = SLAC E49, E61, E87, E89
 *     REACtion        = e deuterium --> e X
 *     Plab            = 4.50-24.50 GeV
 *     Author          = Whitlow
 *     REFerence       = SLAC-357 (1990) (Ph.D)
 *     Additional info : Re-analysed data - BCDMS binning
 *     The nucleon structure functions from a re-analysis of the
 *     SLAC experiments E49, E61, E87, E89, E137 and E140 in electron
 *     deuterium deep inelastic scattering binned in the same
 *     x ranges as the BCDMS data.
 *
 *     x           Q**2             F2             errors
 *     (GeV**2)                     stat.      sys.
 *
 *     The systematic error is taken fully correlated for all the data points p and d
 *
 *     The absolute normalization is taken to be 2.1% for p and 1.7% for d
 *     The relative normalization of 1.1% is taken into account
 *
 *     NB: There are not info about y and since the data points that were taken at
 *     diffferent beam energies are merged together in the same sample
 *     there are no way to evaluate it. Thus it's set arbitrary equal to 0
 *     simply to fill somehow the entry in datawarehouse.res
 */

#include "SLAC.h"

void SLACPFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/slac_p.data";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  double relnor;
  
  relnor = 1.1*0.5;    //relative normalisation of 1.1% between targets
  
  string line;
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    lstream >> fKin1[i];  //x
    lstream >> fKin2[i];  //q2
    lstream >> fData[i];  //obs 
    fKin3[i] = 0;
    
    //  SLAC gives errors in absolute value
    //  and we assume the sys errors to be uncorrelated in order
    //  to overestimate them and give less weight
    //  to these data points in evaluating chi2

    lstream >> fStat[i];
    
    lstream >> fSys[i][0].add;
    fSys[i][0].mult = fSys[i][0].add*100/fData[i];
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    
    fSys[i][1].mult = 2.1;  //absnorm
    fSys[i][1].add = fSys[i][1].mult*fData[i]*1e-2;
    fSys[i][1].type = MULT;
    fSys[i][1].name = "SLACNORM";
    
    fSys[i][2].mult = relnor;     //relnorm
    fSys[i][2].add = fSys[i][2].mult*fData[i]*1e-2;
    fSys[i][2].type = MULT;
    fSys[i][2].name = "SLACRELNORM";
  }
  
  f1.close();
}


/**
 *     EXPeriment      = SLAC E49, E61, E87, E89
 *     REACtion        = e proton --> e X
 *     Author          = Whitlow
 *     REFerence       = PL B282(92)475 and SLAC-357 (1990) (Ph.D)
 *     Additional info : Re-analysed data - BCDMS binning
 *     The proton structure functions from a re-analysis of the
 *     SLAC experiments E49, E61, E87, and E89 in electron proton
 *     deep inelastic scattering binned in the same x ranges as
 *     the BCDMS data.
 *
 *     mean          mean            F2              errors
 *     x           Q**2                       stat.      sys.
 *
 *     EXPeriment      = SLAC E49, E61, E87, E89
 *     REACtion        = e deuterium --> e X
 *     Plab            = 4.50-24.50 GeV
 *     Author          = Whitlow
 *     REFerence       = SLAC-357 (1990) (Ph.D)
 *     Additional info : Re-analysed data - BCDMS binning
 *     The nucleon structure functions from a re-analysis of the
 *     SLAC experiments E49, E61, E87, E89, E137 and E140 in electron
 *     deuterium deep inelastic scattering binned in the same
 *     x ranges as the BCDMS data.
 *
 *     x           Q**2             F2             errors
 *     (GeV**2)                     stat.      sys.
 *
 *     The systematic error is taken fully correlated for all the data points p and d
 *
 *     The absolute normalization is taken to be 2.1% for p and 1.7% for d
 *     The relative normalization of 1.1% is taken into account
 *
 *     NB: There are not info about y and since the data points that were taken at
 *     diffferent beam energies are merged together in the same sample
 *     there are no way to evaluate it. Thus it's set arbitrary equal to 0
 *     simply to fill somehow the entry in datawarehouse.res
 */
void SLACDFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/slac_d.data";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  double relnor;
  
  relnor = 1.1*0.5;    //relative normalisation of 1.1% between targets
  
  string line;
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    lstream >> fKin1[i];  //x
    lstream >> fKin2[i];  //q2
    lstream >> fData[i];  //obs 
    fKin3[i] = 0;
    
    //  SLAC gives errors in absolute value
    //  and we assume the sys errors to be uncorrelated in order
    //  to overestimate them and give less weight
    //  to these data points in evaluating chi2
   
    lstream >> fStat[i];
    
    lstream >> fSys[i][0].add;
    fSys[i][0].mult = fSys[i][0].add*100/fData[i];
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    
    fSys[i][1].mult = 1.7;  //absnorm
    fSys[i][1].add = fSys[i][1].mult*fData[i]*1e-2;
    fSys[i][1].type = MULT;
    fSys[i][1].name = "SLACNORM";
    
    fSys[i][2].mult = -relnor;     //relnorm
    fSys[i][2].add = fSys[i][2].mult*fData[i]*1e-2;
    fSys[i][2].type = MULT;
    fSys[i][2].name = "SLACRELNORM";
  }
  
  f1.close();
}

void SLACD_dwFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/SLACD/slac_d.data";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream nuclearfile("");
  nuclearfile << dataPath() << "rawdata/SLACD/nuclear/output/tables/group_result_table.csv";
  f2.open(nuclearfile.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << nuclearfile.str() << endl;
    exit(-1);
  }
  
  stringstream protonfile("");
  protonfile << dataPath() << "rawdata/SLACD/proton/output/tables/group_result_table.csv";
  f3.open(protonfile.str().c_str(), ios::in);
  
  if (f3.fail()) {
    cerr << "Error opening data file " << protonfile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double relnor;
  int nrep=100;
  int nrealsys=3;
  
  relnor = 1.1*0.5;    //relative normalisation of 1.1% between targets
  
  string line;

  getline(f2,line);
  getline(f3,line);
  
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    lstream >> fKin1[i];  //x
    lstream >> fKin2[i];  //q2
    lstream >> fData[i];  //obs 
    fKin3[i] = 0;
    
    //  SLAC gives errors in absolute value
    //  and we assume the sys errors to be uncorrelated in order
    //  to overestimate them and give less weight
    //  to these data points in evaluating chi2
   
    lstream >> fStat[i];
    
    lstream >> fSys[i][0].add;
    fSys[i][0].mult = fSys[i][0].add*100/fData[i];
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    
    fSys[i][1].mult = 1.7;  //absnorm
    fSys[i][1].add = fSys[i][1].mult*fData[i]*1e-2;
    fSys[i][1].type = MULT;
    fSys[i][1].name = "SLACNORM1";
    
    fSys[i][2].mult = -relnor;     //relnorm
    fSys[i][2].add = fSys[i][2].mult*fData[i]*1e-2;
    fSys[i][2].type = MULT;
    fSys[i][2].name = "SLACRELNORM1";

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
	sysname << "DEUTERON" << l-nrealsys;
	fSys[i][l].name = sysname.str();
      }
    
  }
  
  f1.close();
  f2.close();
  f3.close();
}

void SLACD_shFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/SLACD/slac_d.data";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream nuclearfile("");
  nuclearfile << dataPath() << "rawdata/SLACD/nuclear/output/tables/group_result_table.csv";
  f2.open(nuclearfile.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << nuclearfile.str() << endl;
    exit(-1);
  }
  
  stringstream protonfile("");
  protonfile << dataPath() << "rawdata/SLACD/proton/output/tables/group_result_table.csv";
  f3.open(protonfile.str().c_str(), ios::in);
  
  if (f3.fail()) {
    cerr << "Error opening data file " << protonfile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double relnor;
  int nrep=100;
  int nrealsys=3;
  
  relnor = 1.1*0.5;    //relative normalisation of 1.1% between targets
  
  string line;

  getline(f2,line);
  getline(f3,line);
  
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    lstream >> fKin1[i];  //x
    lstream >> fKin2[i];  //q2
    lstream >> fData[i];  //obs 
    fKin3[i] = 0;
    
    //  SLAC gives errors in absolute value
    //  and we assume the sys errors to be uncorrelated in order
    //  to overestimate them and give less weight
    //  to these data points in evaluating chi2
   
    lstream >> fStat[i];
    
    lstream >> fSys[i][0].add;
    fSys[i][0].mult = fSys[i][0].add*100/fData[i];
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    
    fSys[i][1].mult = 1.7;  //absnorm
    fSys[i][1].add = fSys[i][1].mult*fData[i]*1e-2;
    fSys[i][1].type = MULT;
    fSys[i][1].name = "SLACNORM1";
    
    fSys[i][2].mult = -relnor;     //relnorm
    fSys[i][2].add = fSys[i][2].mult*fData[i]*1e-2;
    fSys[i][2].type = MULT;
    fSys[i][2].name = "SLACRELNORM1";

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
	sysname << "DEUTERON" << l-nrealsys;
	fSys[i][l].name = sysname.str();
      }

    //Compute shifts
    //cout << nuclear/proton_cv << "   " << 0.0 << endl;
    
  }
  
  f1.close();
  f2.close();
  f3.close();
}

void SLACD_dw_iteFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/SLACD/slac_d.data";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream nuclearfile("");
  nuclearfile << dataPath() << "rawdata/SLACD/nuclear_ite/output/tables/group_result_table.csv";
  f2.open(nuclearfile.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << nuclearfile.str() << endl;
    exit(-1);
  }
  
  stringstream protonfile("");
  protonfile << dataPath() << "rawdata/SLACD/proton_ite/output/tables/group_result_table.csv";
  f3.open(protonfile.str().c_str(), ios::in);
  
  if (f3.fail()) {
    cerr << "Error opening data file " << protonfile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double relnor;
  int nrep=100;
  int nrealsys=3;
  
  relnor = 1.1*0.5;    //relative normalisation of 1.1% between targets
  
  string line;

  getline(f2,line);
  getline(f3,line);
  
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    lstream >> fKin1[i];  //x
    lstream >> fKin2[i];  //q2
    lstream >> fData[i];  //obs 
    fKin3[i] = 0;
    
    //  SLAC gives errors in absolute value
    //  and we assume the sys errors to be uncorrelated in order
    //  to overestimate them and give less weight
    //  to these data points in evaluating chi2
   
    lstream >> fStat[i];
    
    lstream >> fSys[i][0].add;
    fSys[i][0].mult = fSys[i][0].add*100/fData[i];
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    
    fSys[i][1].mult = 1.7;  //absnorm
    fSys[i][1].add = fSys[i][1].mult*fData[i]*1e-2;
    fSys[i][1].type = MULT;
    fSys[i][1].name = "SLACNORM1";
    
    fSys[i][2].mult = -relnor;     //relnorm
    fSys[i][2].add = fSys[i][2].mult*fData[i]*1e-2;
    fSys[i][2].type = MULT;
    fSys[i][2].name = "SLACRELNORM1";

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
	sysname << "DEUTERON" << l-nrealsys;
	fSys[i][l].name = sysname.str();
      }
    
  }
  
  f1.close();
  f2.close();
  f3.close();
}

void SLACD_sh_iteFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/SLACD/slac_d.data";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream nuclearfile("");
  nuclearfile << dataPath() << "rawdata/SLACD/nuclear_ite/output/tables/group_result_table.csv";
  f2.open(nuclearfile.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << nuclearfile.str() << endl;
    exit(-1);
  }
  
  stringstream protonfile("");
  protonfile << dataPath() << "rawdata/SLACD/proton_ite/output/tables/group_result_table.csv";
  f3.open(protonfile.str().c_str(), ios::in);
  
  if (f3.fail()) {
    cerr << "Error opening data file " << protonfile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double relnor;
  int nrep=100;
  int nrealsys=3;
  
  relnor = 1.1*0.5;    //relative normalisation of 1.1% between targets
  
  string line;

  getline(f2,line);
  getline(f3,line);
  
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    lstream >> fKin1[i];  //x
    lstream >> fKin2[i];  //q2
    lstream >> fData[i];  //obs 
    fKin3[i] = 0;
    
    //  SLAC gives errors in absolute value
    //  and we assume the sys errors to be uncorrelated in order
    //  to overestimate them and give less weight
    //  to these data points in evaluating chi2
   
    lstream >> fStat[i];
    
    lstream >> fSys[i][0].add;
    fSys[i][0].mult = fSys[i][0].add*100/fData[i];
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    
    fSys[i][1].mult = 1.7;  //absnorm
    fSys[i][1].add = fSys[i][1].mult*fData[i]*1e-2;
    fSys[i][1].type = MULT;
    fSys[i][1].name = "SLACNORM1";
    
    fSys[i][2].mult = -relnor;     //relnorm
    fSys[i][2].add = fSys[i][2].mult*fData[i]*1e-2;
    fSys[i][2].type = MULT;
    fSys[i][2].name = "SLACRELNORM1";

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
	sysname << "DEUTERON" << l-nrealsys;
	fSys[i][l].name = sysname.str();
      }

    //Compute shifts
    //cout << nuclear/proton_cv << "   " << 0.0 << endl;
    
  }
  
  f1.close();
  f2.close();
  f3.close();
}

void SLACD_dw_30Filter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/SLACD/slac_d.data";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream nuclearfile("");
  nuclearfile << dataPath() << "rawdata/SLACD/nuclear_30/output/tables/group_result_table.csv";
  f2.open(nuclearfile.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << nuclearfile.str() << endl;
    exit(-1);
  }
  
  stringstream protonfile("");
  protonfile << dataPath() << "rawdata/SLACD/proton_30/output/tables/group_result_table.csv";
  f3.open(protonfile.str().c_str(), ios::in);
  
  if (f3.fail()) {
    cerr << "Error opening data file " << protonfile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double relnor;
  int nrep=200;
  int nrealsys=3;
  
  relnor = 1.1*0.5;    //relative normalisation of 1.1% between targets
  
  string line;

  getline(f2,line);
  getline(f3,line);
  
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    lstream >> fKin1[i];  //x
    lstream >> fKin2[i];  //q2
    lstream >> fData[i];  //obs 
    fKin3[i] = 0;
    
    //  SLAC gives errors in absolute value
    //  and we assume the sys errors to be uncorrelated in order
    //  to overestimate them and give less weight
    //  to these data points in evaluating chi2
   
    lstream >> fStat[i];
    
    lstream >> fSys[i][0].add;
    fSys[i][0].mult = fSys[i][0].add*100/fData[i];
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    
    fSys[i][1].mult = 1.7;  //absnorm
    fSys[i][1].add = fSys[i][1].mult*fData[i]*1e-2;
    fSys[i][1].type = MULT;
    fSys[i][1].name = "SLACNORM1";
    
    fSys[i][2].mult = -relnor;     //relnorm
    fSys[i][2].add = fSys[i][2].mult*fData[i]*1e-2;
    fSys[i][2].type = MULT;
    fSys[i][2].name = "SLACRELNORM1";

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

void SLACD_sh_30Filter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/SLACD/slac_d.data";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream nuclearfile("");
  nuclearfile << dataPath() << "rawdata/SLACD/nuclear_30/output/tables/group_result_table.csv";
  f2.open(nuclearfile.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << nuclearfile.str() << endl;
    exit(-1);
  }
  
  stringstream protonfile("");
  protonfile << dataPath() << "rawdata/SLACD/proton_30/output/tables/group_result_table.csv";
  f3.open(protonfile.str().c_str(), ios::in);
  
  if (f3.fail()) {
    cerr << "Error opening data file " << protonfile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double relnor;
  int nrep=200;
  int nrealsys=3;
  
  relnor = 1.1*0.5;    //relative normalisation of 1.1% between targets
  
  string line;

  getline(f2,line);
  getline(f3,line);
  
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    lstream >> fKin1[i];  //x
    lstream >> fKin2[i];  //q2
    lstream >> fData[i];  //obs 
    fKin3[i] = 0;
    
    //  SLAC gives errors in absolute value
    //  and we assume the sys errors to be uncorrelated in order
    //  to overestimate them and give less weight
    //  to these data points in evaluating chi2
   
    lstream >> fStat[i];
    
    lstream >> fSys[i][0].add;
    fSys[i][0].mult = fSys[i][0].add*100/fData[i];
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    
    fSys[i][1].mult = 1.7;  //absnorm
    fSys[i][1].add = fSys[i][1].mult*fData[i]*1e-2;
    fSys[i][1].type = MULT;
    fSys[i][1].name = "SLACNORM1";
    
    fSys[i][2].mult = -relnor;     //relnorm
    fSys[i][2].add = fSys[i][2].mult*fData[i]*1e-2;
    fSys[i][2].type = MULT;
    fSys[i][2].name = "SLACRELNORM1";

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
    cout << nuclear/proton_cv << "   " << 0.0 << endl;
    
  }
  
  f1.close();
  f2.close();
  f3.close();
}
