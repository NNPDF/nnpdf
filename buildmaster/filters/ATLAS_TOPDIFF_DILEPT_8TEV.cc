/*
Experiment: CERN-LHC-ATLAS (ATLAS)
Archived as: ARXIV:1607.07281
Published in Phys.Rev. D94 (2016) no.9, 092003

Record in: INSPIRE
Record in: HEPData https://www.hepdata.net/record/ins1477814

Description of the measurement 
Measurements of normalized differential cross-sections of top quark pair (tt¯) 
production are presented as a function of the mass, the transverse momentum and 
the rapidity of the tt¯ system in proton-proton collisions at center-of-mass 
energies of s√ = 8 TeV. The dataset corresponds to an integrated luminosity of 
20.2 fb−1 at 8 TeV, recorded with the ATLAS detector at the Large Hadron 
Collider. Events with top quark pair signatures are selected in the dilepton 
final state, requiring exactly two charged leptons and at least two jets with 
at least one of the jets identified as likely to contain a b-hadron. 
The measured distributions are corrected for detector effects and selection 
efficiency to cross-sections at the parton level.

Description of the buildmaster implementation
Normalized and unnormalised cross sections for the distributions 
(dilepton channel)
differential in the following variables are implemented:
1) top quark pair invariant mass;                   
2) top quark pair rapidity.       

Raw data and full breakdown of systematic uncertainties are from HepData:
https://www.hepdata.net/record/ins1477814
N: normalised
U: unnormalised

1N) TABS 4, 16 HepData; TAB 5 in the paper
2N) TABS 6, 18 HepData; TAB 5 in the paper
1U) TABS 10, 22 HepData
2U) TABS 12, 24 HepData

*/
 
#include "ATLAS_TOPDIFF_DILEPT_8TEV.h"

//N - NORMALISED distributions

//1) Distribution differential in top quark pair invariant mass
void  ATLAS_TOPDIFF_DILEPT_8TEV_TTMNORMFilter::ReadData()
{
  //Opening files
  fstream f1, f2;

  //Central values and statistical uncertainties
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/ATLAS_TOPDIFF_DILEPT_8TEV_TTMNORM/tab4.txt";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/ATLAS_TOPDIFF_DILEPT_8TEV_TTMNORM/tab16.txt";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail()) 
    {
      cerr << "Error opening covariance matrix file " << datafile.str() << endl;
      exit(-1);
    }
  

  //Read central values 
  string line;
  for(int i=0; i<8; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      double adum;
      double mtt;
      lstream >> mtt;
      lstream >> adum;
      lstream >> adum;
      
      fKin1[i] = mtt;   
      fKin2[i] = Mt*Mt;       
      fKin3[i] = 8000;     //sqrt(s)
      fStat[i] = 0;
      
      lstream >> fData[i]; //normalized differential distribution
      fData[i] *= 1e-3;

      for(int j=0; j<4; j++)
	{
	  lstream >> adum;
	}
    }

  //Read covariance matrix
  for(int i=0; i<9; i++)
    {
      getline(f2,line);
    }

  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
  {
    covmat[i] = new double[fNData];
    getline(f2,line);
    istringstream lstream(line);
    double adum;
    lstream >> adum >> adum >> adum;
    for(int j=0; j<fNData; j++)
    {
      lstream >> covmat[i][j];
      covmat[i][j] *= 1e-6;
    }
  }


  //Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
    {
      throw runtime_error("Couldn't generate artificial systematics for " + fSetName);
    }

  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  fSys[i][j].add  = syscor[i][j];
	  fSys[i][j].mult = fSys[i][j].add*1e2/fData[i];
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = "CORR";
	}
    }


  // Clean-up
  for (int i=0; i<fNData; i++)
    delete[] syscor[i];

  delete[] syscor;

  for(int i=0; i<fNData; i++)
    delete[] covmat[i];

  delete[] covmat;


  f1.close();
  f2.close();

}

// -------------------------------------------------------------------

//2) Distribution differential in top quark pair rapidity
void  ATLAS_TOPDIFF_DILEPT_8TEV_TTRAPNORMFilter::ReadData()
{
  //Opening files
  fstream f1, f2;

  //Central values and statistical uncertainties
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/ATLAS_TOPDIFF_DILEPT_8TEV_TTRAPNORM/tab6.txt";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/ATLAS_TOPDIFF_DILEPT_8TEV_TTRAPNORM/tab18.txt";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail()) 
    {
      cerr << "Error opening covariance matrix file " << datafile.str() << endl;
      exit(-1);
    }
  

  //Read central values 
  string line;
  for(int i=0; i<8; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      double adum;
      double raptt;
      lstream >> raptt;
      lstream >> adum;
      lstream >> adum;
      
      fKin1[i] = raptt;   
      fKin2[i] = Mt*Mt;       
      fKin3[i] = 8000;     //sqrt(s)
      fStat[i] = 0;
      
      lstream >> fData[i]; //normalized differential distribution

      for(int j=0; j<4; j++)
	{
	  lstream >> adum;
	}
    }

  //Read covariance matrix
  for(int i=0; i<9; i++)
    {
      getline(f2,line);
    }

  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
  {
    covmat[i] = new double[fNData];
    getline(f2,line);
    istringstream lstream(line);
    double adum;
    lstream >> adum >> adum >> adum;
    for(int j=0; j<fNData; j++)
    {
      lstream >> covmat[i][j];
    }
  }


  //Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
    {
      throw runtime_error("Couldn't generate artificial systematics for " + fSetName);
    }

  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  fSys[i][j].add  = syscor[i][j];
	  fSys[i][j].mult  = fSys[i][j].add*1e2/fData[i];
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = "CORR";
	}
    }


  // Clean-up
  for (int i=0; i<fNData; i++)
    delete[] syscor[i];

  delete[] syscor;

  for(int i=0; i<fNData; i++)
    delete[] covmat[i];

  delete[] covmat;


  f1.close();
  f2.close();

}

// -------------------------------------------------------------------

//U - UNNORMALISED distributions

//1) Distribution differential in top quark pair invariant mass
void  ATLAS_TOPDIFF_DILEPT_8TEV_TTMFilter::ReadData()
{
  //Opening files
  fstream f1, f2;

  //Central values and statistical uncertainties
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/ATLAS_TOPDIFF_DILEPT_8TEV_TTM/tab10.txt";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/ATLAS_TOPDIFF_DILEPT_8TEV_TTM/tab22.txt";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail()) 
    {
      cerr << "Error opening covariance matrix file " << datafile.str() << endl;
      exit(-1);
    }
  

  //Read central values 
  string line;
  for(int i=0; i<8; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      double adum;
      double mtt;
      lstream >> mtt;
      lstream >> adum;
      lstream >> adum;
      
      fKin1[i] = mtt;   
      fKin2[i] = Mt*Mt;       
      fKin3[i] = 8000;     //sqrt(s)
      fStat[i] = 0;
      
      lstream >> fData[i]; 

      for(int j=0; j<4; j++)
	{
	  lstream >> adum;
	}
    }

  //Read covariance matrix
  for(int i=0; i<9; i++)
    {
      getline(f2,line);
    }

  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
  {
    covmat[i] = new double[fNData];
    getline(f2,line);
    istringstream lstream(line);
    double adum;
    lstream >> adum >> adum >> adum;
    for(int j=0; j<fNData; j++)
    {
      lstream >> covmat[i][j];
    }
  }


  //Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
    {
      throw runtime_error("Couldn't generate artificial systematics for " + fSetName);
    }

  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  fSys[i][j].add  = syscor[i][j];
	  fSys[i][j].mult = fSys[i][j].add*1e2/fData[i];
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = "CORR";
	}

      //Luminosity uncertainty
      fSys[i][6].mult = 1.9;  //%
      fSys[i][6].add  = fSys[i][6].mult/100.*fData[i];
      fSys[i][6].type = MULT;
      fSys[i][6].name = "ATLASLUMI12";

    }

  // Clean-up
  for (int i=0; i<fNData; i++)
    delete[] syscor[i];

  delete[] syscor;

  for(int i=0; i<fNData; i++)
    delete[] covmat[i];

  delete[] covmat;


  f1.close();
  f2.close();

}

// -------------------------------------------------------------------

//2) Distribution differential in top quark pair rapidity
void  ATLAS_TOPDIFF_DILEPT_8TEV_TTRAPFilter::ReadData()
{
  //Opening files
  fstream f1, f2;

  //Central values and statistical uncertainties
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/ATLAS_TOPDIFF_DILEPT_8TEV_TTRAP/tab12.txt";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/ATLAS_TOPDIFF_DILEPT_8TEV_TTRAP/tab24.txt";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail()) 
    {
      cerr << "Error opening covariance matrix file " << datafile.str() << endl;
      exit(-1);
    }
  

  //Read central values 
  string line;
  for(int i=0; i<8; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      double adum;
      double raptt;
      lstream >> raptt;
      lstream >> adum;
      lstream >> adum;
      
      fKin1[i] = raptt;   
      fKin2[i] = Mt*Mt;       
      fKin3[i] = 8000;     //sqrt(s)
      fStat[i] = 0;
      
      lstream >> fData[i]; 

      for(int j=0; j<4; j++)
	{
	  lstream >> adum;
	}
    }

  //Read covariance matrix
  for(int i=0; i<9; i++)
    {
      getline(f2,line);
    }

  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
  {
    covmat[i] = new double[fNData];
    getline(f2,line);
    istringstream lstream(line);
    double adum;
    lstream >> adum >> adum >> adum;
    for(int j=0; j<fNData; j++)
    {
      lstream >> covmat[i][j];
    }
  }


  //Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
    {
      throw runtime_error("Couldn't generate artificial systematics for " + fSetName);
    }

  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  fSys[i][j].add  = syscor[i][j];
	  fSys[i][j].mult  = fSys[i][j].add*1e2/fData[i];
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = "CORR";
	}

      //Luminosity uncertainty
      fSys[i][5].mult = 1.9;  //%
      fSys[i][5].add  = fSys[i][5].mult/100.*fData[i];
      fSys[i][5].type = MULT;
      fSys[i][5].name = "ATLASLUMI12";

    }


  // Clean-up
  for (int i=0; i<fNData; i++)
    delete[] syscor[i];

  delete[] syscor;

  for(int i=0; i<fNData; i++)
    delete[] covmat[i];

  delete[] covmat;


  f1.close();
  f2.close();

}

// -------------------------------------------------------------------
