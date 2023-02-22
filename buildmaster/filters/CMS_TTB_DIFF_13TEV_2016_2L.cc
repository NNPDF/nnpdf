/*
Experiment: CERN-LHC-CMS (CMS)
Preprinted as CERN-EP-2018-039
Preprinted as CMS-TOP-17-002
Archived as: ARXIV:1803.08856
Published in Phys. Rev. D97 (2018) no 11, 112003
Measurement of differential cross sections for the production of top quark pairs
at sqrt(s)= 13 TeV

NOTE: for the normalised distributions, the covariance matrix has entries for
N-1 points. For example the TPT distribution has 6 bins, but the covmat is 5x5.
This is reflected in the meta file as the dataset having ndata = 5.

Also for some reason the covmat is NxNxN long (is repeated N times)

differential in the following variables are implemented:
NORMALISED:

1N) normalised top quark transverse momentum
2N) normalised top quark rapidity;
3N) normalised top pair invariant mass;
4N) normalised top pair rapidity;

UNNORMALISED:

1U) absolute top quark transverse momentum
2U) absolute top quark rapidity
3U) absolute top pair invariant mass
4U) absolute top pair rapidity

Raw data and covariance matrices are from HepData:
https://www.hepdata.net/record/ins1703993

1N) TABS 3-4 HepData; Fig 3, TAB 2 1811.06625
2N) TABS 43-44 HepData; Fig 12, TAB 2 1803.08856
3N) TABS 91-92 HepData; Fig 30, TAB 2 1803.08856
4N) TABS 83-84 HepData; Fig 27, TAB 2 1803.08856
*/

#include "CMS_TTB_DIFF_13TEV_2016_2L.h"

//N - NORMALISED distributions

//1N) Distribution differential in top quark transverse momentum
void  CMS_TTB_DIFF_13TEV_2016_2L_TPTNORMFilter::ReadData()
{
  //Opening files
  fstream f1, f2;

  //Data values
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/CMS_TTB_DIFF_13TEV_2016_2L_TPTNORM/HEPData-ins1703993-v1-Table_3.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/CMS_TTB_DIFF_13TEV_2016_2L_TPTNORM/HEPData-ins1703993-v1-Table_4.csv";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read central values - skip row since it doesn't have covmat entries
  string line;
  for(int i=0; i<11; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<fNData; i++)
    {
      double pt_top, dud;
      char comma;

      getline(f1,line);
      istringstream lstream(line);
      lstream >> pt_top >> comma
	      >> dud >> comma
	      >> dud >> comma
	      >> fData[i] >> comma
	      >> dud >> comma
	      >> dud >> comma
	      >> dud >> comma
	      >> dud;

      fKin1[i] = pt_top;  //pTt
      fKin2[i] = Mt*Mt;
      fKin3[i] = 13000;  //sqrt(s) [GeV]
      fStat[i] = 0.;
    }

  //Read covariance matrix - first 10 lines are header
  for(int i=0; i<10; i++)
    {
      getline(f2,line);
    }

  //Create covmat of correct dimensions
  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
    {
      covmat[i] = new double[fNData];

      for(int j=0; j<fNData; j++)
	{
	  double row, col;
	  char comma;

	  getline(f2,line);
	  istringstream lstream(line);
	  lstream >> row >> comma >> col >> comma  >> covmat[i][j];
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

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

// 2N) Distribution differential in top quark rapidity

void  CMS_TTB_DIFF_13TEV_2016_2L_TRAPNORMFilter::ReadData()
{
  //Opening files
  fstream f1, f2;

  //Data values
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/CMS_TTB_DIFF_13TEV_2016_2L_TRAPNORM/HEPData-ins1703993-v1-Table_43.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/CMS_TTB_DIFF_13TEV_2016_2L_TRAPNORM/HEPData-ins1703993-v1-Table_44.csv";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read central values - skip row since it doesn't have ovmat entries
  string line;
  for(int i=0; i<11; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<fNData; i++)
    {
      double y_top, dud;
      char comma;

      getline(f1,line);
      istringstream lstream(line);
      lstream >> y_top >> comma
	      >> dud >> comma
	      >> dud >> comma
	      >> fData[i] >> comma
	      >> dud >> comma
	      >> dud >> comma
	      >> dud >> comma
	      >> dud;

      fKin1[i] = y_top;  //y_t
      fKin2[i] = Mt*Mt;
      fKin3[i] = 13000;  //sqrt(s) [GeV]
      fStat[i] = 0.;
    }

  //Read covariance matrix - first 10 lines are header
  for(int i=0; i<10; i++)
    {
      getline(f2,line);
    }

  //Create covmat of correct dimensions
  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
    {
      covmat[i] = new double[fNData];

      for(int j=0; j<fNData; j++)
	{
	  double row, col;
	  char comma;

	  getline(f2,line);
	  istringstream lstream(line);
	  lstream >> row >> comma >> col >> comma  >> covmat[i][j];
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

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

// 3N) Distribution differential in top pair invariant mass

void  CMS_TTB_DIFF_13TEV_2016_2L_TTMNORMFilter::ReadData()
{
  //Opening files
  fstream f1, f2;

  //Data values
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/CMS_TTB_DIFF_13TEV_2016_2L_TTMNORM/HEPData-ins1703993-v1-Table_91.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/CMS_TTB_DIFF_13TEV_2016_2L_TTMNORM/HEPData-ins1703993-v1-Table_92.csv";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read central values - skip row since it doesn't have ovmat entries
  string line;
  for(int i=0; i<11; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<fNData; i++)
    {
      double m_tt, dud;
      char comma;

      getline(f1,line);
      istringstream lstream(line);
      lstream >> m_tt >> comma
	      >> dud >> comma
	      >> dud >> comma
	      >> fData[i] >> comma
	      >> dud >> comma
	      >> dud >> comma
	      >> dud >> comma
	      >> dud;

      fKin1[i] = m_tt;  //invariant mass tt
      fKin2[i] = Mt*Mt;
      fKin3[i] = 13000;  //sqrt(s) [GeV]
      fStat[i] = 0.;
    }

  //Read covariance matrix - first 10 lines are header
  for(int i=0; i<10; i++)
    {
      getline(f2,line);
    }

  //Create covmat of correct dimensions
  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
    {
      covmat[i] = new double[fNData];

      for(int j=0; j<fNData; j++)
	{
	  double row, col;
	  char comma;

	  getline(f2,line);
	  istringstream lstream(line);
	  lstream >> row >> comma >> col >> comma  >> covmat[i][j];
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

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

// 4N) Distribution differential in top pair rapidity

void  CMS_TTB_DIFF_13TEV_2016_2L_TTRAPNORMFilter::ReadData()
{
  //Opening files
  fstream f1, f2;

  //Data values
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/CMS_TTB_DIFF_13TEV_2016_2L_TTRAPNORM/HEPData-ins1703993-v1-Table_83.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/CMS_TTB_DIFF_13TEV_2016_2L_TTRAPNORM/HEPData-ins1703993-v1-Table_84.csv";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read central values - skip row since it doesn't have ovmat entries
  string line;
  for(int i=0; i<11; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<fNData; i++)
    {
      double y_tt, dud;
      char comma;

      getline(f1,line);
      istringstream lstream(line);
      lstream >> y_tt >> comma
	      >> dud >> comma
	      >> dud >> comma
	      >> fData[i] >> comma
	      >> dud >> comma
	      >> dud >> comma
	      >> dud >> comma
	      >> dud;

      fKin1[i] = y_tt;  // y tt
      fKin2[i] = Mt*Mt;
      fKin3[i] = 13000;  //sqrt(s) [GeV]
      fStat[i] = 0.;
    }

  //Read covariance matrix - first 10 lines are header
  for(int i=0; i<10; i++)
    {
      getline(f2,line);
    }

  //Create covmat of correct dimensions
  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
    {
      covmat[i] = new double[fNData];

      for(int j=0; j<fNData; j++)
	{
	  double row, col;
	  char comma;

	  getline(f2,line);
	  istringstream lstream(line);
	  lstream >> row >> comma >> col >> comma  >> covmat[i][j];
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

//U - UNNORMALISED distributions

//1U) Distribution differential in top quark transverse momentum
void  CMS_TTB_DIFF_13TEV_2016_2L_TPTFilter::ReadData()
{
  //Opening files
  fstream f1, f2;

  //Data values
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/CMS_TTB_DIFF_13TEV_2016_2L_TPT/HEPData-ins1703993-v1-Table_1.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/CMS_TTB_DIFF_13TEV_2016_2L_TPT/HEPData-ins1703993-v1-Table_2.csv";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read central values - skip row since it doesn't have covmat entries
  string line;
  for(int i=0; i<11; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<fNData; i++)
    {
      double pt_top, dud;
      char comma;

      getline(f1,line);
      istringstream lstream(line);
      lstream >> pt_top >> comma
	      >> dud >> comma
	      >> dud >> comma
	      >> fData[i] >> comma
	      >> dud >> comma
	      >> dud >> comma
	      >> dud >> comma
	      >> dud;

      fKin1[i] = pt_top;  //pTt
      fKin2[i] = Mt*Mt;
      fKin3[i] = 13000;  //sqrt(s) [GeV]
      fStat[i] = 0.;
    }

  //Read covariance matrix - first 10 lines are header
  for(int i=0; i<10; i++)
    {
      getline(f2,line);
    }

  //Create covmat of correct dimensions
  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
    {
      covmat[i] = new double[fNData];

      for(int j=0; j<fNData; j++)
	{
	  double row, col;
	  char comma;

	  getline(f2,line);
	  istringstream lstream(line);
	  lstream >> row >> comma >> col >> comma  >> covmat[i][j];
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

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

// 2U) Distribution differential in top quark rapidity

void  CMS_TTB_DIFF_13TEV_2016_2L_TRAPFilter::ReadData()
{
  //Opening files
  fstream f1, f2;

  //Data values
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/CMS_TTB_DIFF_13TEV_2016_2L_TRAP/HEPData-ins1703993-v1-Table_41.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/CMS_TTB_DIFF_13TEV_2016_2L_TRAP/HEPData-ins1703993-v1-Table_42.csv";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read central values - skip row since it doesn't have ovmat entries
  string line;
  for(int i=0; i<11; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<fNData; i++)
    {
      double y_top, dud;
      char comma;

      getline(f1,line);
      istringstream lstream(line);
      lstream >> y_top >> comma
	      >> dud >> comma
	      >> dud >> comma
	      >> fData[i] >> comma
	      >> dud >> comma
	      >> dud >> comma
	      >> dud >> comma
	      >> dud;

      fKin1[i] = y_top;  //y_t
      fKin2[i] = Mt*Mt;
      fKin3[i] = 13000;  //sqrt(s) [GeV]
      fStat[i] = 0.;
    }

  //Read covariance matrix - first 10 lines are header
  for(int i=0; i<10; i++)
    {
      getline(f2,line);
    }

  //Create covmat of correct dimensions
  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
    {
      covmat[i] = new double[fNData];

      for(int j=0; j<fNData; j++)
	{
	  double row, col;
	  char comma;

	  getline(f2,line);
	  istringstream lstream(line);
	  lstream >> row >> comma >> col >> comma  >> covmat[i][j];
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

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

// 3U) Distribution differential in top pair invariant mass

void  CMS_TTB_DIFF_13TEV_2016_2L_TTMFilter::ReadData()
{
  //Opening files
  fstream f1, f2;

  //Data values
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/CMS_TTB_DIFF_13TEV_2016_2L_TTM/HEPData-ins1703993-v1-Table_89.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/CMS_TTB_DIFF_13TEV_2016_2L_TTM/HEPData-ins1703993-v1-Table_90.csv";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read central values - skip row since it doesn't have ovmat entries
  string line;
  for(int i=0; i<11; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<fNData; i++)
    {
      double m_tt, dud;
      char comma;

      getline(f1,line);
      istringstream lstream(line);
      lstream >> m_tt >> comma
	      >> dud >> comma
	      >> dud >> comma
	      >> fData[i] >> comma
	      >> dud >> comma
	      >> dud >> comma
	      >> dud >> comma
	      >> dud;

      fKin1[i] = m_tt;  //invariant mass tt
      fKin2[i] = Mt*Mt;
      fKin3[i] = 13000;  //sqrt(s) [GeV]
      fStat[i] = 0.;
    }

  //Read covariance matrix - first 10 lines are header
  for(int i=0; i<10; i++)
    {
      getline(f2,line);
    }

  //Create covmat of correct dimensions
  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
    {
      covmat[i] = new double[fNData];

      for(int j=0; j<fNData; j++)
	{
	  double row, col;
	  char comma;

	  getline(f2,line);
	  istringstream lstream(line);
	  lstream >> row >> comma >> col >> comma  >> covmat[i][j];
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

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

// 4U) Distribution differential in top pair rapidity

void  CMS_TTB_DIFF_13TEV_2016_2L_TTRAPFilter::ReadData()
{
  //Opening files
  fstream f1, f2;

  //Data values
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/CMS_TTB_DIFF_13TEV_2016_2L_TTRAP/HEPData-ins1703993-v1-Table_81.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/CMS_TTB_DIFF_13TEV_2016_2L_TTRAP/HEPData-ins1703993-v1-Table_82.csv";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read central values - skip row since it doesn't have ovmat entries
  string line;
  for(int i=0; i<11; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<fNData; i++)
    {
      double y_tt, dud;
      char comma;

      getline(f1,line);
      istringstream lstream(line);
      lstream >> y_tt >> comma
	      >> dud >> comma
	      >> dud >> comma
	      >> fData[i] >> comma
	      >> dud >> comma
	      >> dud >> comma
	      >> dud >> comma
	      >> dud;

      fKin1[i] = y_tt;  // y tt
      fKin2[i] = Mt*Mt;
      fKin3[i] = 13000;  //sqrt(s) [GeV]
      fStat[i] = 0.;
    }

  //Read covariance matrix - first 10 lines are header
  for(int i=0; i<10; i++)
    {
      getline(f2,line);
    }

  //Create covmat of correct dimensions
  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
    {
      covmat[i] = new double[fNData];

      for(int j=0; j<fNData; j++)
	{
	  double row, col;
	  char comma;

	  getline(f2,line);
	  istringstream lstream(line);
	  lstream >> row >> comma >> col >> comma  >> covmat[i][j];
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
