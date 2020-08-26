/*
Experiment: CERN-LHC-CMS (CMS)
Preprinted as CERN-EP-2018-039
Preprinted as CMS-TOP-17-002
Archived as: ARXIV:1803.08856
Published in Phys. Rev. D97 (2018) no 11, 112003
Measurement of differential cross sections for the production of top quark pairs
and of additional jets in lepton+jets events from pp collisions at sqrt(s)= 13
TeV
differential in the following variables are implemented:
1U) unnormalised top quark transverse momentum
2U) unnormalised top quark rapidity;
3U) unnormalised top pair invariant mass;
4U) unnormalised top pair rapidity;
1N) normalised top quark transverse momentum
2N) normalised top quark rapidity;
3N) normalised top pair invariant mass;
4N) normalised top pair rapidity;
Raw data and covariance matrices are from HepData:
https://www.hepdata.net/record/ins1663958
1U) TABS 182-183 HepData; Fig 11, TAB 3 1803.08856
2U) TABS 184-185 HepData; Fig 12, TAB 3 1803.08856
3U) TABS 186-187 HepData; Fig 16, TAB 5 1803.08856
4U) TABS 190-191 HepData; Fig 16, TAB 3 1803.08856
1N) TABS 215-216 HepData; Fig 11, TAB 4 1803.08856
2N) TABS 217-218 HepData; Fig 12, TAB 4 1803.08856
3N) TABS 219-220 HepData; Fig 16, TAB 6 1803.08856
4N) TABS 223-224 HepData; Fig 16, TAB 4 1803.08856
NOTES:
The Normalised distributions are normalised to themselves and so the final bin
is a linear combination of the other bins. For the normalised distributions here
the meta file was changed to set NData: N_bins - 1 and then the loop which reads
the covmat was edited to loop over i for i<(N_bins - 1) and for j<N_bins and then
only read the line into the covmat if j<(N_bins - 1). That way the rawdata lines
which correspond to j==N_bins are skipped.
*/

#include "CMS_TTB_DIFF_13TEV_2016_LJ.h"

//U - UNNORMALISED distributions

//1U) Distribution differential in top quark transverse momentum
void CMS_TTB_DIFF_13TEV_2016_LJ_TPTFilter::ReadData()
{
  //Opening files
  fstream f1, f2;

  //Data values
  stringstream datafile("");
  datafile << dataPath()
           << "rawdata/CMS_TTB_DIFF_13TEV_2016_LJ_TPT/HEPData-ins1663958-v2-Table_182.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  //Covariance matrix
  stringstream covfile("");
  covfile << dataPath()
          << "rawdata/CMS_TTB_DIFF_13TEV_2016_LJ_TPT/HEPData-ins1663958-v2-Table_183.csv";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  //Read central values
  string line;
  for (int i = 0; i < 11; i++)
  {
    getline(f1, line);
  }

  for (int i = 0; i < fNData; i++)
  {
    double pt_top, dud;
    char comma;

    getline(f1, line);
    istringstream lstream(line);
    lstream >> pt_top >> comma
            >> dud >> comma
            >> dud >> comma
            >> fData[i] >> comma
            >> dud >> comma
            >> dud >> comma
            >> dud >> comma
            >> dud;

    fKin1[i] = pt_top; //pTt
    fKin2[i] = Mt * Mt;
    fKin3[i] = 13000; //sqrt(s) [GeV]
    fStat[i] = 0.;
  }

  //Read covariance matrix
  for (int i = 0; i < 11; i++)
  {
    getline(f2, line);
  }

  //Create covmat of correct dimensions
  double **covmat = new double *[fNData];
  for (int i = 0; i < fNData; i++)
  {
    covmat[i] = new double[fNData];

    for (int j = 0; j < fNData; j++)
    {
      double row, col;
      char comma;

      getline(f2, line);
      istringstream lstream(line);
      lstream >> row >> comma >> col >> comma >> covmat[i][j];
    }
  }

  //Generate artificial systematics
  double **syscor = new double *[fNData];
  for (int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if (!genArtSys(fNData, covmat, syscor))
  {
    throw runtime_error("Couldn't generate artificial systematics for " + fSetName);
  }

  for (int i = 0; i < fNData; i++)
  {
    for (int j = 0; j < fNData; j++)
    {
      fSys[i][j].add = syscor[i][j];
      fSys[i][j].mult = fSys[i][j].add * 1e2 / fData[i];
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

//2U) Distribution differential in top quark rapidity
void CMS_TTB_DIFF_13TEV_2016_LJ_TRAPFilter::ReadData()
{
  //Opening files
  fstream f1, f2;

  //Data values
  stringstream datafile("");
  datafile << dataPath()
           << "rawdata/CMS_TTB_DIFF_13TEV_2016_LJ_TRAP/HEPData-ins1663958-v2-Table_184.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  //Covariance matrix
  stringstream covfile("");
  covfile << dataPath()
          << "rawdata/CMS_TTB_DIFF_13TEV_2016_LJ_TRAP/HEPData-ins1663958-v2-Table_185.csv";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  //Read central values
  string line;
  for (int i = 0; i < 11; i++)
  {
    getline(f1, line);
  }

  for (int i = 0; i < fNData; i++)
  {
    double y_top, dud;
    char comma;

    getline(f1, line);
    istringstream lstream(line);
    lstream >> y_top >> comma
            >> dud >> comma
            >> dud >> comma
            >> fData[i] >> comma
            >> dud >> comma
            >> dud >> comma
            >> dud >> comma
            >> dud;

    fKin1[i] = y_top; //top rapidity
    fKin2[i] = Mt * Mt;
    fKin3[i] = 13000; //sqrt(s) [GeV]
    fStat[i] = 0.;
  }

  //Read covariance matrix
  for (int i = 0; i < 11; i++)
  {
    getline(f2, line);
  }

  //Create covmat of correct dimensions
  double **covmat = new double *[fNData];
  for (int i = 0; i < fNData; i++)
  {
    covmat[i] = new double[fNData];

    for (int j = 0; j < fNData; j++)
    {
      double row, col;
      char comma;

      getline(f2, line);
      istringstream lstream(line);
      lstream >> row >> comma >> col >> comma >> covmat[i][j];
    }
  }

  //Generate artificial systematics
  double **syscor = new double *[fNData];
  for (int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if (!genArtSys(fNData, covmat, syscor))
  {
    throw runtime_error("Couldn't generate artificial systematics for " + fSetName);
  }

  for (int i = 0; i < fNData; i++)
  {
    for (int j = 0; j < fNData; j++)
    {
      fSys[i][j].add = syscor[i][j];
      fSys[i][j].mult = fSys[i][j].add * 1e2 / fData[i];
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

//3U) Distribution differential in top pair invariant mass

void CMS_TTB_DIFF_13TEV_2016_LJ_TTMFilter::ReadData()
{
  //Opening files
  fstream f1, f2;

  //Data values
  stringstream datafile("");
  datafile << dataPath()
           << "rawdata/CMS_TTB_DIFF_13TEV_2016_LJ_TTM/HEPData-ins1663958-v2-Table_186.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  //Covariance matrix
  stringstream covfile("");
  covfile << dataPath()
          << "rawdata/CMS_TTB_DIFF_13TEV_2016_LJ_TTM/HEPData-ins1663958-v2-Table_187.csv";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  //Read central values
  string line;
  for (int i = 0; i < 11; i++)
  {
    getline(f1, line);
  }

  for (int i = 0; i < fNData; i++)
  {
    double m_tt, dud;
    char comma;

    getline(f1, line);
    istringstream lstream(line);
    lstream >> m_tt >> comma
            >> dud >> comma
            >> dud >> comma
            >> fData[i] >> comma
            >> dud >> comma
            >> dud >> comma
            >> dud >> comma
            >> dud;

    fKin1[i] = m_tt; //M tt
    fKin2[i] = Mt * Mt;
    fKin3[i] = 13000; //sqrt(s) [GeV]
    fStat[i] = 0.;
  }

  //Read covariance matrix
  for (int i = 0; i < 11; i++)
  {
    getline(f2, line);
  }

  //Create covmat of correct dimensions
  double **covmat = new double *[fNData];
  for (int i = 0; i < fNData; i++)
  {
    covmat[i] = new double[fNData];

    for (int j = 0; j < fNData; j++)
    {
      double row, col;
      char comma;

      getline(f2, line);
      istringstream lstream(line);
      lstream >> row >> comma >> col >> comma >> covmat[i][j];
    }
  }

  //Generate artificial systematics
  double **syscor = new double *[fNData];
  for (int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if (!genArtSys(fNData, covmat, syscor))
  {
    throw runtime_error("Couldn't generate artificial systematics for " + fSetName);
  }

  for (int i = 0; i < fNData; i++)
  {
    for (int j = 0; j < fNData; j++)
    {
      fSys[i][j].add = syscor[i][j];
      fSys[i][j].mult = fSys[i][j].add * 1e2 / fData[i];
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

//4U) Distribution differential in top pair rapidity

void CMS_TTB_DIFF_13TEV_2016_LJ_TTRAPFilter::ReadData()
{
  //Opening files
  fstream f1, f2;

  //Data values
  stringstream datafile("");
  datafile << dataPath()
           << "rawdata/CMS_TTB_DIFF_13TEV_2016_LJ_TTRAP/HEPData-ins1663958-v2-Table_190.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  //Covariance matrix
  stringstream covfile("");
  covfile << dataPath()
          << "rawdata/CMS_TTB_DIFF_13TEV_2016_LJ_TTRAP/HEPData-ins1663958-v2-Table_191.csv";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  //Read central values
  string line;
  for (int i = 0; i < 11; i++)
  {
    getline(f1, line);
  }

  for (int i = 0; i < fNData; i++)
  {
    double y_tt, dud;
    char comma;

    getline(f1, line);
    istringstream lstream(line);
    lstream >> y_tt >> comma
            >> dud >> comma
            >> dud >> comma
            >> fData[i] >> comma
            >> dud >> comma
            >> dud >> comma
            >> dud >> comma
            >> dud;

    fKin1[i] = y_tt; //rapidity tt
    fKin2[i] = Mt * Mt;
    fKin3[i] = 13000; //sqrt(s) [GeV]
    fStat[i] = 0.;
  }

  //Read covariance matrix
  for (int i = 0; i < 11; i++)
  {
    getline(f2, line);
  }

  //Create covmat of correct dimensions
  double **covmat = new double *[fNData];
  for (int i = 0; i < fNData; i++)
  {
    covmat[i] = new double[fNData];

    for (int j = 0; j < fNData; j++)
    {
      double row, col;
      char comma;

      getline(f2, line);
      istringstream lstream(line);
      lstream >> row >> comma >> col >> comma >> covmat[i][j];
    }
  }

  //Generate artificial systematics
  double **syscor = new double *[fNData];
  for (int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if (!genArtSys(fNData, covmat, syscor))
  {
    throw runtime_error("Couldn't generate artificial systematics for " + fSetName);
  }

  for (int i = 0; i < fNData; i++)
  {
    for (int j = 0; j < fNData; j++)
    {
      fSys[i][j].add = syscor[i][j];
      fSys[i][j].mult = fSys[i][j].add * 1e2 / fData[i];
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

// N - NORMALISED distributions

//1N) Distribution differential in top transverse momentum

void CMS_TTB_DIFF_13TEV_2016_LJ_TPTNORMFilter::ReadData()
{
  //Opening files
  fstream f1, f2;

  //Data values
  stringstream datafile("");
  datafile << dataPath()
           << "rawdata/CMS_TTB_DIFF_13TEV_2016_LJ_TPTNORM/HEPData-ins1663958-v2-Table_215.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  //Covariance matrix
  stringstream covfile("");
  covfile << dataPath()
          << "rawdata/CMS_TTB_DIFF_13TEV_2016_LJ_TPTNORM/HEPData-ins1663958-v2-Table_216.csv";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  //Read central values
  string line;
  for (int i = 0; i < 11; i++)
  {
    getline(f1, line);
  }

  for (int i = 0; i < fNData; i++)
  {
    double pt_top, dud;
    char comma;

    getline(f1, line);
    istringstream lstream(line);
    lstream >> pt_top >> comma
            >> dud >> comma
            >> dud >> comma
            >> fData[i] >> comma
            >> dud >> comma
            >> dud >> comma
            >> dud >> comma
            >> dud;

    fKin1[i] = pt_top; //pTt
    fKin2[i] = Mt * Mt;
    fKin3[i] = 13000; //sqrt(s) [GeV]
    fStat[i] = 0.;
  }

  //Read covariance matrix
  for (int i = 0; i < 11; i++)
  {
    getline(f2, line);
  }

  //Create covmat of correct dimensions
  // NOTE: fNData is set to N_bins - 1 so we need to loop with j<fNData+1 but
  // only fill in covmat if j<fNData so last bins are skipped.
  double **covmat = new double *[fNData];
  for (int i = 0; i < fNData; i++)
  {
    covmat[i] = new double[fNData];

    for (int j = 0; j < (fNData + 1); j++)
    {
      double row, col;
      char comma;

      getline(f2, line);
      if (j < fNData)
      {
        istringstream lstream(line);
        lstream >> row >> comma >> col >> comma >> covmat[i][j];
      }
    }
  }

  //Generate artificial systematics
  double **syscor = new double *[fNData];
  for (int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if (!genArtSys(fNData, covmat, syscor))
  {
    throw runtime_error("Couldn't generate artificial systematics for " + fSetName);
  }

  for (int i = 0; i < fNData; i++)
  {
    for (int j = 0; j < fNData; j++)
    {
      fSys[i][j].add = syscor[i][j];
      fSys[i][j].mult = fSys[i][j].add * 1e2 / fData[i];
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

//2N) Distribution differential in top rapidity

void CMS_TTB_DIFF_13TEV_2016_LJ_TRAPNORMFilter::ReadData()
{
  //Opening files
  fstream f1, f2;

  //Data values
  stringstream datafile("");
  datafile << dataPath()
           << "rawdata/CMS_TTB_DIFF_13TEV_2016_LJ_TRAPNORM/HEPData-ins1663958-v2-Table_217.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  //Covariance matrix
  stringstream covfile("");
  covfile << dataPath()
          << "rawdata/CMS_TTB_DIFF_13TEV_2016_LJ_TRAPNORM/HEPData-ins1663958-v2-Table_218.csv";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  //Read central values
  string line;
  for (int i = 0; i < 11; i++)
  {
    getline(f1, line);
  }

  for (int i = 0; i < fNData; i++)
  {
    double y_top, dud;
    char comma;

    getline(f1, line);
    istringstream lstream(line);
    lstream >> y_top >> comma
            >> dud >> comma
            >> dud >> comma
            >> fData[i] >> comma
            >> dud >> comma
            >> dud >> comma
            >> dud >> comma
            >> dud;

    fKin1[i] = y_top; //top rapidity
    fKin2[i] = Mt * Mt;
    fKin3[i] = 13000; //sqrt(s) [GeV]
    fStat[i] = 0.;
  }

  //Read covariance matrix
  for (int i = 0; i < 11; i++)
  {
    getline(f2, line);
  }

  //Create covmat of correct dimensions
  // NOTE: fNData is set to N_bins - 1 so we need to loop with j<fNData+1 but
  // only fill in covmat if j<fNData so last bins are skipped.
  double **covmat = new double *[fNData];
  for (int i = 0; i < fNData; i++)
  {
    covmat[i] = new double[fNData];

    for (int j = 0; j < (fNData + 1); j++)
    {
      double row, col;
      char comma;

      getline(f2, line);
      if (j < fNData)
      {
        istringstream lstream(line);
        lstream >> row >> comma >> col >> comma >> covmat[i][j];
      }
    }
  }

  //Generate artificial systematics
  double **syscor = new double *[fNData];
  for (int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if (!genArtSys(fNData, covmat, syscor))
  {
    throw runtime_error("Couldn't generate artificial systematics for " + fSetName);
  }

  for (int i = 0; i < fNData; i++)
  {
    for (int j = 0; j < fNData; j++)
    {
      fSys[i][j].add = syscor[i][j];
      fSys[i][j].mult = fSys[i][j].add * 1e2 / fData[i];
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

//3N) Distribution differential in top pair invariant mass

void CMS_TTB_DIFF_13TEV_2016_LJ_TTMNORMFilter::ReadData()
{
  //Opening files
  fstream f1, f2;

  //Data values
  stringstream datafile("");
  datafile << dataPath()
           << "rawdata/CMS_TTB_DIFF_13TEV_2016_LJ_TTMNORM/HEPData-ins1663958-v2-Table_219.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  //Covariance matrix
  stringstream covfile("");
  covfile << dataPath()
          << "rawdata/CMS_TTB_DIFF_13TEV_2016_LJ_TTMNORM/HEPData-ins1663958-v2-Table_220.csv";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  //Read central values
  string line;
  for (int i = 0; i < 11; i++)
  {
    getline(f1, line);
  }

  for (int i = 0; i < fNData; i++)
  {
    double m_tt, dud;
    char comma;

    getline(f1, line);
    istringstream lstream(line);
    lstream >> m_tt >> comma
            >> dud >> comma
            >> dud >> comma
            >> fData[i] >> comma
            >> dud >> comma
            >> dud >> comma
            >> dud >> comma
            >> dud;

    fKin1[i] = m_tt; //M tt
    fKin2[i] = Mt * Mt;
    fKin3[i] = 13000; //sqrt(s) [GeV]
    fStat[i] = 0.;
  }

  //Read covariance matrix
  for (int i = 0; i < 11; i++)
  {
    getline(f2, line);
  }

  //Create covmat of correct dimensions
  // NOTE: fNData is set to N_bins - 1 so we need to loop with j<fNData+1 but
  // only fill in covmat if j<fNData so last bins are skipped.
  double **covmat = new double *[fNData];
  for (int i = 0; i < fNData; i++)
  {
    covmat[i] = new double[fNData];

    for (int j = 0; j < (fNData + 1); j++)
    {
      double row, col;
      char comma;

      getline(f2, line);
      if (j < fNData)
      {
        istringstream lstream(line);
        lstream >> row >> comma >> col >> comma >> covmat[i][j];
      }
    }
  }

  //Generate artificial systematics
  double **syscor = new double *[fNData];
  for (int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if (!genArtSys(fNData, covmat, syscor))
  {
    throw runtime_error("Couldn't generate artificial systematics for " + fSetName);
  }

  for (int i = 0; i < fNData; i++)
  {
    for (int j = 0; j < fNData; j++)
    {
      fSys[i][j].add = syscor[i][j];
      fSys[i][j].mult = fSys[i][j].add * 1e2 / fData[i];
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

//4N) Distribution differential in top pair rapidity

void CMS_TTB_DIFF_13TEV_2016_LJ_TTRAPNORMFilter::ReadData()
{
  //Opening files
  fstream f1, f2;

  //Data values
  stringstream datafile("");
  datafile << dataPath()
           << "rawdata/CMS_TTB_DIFF_13TEV_2016_LJ_TTRAPNORM/HEPData-ins1663958-v2-Table_223.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  //Covariance matrix
  stringstream covfile("");
  covfile << dataPath()
          << "rawdata/CMS_TTB_DIFF_13TEV_2016_LJ_TTRAPNORM/HEPData-ins1663958-v2-Table_224.csv";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  //Read central values
  string line;
  for (int i = 0; i < 11; i++)
  {
    getline(f1, line);
  }

  for (int i = 0; i < fNData; i++)
  {
    double y_tt, dud;
    char comma;

    getline(f1, line);
    istringstream lstream(line);
    lstream >> y_tt >> comma
            >> dud >> comma
            >> dud >> comma
            >> fData[i] >> comma
            >> dud >> comma
            >> dud >> comma
            >> dud >> comma
            >> dud;

    fKin1[i] = y_tt; //rapidity tt
    fKin2[i] = Mt * Mt;
    fKin3[i] = 13000; //sqrt(s) [GeV]
    fStat[i] = 0.;
  }

  //Read covariance matrix
  for (int i = 0; i < 11; i++)
  {
    getline(f2, line);
  }

  //Create covmat of correct dimensions
  // NOTE: fNData is set to N_bins - 1 so we need to loop with j<fNData+1 but
  // only fill in covmat if j<fNData so last bins are skipped.
  double **covmat = new double *[fNData];
  for (int i = 0; i < fNData; i++)
  {
    covmat[i] = new double[fNData];

    for (int j = 0; j < (fNData + 1); j++)
    {
      double row, col;
      char comma;

      getline(f2, line);
      if (j < fNData)
      {
        istringstream lstream(line);
        lstream >> row >> comma >> col >> comma >> covmat[i][j];
      }
    }
  }

  //Generate artificial systematics
  double **syscor = new double *[fNData];
  for (int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if (!genArtSys(fNData, covmat, syscor))
  {
    throw runtime_error("Couldn't generate artificial systematics for " + fSetName);
  }

  for (int i = 0; i < fNData; i++)
  {
    for (int j = 0; j < fNData; j++)
    {
      fSys[i][j].add = syscor[i][j];
      fSys[i][j].mult = fSys[i][j].add * 1e2 / fData[i];
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
