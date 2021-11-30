/*
Differential Drell-Yan cross section measurements in the dielectron, dimuon and dilepton invariant
masses from CMS at the 13 TeV LHC
LHC-CMS 13 TeV
-------------

Integrated luminosity: 5.1 fb^-1

Archived as: https://arxiv.org/abs/1812.10529v2
Published in: https://link.springer.com/article/10.1007%2FJHEP12%282019%29061
HEPData: https://www.hepdata.net/record/ins1711625

Distributions are in pb/GeV.

Notes:
- Uncertainties are included via a covariance matrix, which is provided for each dataset
- For the combined channel, the covariance matrix includes all sources of uncertainty (including
  the statistical and luminosity uncertainties!)
- For each of the electron and muon channels, the covariance matrix includes all sources of
  uncertainty, including the statistical uncertainty, but EXCLUDING the luminosity uncertainty.
  This is therefore accounted for explicitly via the corresponding filters.
*/


/*
This is in the combined electron-muon channel
*/
#include "CMS_HMDY_13TEV.h"

void CMS_HMDY_13TEVFilter::ReadData()
{
  // Opening files
  fstream file1, file2;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << "/" << fSetName << ".data";
  file1.open(datafile.str().c_str(), ios::in);

  if (file1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream covfile("");
  covfile << dataPath() << "rawdata/" << fSetName << "/" << fSetName << ".cov";
  file2.open(covfile.str().c_str(), ios::in);

  if (file2.fail()) {
    cerr << "Error opening data file " << covfile.str() << endl;
    exit(-1);
  }

  string line;

  // Skip comments at top of file
  for (int i = 0; i < 7; i++)
  {
    getline(file1, line);
  }

  double Mll, Mll_low, Mll_high;

  // Filter data file
  for (int i = 0; i < fNData; i++)
  {
    getline(file1, line);
    istringstream lstream(line);

    lstream >> Mll_low >> Mll_high >> fData[i];
    Mll = (Mll_low + Mll_high)/2;
    fKin1[i] = Mll;                               // Mll
    fKin2[i] = pow(Mll,2);                        // Mll^2
    fKin3[i] = 13000;                             // sqrt(s) (GeV)
    fStat[i] = 0;                                 // Statistical uncertainty is included in covariance matrix so set this to 0
  }

  // Skip comments at top of file
  for (int i = 0; i < 4; i++)
  {
    getline(file2, line);
  }

  // Filter covariance matrix file
  double** covmat = new double*[fNData];
  for (int i = 0; i < fNData; i++)
  {
    covmat[i] = new double[fNData];
    getline(file2, line);
    istringstream lstream(line);
    for (int j = 0; j < fNData; j++)
    {
      lstream >> covmat[i][j];
    }
  }

  // Generate artificial systematics
  double** syscor = new double*[fNData];
  for (int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if (!genArtSys(fNData,covmat,syscor))
   {
     cerr << " in " << fSetName << endl;
     exit(-1);
   }

  // Assign artificial systematics
  for (int i = 0; i < fNData; i++)
  {
    for (int j = 0; j < fNSys; j++)
    {
      fSys[i][j].add = syscor[i][j];
      fSys[i][j].mult = fSys[i][j].add*100/fData[i];
      fSys[i][j].type = ADD;
      fSys[i][j].name = "CORR";
    }
  }

  file1.close();
  file2.close();
  
  // Clean-up
  for (int i=0; i<fNData; i++)
    delete[] syscor[i];
  
  delete[] syscor;
  
  for (int i=0; i<fNData; i++)
    delete[] covmat[i];
  
  delete[] covmat;
  
}

/*
Same as CMS_HMDY_13TEV but in the electron channel
*/
void CMS_HMDY_DE_13TEVFilter::ReadData()
{
  // Opening files
  fstream file1, file2;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << "/" << fSetName << ".data";
  file1.open(datafile.str().c_str(), ios::in);

  if (file1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream covfile("");
  covfile << dataPath() << "rawdata/" << fSetName << "/" << fSetName << ".cov";
  file2.open(covfile.str().c_str(), ios::in);

  if (file2.fail()) {
    cerr << "Error opening data file " << covfile.str() << endl;
    exit(-1);
  }

  string line;

  // Skip comments at top of file
  for (int i = 0; i < 11; i++)
  {
    getline(file1, line);
  }

  double temp, Mll;
  // Filter data file
  for (int i = 0; i < fNData; i++)
  {
    getline(file1, line);
    istringstream lstream(line);

    lstream >> Mll >> temp >> temp >> fData[i] >> temp;
    fKin1[i] = Mll;                               // Mll
    fKin2[i] = pow(Mll,2);                        // Mll^2
    fKin3[i] = 13000;                             // sqrt(s) (GeV)
    fStat[i] = 0;                                 // Statistical uncertainty is included in covariance matrix so set this to 0
  }

  // Skip comments at top of file
  for (int i = 0; i < 5; i++)
  {
    getline(file2, line);
  }

  // Filter covariance matrix file
  double** covmat = new double*[fNData];
  for (int i = 0; i < fNData; i++)
  {
    covmat[i] = new double[fNData];
    getline(file2, line);
    istringstream lstream(line);
    for (int j = 0; j < fNData; j++)
    {
      lstream >> covmat[i][j];
    }
  }

  // Generate artificial systematics
  double** syscor = new double*[fNData];
  for (int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if (!genArtSys(fNData,covmat,syscor))
   {
     cerr << " in " << fSetName << endl;
     exit(-1);
   }

  // Assign artificial systematics
  for (int i = 0; i < fNData; i++)
  {
    for (int j = 0; j < fNSys - 1; j++)
    {
      fSys[i][j].add = syscor[i][j];
      fSys[i][j].mult = fSys[i][j].add*100/fData[i];
      fSys[i][j].type = ADD;
      fSys[i][j].name = "CORR";
    }
  }

  // Assign luminosity uncertainty
  for (int i = 0; i < fNData; i++)
  {
    fSys[i][fNSys-1].mult = 2.3;
    fSys[i][fNSys-1].add = fSys[i][fNSys-1].mult*fData[i]/100;
    fSys[i][fNSys-1].type = MULT;
    fSys[i][fNSys-1].name = "CMSLUMI13";
  }

  file1.close();
  file2.close();

  // Clean-up
  for (int i=0; i<fNData; i++)
    delete[] syscor[i];
  
  delete[] syscor;
  
  for (int i=0; i<fNData; i++)
    delete[] covmat[i];
  
  delete[] covmat;
  
}

/*
Same as CMS_HMDY_13TEV but in the muon channel
*/
void CMS_HMDY_DM_13TEVFilter::ReadData()
{
  // Opening files
  fstream file1, file2;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/" << fSetName << "/" << fSetName << ".data";
  file1.open(datafile.str().c_str(), ios::in);

  if (file1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream covfile("");
  covfile << dataPath() << "rawdata/" << fSetName << "/" << fSetName << ".cov";
  file2.open(covfile.str().c_str(), ios::in);

  if (file2.fail()) {
    cerr << "Error opening data file " << covfile.str() << endl;
    exit(-1);
  }

  string line;

  // Skip comments at top of file
  for (int i = 0; i < 11; i++)
  {
    getline(file1, line);
  }

  double temp, Mll;
  // Filter data file
  for (int i = 0; i < fNData; i++)
  {
    getline(file1, line);
    istringstream lstream(line);

    lstream >> Mll >> temp >> temp >> fData[i] >> temp;
    fKin1[i] = Mll;                               // Mll
    fKin2[i] = pow(Mll,2);                        // Mll^2
    fKin3[i] = 13000;                             // sqrt(s) (GeV)
    fStat[i] = 0;                                 // Statistical uncertainty is included in covariance matrix so set this to 0
  }

  // Skip comments at top of file
  for (int i = 0; i < 5; i++)
  {
    getline(file2, line);
  }

  // Filter covariance matrix file
  double** covmat = new double*[fNData];
  for (int i = 0; i < fNData; i++)
  {
    covmat[i] = new double[fNData];
    getline(file2, line);
    istringstream lstream(line);
    for (int j = 0; j < fNData; j++)
    {
      lstream >> covmat[i][j];
    }
  }

  // Generate artificial systematics
  double** syscor = new double*[fNData];
  for (int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if (!genArtSys(fNData,covmat,syscor))
   {
     cerr << " in " << fSetName << endl;
     exit(-1);
   }

  // Assign artificial systematics
  for (int i = 0; i < fNData; i++)
  {
    for (int j = 0; j < fNSys - 1; j++)
    {
      fSys[i][j].add = syscor[i][j];
      fSys[i][j].mult = fSys[i][j].add*100/fData[i];
      fSys[i][j].type = ADD;
      fSys[i][j].name = "CORR";
    }
  }

  // Assign luminosity uncertainty
  for (int i = 0; i < fNData; i++)
  {
    fSys[i][fNSys-1].mult = 2.3;
    fSys[i][fNSys-1].add = fSys[i][fNSys-1].mult*fData[i]/100;
    fSys[i][fNSys-1].type = MULT;
    fSys[i][fNSys-1].name = "CMSLUMI13";
  }

  file1.close();
  file2.close();

  // Clean-up
  for (int i=0; i<fNData; i++)
    delete[] syscor[i];
  
  delete[] syscor;
  
  for (int i=0; i<fNData; i++)
    delete[] covmat[i];
  
  delete[] covmat;
  
}
