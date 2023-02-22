#include "CMSWMU8TEV.h"

/**
 *  1603.01803
 */
void CMSWMU8TEVFilter::ReadData()
{
  // Opening files
  fstream f1, f2;

  stringstream datafileWp("");
  datafileWp << dataPath() << "rawdata/"<< fSetName << "/CMSWMU8TEV-WP.data";
  f1.open(datafileWp.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafileWp.str() << endl;
    exit(-1);
  }

  stringstream datafileWm("");
  datafileWm << dataPath() << "rawdata/"<< fSetName << "/CMSWMU8TEV-WM.data";
  f2.open(datafileWm.str().c_str(), ios::in);

  if (f2.fail()) {
    cerr << "Error opening data file " << datafileWp.str() << endl;
    exit(-1);
  }

  // Total systematic uncertainty
  double systot[fNData];

  // Filter data file
  for (int idat = 0; idat < fNData/2; idat++)
  {
    string line; getline(f1,line);
    istringstream lstream(line);

    double fdum;
    lstream >> fKin1[idat] >> fdum >> fdum >> fData[idat] >> fStat[idat] >> fdum >> systot[idat];
    fKin2[idat] = pow(MW,2);
    fKin3[idat] = 8e3;         // LHC ( TeV)
    fData[idat] *= 1000.;      // Convert from pb (datafiles) to fb (APPLgrid)
    fStat[idat] *= 1000.;
    systot[idat] *= 1000.;
  }

  for (int idat = fNData/2; idat < fNData; idat++)
  {
    string line; getline(f2,line);
    istringstream lstream(line);

    double fdum;
    lstream >> fKin1[idat] >> fdum >> fdum >> fData[idat] >> fStat[idat] >> fdum >> systot[idat];
    fKin2[idat] = pow(MW,2);
    fKin3[idat] = 8e3;         // LHC ( TeV)
    fData[idat] *= 1000.;      // Convert from pb (datafiles) to fb (APPLgrid)
    fStat[idat] *= 1000.;
    systot[idat] *= 1000.;
  }

  f1.close();
  f2.close();

  // Generate artificial systematics for systematic covmat
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++) syscor[i] = new double[fNData];  
  GenArtSys("sys", systot, syscor);

  // Generate artificial systematics for statistical covmat
  double** statcor = new double*[fNData];
  for(int i = 0; i < fNData; i++) statcor[i] = new double[fNData];  
  GenArtSys("sys", fStat, statcor);

  for (int i = 0; i < fNData; i++)
  {
    for (int l = 0; l < fNData; l++)
    {
      // Systematic covariance matrix
      fSys[i][l].add = syscor[i][l];
      fSys[i][l].mult = fSys[i][l].add*100/fData[i];
      fSys[i][l].type = MULT;
      fSys[i][l].name = "CORR";
      // Statistical covariance matrix
      fSys[i][l+fNData].add = statcor[i][l];
      fSys[i][l+fNData].mult = fSys[i][l+fNData].add*100/fData[i];
      fSys[i][l+fNData].type = ADD;
      fSys[i][l+fNData].name = "CORR";
    }

    // Clear stat
    fStat[i] = 0.0;

    // Luminosity Uncertainty
    // CMS Luminosity Uncertainty, 2012 data set: 2.6%
    // http://cds.cern.ch/record/1598864?ln=en
    fSys[i][fNSys-1].mult = 2.6;
    fSys[i][fNSys-1].add  = fSys[i][fNSys-1].mult*fData[i]*1e-2;
    fSys[i][fNSys-1].type = MULT;
    fSys[i][fNSys-1].name = "CMSLUMI12";
  }

  // Cleanup
  for(int i = 0; i < fNData; i++)
  {
    delete[] syscor[i];
    delete[] statcor[i];
  }

  delete[] syscor;
  delete[] statcor;
}

// Process artificial systematics
void CMSWMU8TEVFilter::GenArtSys(std::string const& type, double* std, double** artsys)
{
  // Read Matrix
  stringstream matfile; fstream ms;
  matfile << dataPath() << "rawdata/"<< fSetName << "/CMSWMU8TEV."<<type<<"corr";
  ms.open(matfile.str().c_str(), ios::in);

  if (ms.fail()) {
    cerr << "Error opening data file " << matfile.str() << endl;
    exit(-1);
  }

  // Filter Correlation Matrix
  double** covmat = new double*[fNData];
  for(int i = 0; i < fNData; i++)
  {
    covmat[i] = new double[fNData];
    string line; getline(ms,line);
    istringstream lstream(line);
    for(int j = 0; j < fNData; j++)
    {
      double entry = 0.0; lstream >> entry;             // Read correlation
      covmat[i][j] = (entry/100.)*std[i]*std[j];  // Compute covariance mat. from correlation mat.
    }
  }

  // Generate artificial systematics
  if(!genArtSys(fNData,covmat,artsys))
   {
     cerr << " in " << fSetName << endl;
     exit(-1);
   }

  // Free covmat
  for(int i = 0; i < fNData; i++)
    delete[] covmat[i];
  delete[] covmat;
}

