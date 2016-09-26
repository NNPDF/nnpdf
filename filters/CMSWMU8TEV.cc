#include "CMSWMU8TEV.h"

/**
 *
 */
void CMSWMU8TEVFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3, f4;

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

  stringstream syscorfile("");
  syscorfile << dataPath() << "rawdata/"<< fSetName << "/CMSWMU8TEV.syscorr";
  f3.open(syscorfile.str().c_str(), ios::in);

  if (f3.fail()) {
    cerr << "Error opening data file " << syscorfile.str() << endl;
    exit(-1);
  }

  stringstream statcorfile("");
  statcorfile << dataPath() << "rawdata/"<< fSetName << "/CMSWMU8TEV.statcorr";
  f4.open(statcorfile.str().c_str(), ios::in);

  if (f4.fail()) {
    cerr << "Error opening data file " << statcorfile.str() << endl;
    exit(-1);
  }

  string line;
  float fdum;
  float systot[fNData];

  // Filter data file
  for (int idat = 0; idat < fNData/2; idat++)
  {
    getline(f1,line);
    istringstream lstream(line);

    lstream >> fKin1[idat] >> fdum >> fdum >> fData[idat] >> fStat[idat] >> fdum >> systot[idat];
    fKin2[idat] = pow(MW,2);
    fKin3[idat] = 8e3;         // LHC ( TeV)
    fData[idat] *= 1000.;      // Convert from pb (datafiles) to fb (APPLgrid)
    fStat[idat] *= 1000.;
    systot[idat] *= 1000.;
  }

  for (int idat = fNData/2; idat < fNData; idat++)
  {
    getline(f2,line);
    istringstream lstream(line);

    lstream >> fKin1[idat] >> fdum >> fdum >> fData[idat] >> fStat[idat] >> fdum >> systot[idat];
    fKin2[idat] = pow(MW,2);
    fKin3[idat] = 8e3;         // LHC ( TeV)
    fData[idat] *= 1000.;      // Convert from pb (datafiles) to fb (APPLgrid)
    fStat[idat] *= 1000.;
    systot[idat] *= 1000.;
  }

  // Filter Correlation Matrix

  double** corrmat = new double*[fNData];
  double** covmat = new double*[fNData];
  for(int i = 0; i < fNData; i++)
  {
    corrmat[i] = new double[fNData];
    covmat[i] = new double[fNData];
    getline(f3,line);
    istringstream lstream(line);
    for(int j = 0; j < fNData; j++)
    {
      lstream >> corrmat[i][j];
      corrmat[i][j] /=100.;                              // Correlation matrix entries given percentage
      covmat[i][j] = corrmat[i][j]*systot[i]*systot[j];  // Compute covariance mat. from correlation mat.
    }
  }

  // Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
   {
     cerr << " in " << fSetName << endl;
     exit(-1);
   }

  for (int i = 0; i < fNData; i++)
  {
    for (int l = 0; l < fNSys-1; l++)
    {
      fSys[i][l].add = syscor[i][l];
      fSys[i][l].mult = fSys[i][l].add*100/fData[i];
      fSys[i][l].type = MULT;
      fSys[i][l].name = "CORR";
    }
    // Luminosity Uncertainty
    // CMS Luminosity Uncertainty, 2012 data set: 2.6%
    // http://cds.cern.ch/record/1598864?ln=en

    fSys[i][fNSys-1].mult = 2.6;
    fSys[i][fNSys-1].add  = fSys[i][fNSys-1].mult*fData[i]*1e-2;
    fSys[i][fNSys-1].type = MULT;
    fSys[i][fNSys-1].name = "CMSLUMI12";
  }

  f1.close();
  f2.close();
  f3.close();
  f4.close();
}
