#include "CMSDY2D12.h"

/**
 *
 */
void CMSDY2D12Filter::ReadData()
{
  // Opening files
  fstream f1, f2;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"<< fSetName << "/CMSDY2D12_dsigdy.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream covfile("");
  covfile << dataPath() << "rawdata/"<< fSetName << "/CMSDY2D12_dsigdy.cov";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail()) {
    cerr << "Error opening data file " << covfile.str() << endl;
    exit(-1);
  }


  string line;
  int idum;
  double Mll;

  // Filter data file
  for (int idat = 0; idat < fNData; idat++)
  {
    getline(f1,line);
    istringstream lstream(line);

    lstream >> idum >> Mll >> fKin1[idat] >> fData[idat];
    fKin2[idat] = pow(Mll,2);
    fKin3[idat] = 8e3;         // LHC ( TeV)
    fData[idat] *= 1000.;      // Convert from pb (datafiles) to fb (APPLgrid)
    fStat[idat] = 0; // Statistical errors are correlated

  }

  // Filter Covariance Matrix

  double** covmat = new double*[fNData];
  for(int i = 0; i < fNData; i++)
  {
    covmat[i] = new double[fNData];
    getline(f2,line);
    istringstream lstream(line);
    for(int j = 0; j < fNData; j++)
    {
      lstream >> covmat[i][j];
      covmat[i][j] *= 1e6;              // correct for units
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

  for(int i = 0; i < fNData; i++) 
    {
      delete[] syscor[i];
      delete[] covmat[i];
    }
  delete[] syscor;
  delete[] covmat;
  
  f1.close();
  f2.close();
}
