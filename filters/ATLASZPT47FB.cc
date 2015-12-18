/*
 *
 * ATLASZPT47FB - ATLAS Z pt distribution, 7 TeV, 4.7 fb^-1
 * 
 * Reference:
 * HepData:
 *
 */

#include "ATLAS.h"

void ATLASZPT47FBFilter::ReadData()
{
  // Opening files
  fstream f1, f2;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/ATLASZPT47FB.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream datafile2("");

  datafile2 << dataPath() << "rawdata/"
  << fSetName << "/ATLASZPT47FB.cov";
  f2.open(datafile2.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double corr[fNData],toter[fNData],binmin, binmax;

  // Reading data
  string tmp;
  getline(f1,tmp);

  // Filtering data
  for (int i = 0; i < fNData; i++)
  {
    // Data are normalized to the total cross section, Units are GeV-1
    f1 >> binmin >> binmax >> fData[i] >> fStat[i] >> fSys[i][0].add >> corr[i];
    fKin1[i] = binmin + 0.5*(binmax-binmin);
    fKin2[i] = pow(MZ,2.0);
    fKin3[i] = 7E3;
    fStat[i] = fStat[i] * fData[i] * 1e-2; //convert to absolute error
    fSys[i][0].add = fSys[i][0].add * fData[i] * 1e-2; // uncorrelated systematics - convert to absolute error
    corr[i]   = corr[i]   * fData[i] * 1e-2; //convert to absolute error
    toter[i]  = pow(fStat[i]*fStat[i]+fSys[i][0].add*fSys[i][0].add+corr[i]*corr[i],0.5);
  }

  // Reading covmat
  string line;
  double** covmat = new double*[fNData];
  double corrmat[fNData][fNData];
  for(int i = 0; i < fNData; i++)
    {
      covmat[i] = new double[fNData];
      getline(f2,line);
      istringstream lstream(line);
    for(int j = 0; j < fNData; j++)
      {
	lstream >> corrmat[i][j];
	//	covmat[i][j] = corrmat[i][j] * corr[i] * corr[j];   // convert from corr. to cov.
	covmat[i][j] = corrmat[i][j] * toter[i] * toter[j];   // convert from corr. to cov.
      }
    }


  // Generating artificial systematics
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
      for (int l = 1; l < fNSys; l++)
	{
	  fSys[i][l].add  = syscor[i][l-1];
	  fSys[i][l].mult = fSys[i][l].add*100/fData[i];
	  fSys[i][l].type = ADD ;
	  fSys[i][l].name = "CORR";
	}
      fSys[i][0].mult = fSys[i][0].add*100/fData[i];
      fSys[i][0].type = ADD;
      fSys[i][0].name = "UNCORR";
    }
  f1.close();
  f2.close();
}
