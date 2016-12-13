/*
 * ATLASZPT7TEV - ATLAS Z pt distribution, 7 TeV, 4.7 fb^-1
 * This is the implementation of a three separate bins in rapidity
 * Reference: https://arxiv.org/abs/1406.3660
 * HepData: http://hepdata.cedar.ac.uk/view/ins1300647
 *
 */

#include "ATLASZPT7TEV.h"

void ATLASZPT7TEVFilter::ReadData()
{
  // Opening data file
  fstream f1, c1, c2, c3;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/ATLASZPT47FB.data";
  f1.open(datafile.str().c_str(), ios::in);
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Opening correlation files
  stringstream corrbin1("");
  corrbin1 << dataPath() << "rawdata/"
  << fSetName << "/CombinedCorrelations_yb1_qedBorn.txt";
  c1.open(corrbin1.str().c_str(), ios::in);
  if (c1.fail()) {
    cerr << "Error opening data file " << corrbin1.str() << endl;
    exit(-1);
  }
  stringstream corrbin2("");
  corrbin2 << dataPath() << "rawdata/"
  << fSetName << "/CombinedCorrelations_yb2_qedBorn.txt";
  c2.open(corrbin2.str().c_str(), ios::in);
  if (c2.fail()) {
    cerr << "Error opening data file " << corrbin2.str() << endl;
    exit(-1);
  }
  stringstream corrbin3("");
  corrbin3 << dataPath() << "rawdata/"
  << fSetName << "/CombinedCorrelations_yb3_qedBorn.txt";
  c3.open(corrbin3.str().c_str(), ios::in);
  if (c3.fail()) {
    cerr << "Error opening data file " << corrbin3.str() << endl;
    exit(-1);
  }

  // Reading data
  string tmp;
  getline(f1,tmp);

  // Filtering data
  int DataBin = fNData/3;
  double corr[fNData],toter[fNData],binmin, binmax,dummy;
  double DataDressed[fNData],StatDressed[fNData],SysUncDressed[fNData],SysCorDressed[fNData];

  for (int i = 0; i < DataBin; i++)
  {
    // Data are normalized to the total cross section, units are GeV-1, read only data
    // Errors are already given as absolute numbers in the HEPDATA table format
    f1 >> fKin2[i] >> binmin >> binmax >> fData[i] >> fStat[i] >> dummy >> fSys[i][0].add >> dummy >> corr[i] >> dummy >> DataDressed[i] >> StatDressed[i] >> dummy >> SysUncDressed[i] >> dummy >> SysCorDressed[i] >> dummy >> fData[i+DataBin] >> fStat[i+DataBin] >> dummy >> fSys[i+DataBin][0].add >> dummy >> corr[i+DataBin] >> dummy >> DataDressed[i+DataBin] >> StatDressed[i+DataBin] >> dummy >> SysUncDressed[i+DataBin] >> dummy >> SysCorDressed[i+DataBin] >> dummy >> fData[i+2 * DataBin] >> fStat[i+ 2* DataBin] >> dummy >> fSys[i+ 2*DataBin][0].add >> dummy >> corr[i+2*DataBin] >> dummy >> DataDressed[i+2*DataBin] >> StatDressed[i+2*DataBin] >> dummy >> SysUncDressed[i+2*DataBin] >> dummy >> SysCorDressed[i+2*DataBin] >> dummy;
    fKin2[i] *= fKin2[i];
    fKin2[i+DataBin] = fKin2[i];
    fKin2[i+2*DataBin] = fKin2[i];
  }

  for (int i = 0; i < fNData; i++)
    {
      fKin1[i] = i < 26 ? 0.5 : (i < 52 ? 1.5 : 2.25);
      fKin3[i] = 7E3;
      toter[i]  = pow(fStat[i]*fStat[i]+fSys[i][0].add*fSys[i][0].add+corr[i]*corr[i],0.5);
    }

  // Initialize covariance matrix
  double** covmat = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    {
      covmat[i] = new double[fNData];
      for(int j = 0; j < fNData; j++)
	{
	  covmat[i][j] = 0.;
	}
    }

  // Reading covariance matrix of the first bin (top left 26 x 26 elements)
  string line1,line2,line3;
  double corrmat1[DataBin][DataBin];
  double corrmat2[DataBin][DataBin];
  double corrmat3[DataBin][DataBin];
  for(int i = 0; i < DataBin; i++)
    {
      covmat[i] = new double[fNData];
      getline(c1,line1);
      getline(c2,line2);
      getline(c3,line3);
      istringstream lstream1(line1);
      istringstream lstream2(line2);
      istringstream lstream3(line3);
      for(int j = 0; j < DataBin; j++)
	{
	  lstream1 >> corrmat1[i][j];
	  lstream2 >> corrmat2[i][j];
	  lstream3 >> corrmat3[i][j];
	  covmat[i][j] = corrmat1[i][j] * toter[i] * toter[j];   // convert from corr. to cov. multiplying by total error
	  covmat[i+DataBin][j+DataBin] = corrmat2[i][j] * toter[i+DataBin] * toter[j+DataBin];
	  covmat[i+2*DataBin][j+2*DataBin] = corrmat3[i][j] * toter[i+2*DataBin] * toter[j+2*DataBin];
	}
    }

  /* CHECK
     for(int i = 0; i < DataBin; i++)
     {
     for(int j = 0; j < DataBin; j++)
     {
     cout << i << "  " << j << "  " << covmat[i][j] << " " << toter[i]*toter[j] << endl;
     }
     }
     exit(-1);
  */

  // Generating artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
   {
     cerr << " in " << fSetName << endl;
     exit(-1);
   }

  // The uncorrelated systematics is assigned to fSys[i][0]
  for (int i = 0; i < fNData; i++)
    {
      for (int l = 0; l < fNSys; l++)
	{
	  fSys[i][l].add  = syscor[i][l-1];
	  fSys[i][l].mult = fSys[i][l].add*100/fData[i];
	  fSys[i][l].type = ADD;
	  fSys[i][l].name = (l == 0 ? "UNCORR" : "CORR");
	}
      /*
      fSys[i][0].mult = fSys[i][0].add*100/fData[i];
      fSys[i][0].type = ADD;
      fSys[i][0].name = "UNCORR";
      */
      fStat[i] = 0.0;
    }

  f1.close();
  c1.close();
  c2.close();
  c3.close();
}

/*
 *
 * ATLASZPT47FB - ATLAS Z pt distribution, 7 TeV, 4.7 fb^-1
 * This is the implementation of a single bin in rapidity
 * Reference:
 * HepData:
 *
 */

/*
#include "ATLASZPT47FB.h"

void ATLASZPT47FBFilter::ReadData()
{
  // Opening files
  fstream f1, f2;

  // From HEPDATA, Born colums of table 2 in http://hepdata.cedar.ac.uk/view/ins1300647
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/ATLASZPT47FB_inclusive.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream datafile2("");

  datafile2 << dataPath() << "rawdata/"
  << fSetName << "/ATLASZPT47FB_inclusive.cov";
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
*/
