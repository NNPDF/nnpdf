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
  << fSetName << "/ATLASZPT47FB_yoda.dat";
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

  // Data are normalized to the total cross section, units are GeV-1, read only data
  // Errors are already given as absolute numbers in the HEPDATA table format
  double toter[fNData], dummy;
  for (int i = 0; i < fNData; i++)
  {
    f1 >> fKin2[i] >> dummy >> dummy >> fData[i] >> dummy >> toter[i];

    fKin1[i] = i < 26 ? 0.5 : (i < 52 ? 1.5 : 2.25);
    fKin2[i] *= fKin2[i];
    fKin3[i] = 7E3;

    fStat[i] = 0;
  }

  // Initialize covariance matrix
  double** covmat = new double*[fNData];
  for(int i = 0; i < fNData; i++)
  {
    covmat[i] = new double[fNData];
    for(int j = 0; j < fNData; j++)
      covmat[i][j] = 0.;
  }

  // Reading covariance matrix of the first bin (top left 26 x 26 elements)
  string line1,line2,line3;
  const int DataBin = fNData/3;
  double corrmat1[DataBin][DataBin];
  double corrmat2[DataBin][DataBin];
  double corrmat3[DataBin][DataBin];
  for(int i = 0; i < DataBin; i++)
    {
      //covmat[i] = new double[fNData];
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

  // Generating artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];
  if(!genArtSys(fNData,covmat,syscor))
  {
    cerr << " in " << fSetName << endl;
    exit(-1);
  }

  // Assign artificial systematics
  for (int i = 0; i < fNData; i++)
    for (int l = 0; l < fNSys; l++)
  	{
  	  fSys[i][l].add  = syscor[i][l];
  	  fSys[i][l].mult = fSys[i][l].add*100/fData[i];
  	  fSys[i][l].type = MULT;
  	  fSys[i][l].name = "CORR";
  	}
  
  for(int i = 0; i < fNData; i++) 
    {
      delete[] syscor[i];
      delete[] covmat[i];
    }
  delete[] syscor;
  delete[] covmat;

  f1.close();
  c1.close();
  c2.close();
  c3.close();
}
