/*
Reference:
   [1909.02845]
   Combined measurements of Higgs boson production and decay using up to
   80 fb-1 of proton–proton collision data at √s = 13 TeV collected with 
   the ATLAS experiment
   Phys. Rev. D101 012002

There are 16 signal strengths in the following order
ggF gamma gamma
    Z Z
    W W
    tau tau
VBF gamma gamma
    Z Z
    W W
    tau tau
    b bbar
VH  gamma gamma
    Z Z
    b bar
tth gamma gamma
    V V
    tau tau
    b b

The implementation is based on Fig.5 (including uncertainties) and on Fig.6 
(the correlation matrix for the total uncertainty).
*/

#include "ATLAS_hxsec_RunII.h"

void ATLAS_hxsec_RunIIFilter::ReadData()
{
  fstream f1;
  fstream f2;

  //Central values and total uncertainty
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/ATLAS_hxsec_RunII/data.txt";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Correlations between the 16 data points
  stringstream datafile_corr("");
  datafile_corr << dataPath()
		<< "rawdata/ATLAS_hxsec_RunII/corr.txt";
  f2.open(datafile_corr.str().c_str(), ios::in);

  if (f2.fail())
    {
      cerr << "Error opening data file " << datafile_corr.str() << endl;
      exit(-1);
    }

  //Read central values and total uncertainties
  string line;

  double* Sys = new double[fNData];
  double** corrmat = new double*[fNData];
  double** syscor  = new double*[fNData];
  
  for(int i=0; i<fNData; i++)
    {
   
      getline(f1,line);
      istringstream lstream(line);
      fKin1[i] = 0.;
      fKin2[i] = 0.;
      fKin3[i] = 0.;
      lstream >> fData[i] >> fStat[i] >> Sys[i];
	     
      corrmat[i] = new double[fNData];
      syscor[i]  = new double[fNData];
      getline(f2,line);
      istringstream kstream(line);

      for(int j=0; j<fNData; j++)
	{
	  kstream >> corrmat[i][j];
	}
    }

  //Generate covariance matrix from correlation matrix
  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  corrmat[i][j] = corrmat[i][j]*Sys[i]*Sys[j];
	}
    }

  //Generate artificial systematics from covariance matrix
  if(!genArtSys(fNData,corrmat,syscor))
    {
      throw runtime_error("Couldn't generate artificial systematics for " + fSetName);
    }

  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNSys; j++)
	{
	  fSys[i][j].add  = syscor[i][j];
	  fSys[i][j].mult = fSys[i][j].add*1e2/fData[i];
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = "CORR";
	}
    } 

  for(int i = 0; i < fNData; i++) 
    {
      delete[] syscor[i];
      delete[] corrmat[i];
    }
  delete[] Sys;
  delete[] syscor;
  delete[] corrmat;
  
  f1.close();
  f2.close();

}
