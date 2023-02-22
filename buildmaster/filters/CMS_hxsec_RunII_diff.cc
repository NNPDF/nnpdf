/*
Reference:
   arXiv:     [1812.06504]
   hepdata:   https://www.hepdata.net/record/ins1709330
   published: Phys.Lett. B792 (2019) 369-396 2019 
Description:
   Differential spectra for the combined Higgs boson production and decay in the
   H -> gamma gamma, H -> Z Z and H -> b bbar channels, measured as a function 
   of the transverse momentum or rapidity of the Higgs. An additional 
   measurement related to the ggH process is implemented separately. Cross
   sections are measured by the CMS detector at the LHC, with a c.m. energy of
   13 TeV and an integrated luminosity of 35.9 fb-1.
   
   The implementation is based on Tabs. A1, A2 and A4 (central value and 
   uncertainties) and Tabs. B1 and B3 (correlations). Three data sets are 
   defined separately:
    - CMS_hxsec_RunII_pTH: combined pT spectrum
    - CMS_hxsec_RunII_pTH_ggH: ggH pT spectrum
    - CMS_hxsec_RunII_yH: combined yH spectrum
   The information on the signal strengths is retrieved from the hepdata entry
*/

#include "CMS_hxsec_RunII_diff.h"

//CMS_hxsec_RunII_pTH: combined pT spectrum
void CMS_hxsec_RunII_diff_pTHFilter::ReadData()
{
  fstream f1;
  fstream f2;

  //Central values and total uncertainty
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/CMS_hxsec_RunII_diff/data_pTH.txt";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Correlations between data points
  stringstream datafile_corr("");
  datafile_corr << dataPath()
		<< "rawdata/CMS_hxsec_RunII_diff/corr_pTH.txt";
  f2.open(datafile_corr.str().c_str(), ios::in);

  if (f2.fail())
    {
      cerr << "Error opening data file " << datafile_corr.str() << endl;
      exit(-1);
    }

  //Read central values and total uncertainties
  string line;

  double* Sys = new double[fNData];
  double* SM_cov = new double[fNData];
  double** corrmat = new double*[fNData];
  double** syscor  = new double*[fNData];
  
  for(int i=0; i<fNData; i++)
    {
      double ddum;
      double SS_cv, SS_er;
      getline(f1,line);
      istringstream lstream(line);
      fKin1[i] = 0.;
      fKin2[i] = 0.;
      fKin3[i] = 0.;
      lstream >> ddum >> ddum >> fData[i] >> Sys[i]
	      >> SS_cv >> SS_er;
      
      //cout << fData[i]/SS_cv << endl;
      SM_cov[i] = fData[i]/SS_cv * SS_er/SS_cv;

      fStat[i] = 0.0;
	     
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

  /*
  //Generate SM covariance matrix
  for(int i=0; i< fNData; i++)
    {
      for (int j=0; j< fNData; j++)
	{
	  cout << SM_cov[i]*SM_cov[j] << "   ";
	}
      cout << endl;
    }
  */

  for(int i = 0; i < fNData; i++) 
    {
      delete[] syscor[i];
      delete[] corrmat[i];
    }
  
  delete[] SM_cov;
  delete[] Sys;
  delete[] syscor;
  delete[] corrmat;
  
  f1.close();
  f2.close();
}

//CMS_hxsec_RunII_pTH_ggH: ggH pT spectrum
void CMS_hxsec_RunII_diff_pTH_ggHFilter::ReadData()
{
  fstream f1;
  fstream f2;

  //Central values and total uncertainty
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/CMS_hxsec_RunII_diff/data_pTH_ggH.txt";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Correlations between data points
  stringstream datafile_corr("");
  datafile_corr << dataPath()
		<< "rawdata/CMS_hxsec_RunII_diff/corr_pTH_ggH.txt";
  f2.open(datafile_corr.str().c_str(), ios::in);

  if (f2.fail())
    {
      cerr << "Error opening data file " << datafile_corr.str() << endl;
      exit(-1);
    }

  //Read central values and total uncertainties
  string line;

  double* Sys = new double[fNData];
  double* SM_cov = new double[fNData];
  double** corrmat = new double*[fNData];
  double** syscor  = new double*[fNData];
  
  for(int i=0; i<fNData; i++)
    {
      double ddum;
      double SS_cv, SS_er;
      getline(f1,line);
      istringstream lstream(line);
      fKin1[i] = 0.;
      fKin2[i] = 0.;
      fKin3[i] = 0.;
      lstream >> ddum >> ddum >> fData[i] >> Sys[i]
	      >> SS_cv >> SS_er;

      //cout << fData[i]/SS_cv << endl;
      SM_cov[i] = fData[i]/SS_cv * SS_er/SS_cv;

      fStat[i] = 0.0;
	     
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

  /*
  //Generate SM covariance matrix
  for(int i=0; i< fNData; i++)
    {
      for (int j=0; j< fNData; j++)
	{
	  cout << SM_cov[i]*SM_cov[j] << "   ";
	}
      cout << endl;
    }
  */
  
  for(int i = 0; i < fNData; i++) 
    {
      delete[] syscor[i];
      delete[] corrmat[i];
    }
  
  delete[] SM_cov;
  delete[] Sys;
  delete[] syscor;
  delete[] corrmat;
  
  f1.close();
  f2.close();

}

//CMS_hxsec_RunII_yH: combined yH spectrum
void CMS_hxsec_RunII_diff_yHFilter::ReadData()
{
  fstream f1;
  fstream f2;

  //Central values and total uncertainty
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/CMS_hxsec_RunII_diff/data_yH.txt";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Correlations between data points
  stringstream datafile_corr("");
  datafile_corr << dataPath()
		<< "rawdata/CMS_hxsec_RunII_diff/corr_yH.txt";
  f2.open(datafile_corr.str().c_str(), ios::in);

  if (f2.fail())
    {
      cerr << "Error opening data file " << datafile_corr.str() << endl;
      exit(-1);
    }

  //Read central values and total uncertainties
  string line;

  double* Sys = new double[fNData];
  double* SM_cov = new double[fNData];
  double** corrmat = new double*[fNData];
  double** syscor  = new double*[fNData];
  
  for(int i=0; i<fNData; i++)
    {
      double ddum;
      double SS_cv, SS_er;
      getline(f1,line);
      istringstream lstream(line);
      fKin1[i] = 0.;
      fKin2[i] = 0.;
      fKin3[i] = 0.;
      lstream >> ddum >> ddum >> fData[i] >> Sys[i]
	      >> SS_cv >> SS_er;

      //cout << fData[i]/SS_cv << endl;
      SM_cov[i] = fData[i]/SS_cv * SS_er/SS_cv;

      fStat[i] = 0.0;
	     
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

  /*
  //Generate SM covariance matrix
  for(int i=0; i< fNData; i++)
    {
      for (int j=0; j< fNData; j++)
	{
	  cout << SM_cov[i]*SM_cov[j] << "   ";
	}
      cout << endl;
    }
  */

  for(int i = 0; i < fNData; i++) 
    {
      delete[] syscor[i];
      delete[] corrmat[i];
    }
  delete[] SM_cov;
  delete[] Sys;
  delete[] syscor;
  delete[] corrmat;
  
  f1.close();
  f2.close();
}

