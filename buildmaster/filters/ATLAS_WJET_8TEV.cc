/*
Name_exp : ATLAS_WJET_8TEV
Reference: Measurement of differential cross sections and 𝑊+/𝑊−
           cross-section ratios for 𝑊 boson production in association 
           with jets at 𝑠√=8 TeV with the ATLAS detector
ArXiv    : arxiv:1711.03296
Published: JHEP 1805 (2018) 077
Hepdata  : https://hepdata.net/record/ins1635273

measurement of the W boson production cross section and the W+/W- cross-section 
ratio, both in association with jets, in proton–proton collisions at s=sqrt{8} 
TeV with the ATLAS experiment at the Large Hadron Collider. The measurement is 
performed in final states containing one electron and missing transverse 
momentum using data corresponding to an integrated luminosity of 20.2 fb-1. 
Differential cross sections for events with at least one or two jets are 
presented for a range of observables, including jet transverse momenta 
and the transverse momentum of the W boson. The differential cross sections of 
positively and negatively charged W bosons are measured separately.

Four distributions are considered in the following:
1) W+ + jet, distribution differential in the transverse momentum 
   of the W boson;
2) W- + jet, distribution differential in the transverse momentum 
   of the W boson;
3) W+ + jet, distribution differential in the transverse momentum 
   of the leading jet;
4) W- + jet, distribution differential in the transverse momentum 
   of the leading jet;

The information on experimental uncertainties is retrieved from the hepdata 
entry. Each source of systematic uncertainty (except unfolding uncertainties) 
is assumed to be bin-by-bin correlated within each distribution and between W+ 
and W- production within the same distribution. A statistical correlation 
matrix is implemented to account for correlations of the statistical 
uncertainty within each distribution. Non perturbative corrections (optional) 
are implemented as an extracorrelated source of systematic uncertainty that 
deweights the data.
*/

#include "ATLAS_WJET_8TEV.h"

//1)W+ distribution differential in the transverse momentum of the W boson

void ATLAS_WP_JET_8TEV_PTFilter::ReadData()
{
  fstream f1;
  fstream f2;
  fstream f3;

  //Full breakdown of statistical uncertainties
  stringstream datafile_sys("");
  datafile_sys << dataPath()
	   << "rawdata/ATLAS_WJET_8TEV/HEPData-ins1635273-v1-Table_13.csv";
  f1.open(datafile_sys.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile_sys.str() << endl;
      exit(-1);
    }

  //Correlation matrix of statistical ucnertainties
  stringstream datafile_stat("");
  datafile_stat << dataPath()
		<< "rawdata/ATLAS_WJET_8TEV/HEPData-ins1635273-v1-Table_14.csv";

  f2.open(datafile_stat.str().c_str(), ios::in);

  if (f2.fail())
    {
      cerr << "Error opening data file " << datafile_stat.str() << endl;
      exit(-1);
    }

  //Non-perturbative corrections
  stringstream datafile_np("");
  datafile_np << dataPath()
		<< "rawdata/ATLAS_WJET_8TEV/HEPData-ins1635273-v1-Table_67.csv";

  f3.open(datafile_np.str().c_str(), ios::in);

  if (f3.fail())
    {
      cerr << "Error opening data file " << datafile_np.str() << endl;
      exit(-1);
    }

  //Read central value
  string line;
  for(int i=0; i<16; i++)
    {
      getline(f1,line);
      getline(f2,line);
      getline(f3,line);
    }

  double** corrmat = new double*[fNData];
  double** syscor  = new double*[fNData];
  double** sysR    = new double*[fNData];
  double** sysL    = new double*[fNData];
  double npcorr[fNData];
  const int nrealsys = 54;
  double ddum;
  char comma;

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      lstream >> fKin1[i] >> comma
	      >> ddum     >> comma
	      >> ddum     >> comma
	      >> fData[i] >> comma
	      >> fStat[i] >> comma
	      >> ddum     >> comma;

      fKin2[i] = 0;
      fKin3[i] = 8000; //GeV

      sysR[i]    = new double[nrealsys];
      sysL[i]    = new double[nrealsys];

      for(int k=0; k<nrealsys; k++)
	{ 
	  lstream >> sysR[i][k] >> comma >> sysL[i][k] >> comma;
	}

      corrmat[i] = new double[fNData];
      syscor[i]  = new double[fNData];

      for(int j=0; j<fNData; j++)
	{
	  getline(f2,line);
	  istringstream istream(line);
	  istream >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> corrmat[i][j];
	}

      getline(f3,line);
      istringstream kstream(line);
      kstream >> ddum >> comma
	      >> ddum >> comma
	      >> ddum >> comma
	      >> npcorr[i];

    }

  //Generate artificial systematics
  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  corrmat[i][j] *= fStat[i]*fStat[j];
	}
    }
  
  if(!genArtSys(fNData,corrmat,syscor))
    {
      throw runtime_error("Couldn't generate artificial systematics for " + fSetName);
    }

  for(int i=0; i<fNData; i++)
    {
      fStat[i]=0.;
      for(int j=0; j<fNData; j++)
	{
	  fSys[i][j].add  = syscor[i][j];
	  fSys[i][j].mult = fSys[i][j].add*1e2/fData[i];
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = "CORR";
	}
    }

  //Real systematics (correlated across bins and W+/W- distributions)
  for(int i=0; i<fNData; i++)
    {
      
      for(int k=0; k<50; k++)
	{
	  sysR[i][k] /= sqrt(2.);
	  sysL[i][k] /= sqrt(2.);

	  double tmp1, tmp2;
	  tmp1=sysR[i][k];
	  tmp2=sysL[i][k];
	  
	  //Case 1: sysL and sysR are both negative
	  if(tmp1<0 && tmp2<0)
	    {
	      if(tmp2<tmp1 || tmp2==tmp1)
		{
		  sysR[i][k] = 0.0;
		  sysL[i][k] = tmp2;
		}
	      if(tmp2>tmp1)
		{
		  sysR[i][k] = 0.0;
		  sysL[i][k] = tmp1;
		}
	    }

	  //Case 2: sysL and sysR are both positive
	  if(tmp1>0.0 && tmp2>0.0)
	    {
	      if(tmp1>tmp2 || tmp1==tmp2)
		{
		  sysR[i][k] = tmp1;
		  sysL[i][k] = 0.0;
		}
	      if(tmp1<tmp2)
		{
		  sysR[i][k] = tmp2;
		  sysL[i][k] = 0.0;
		}
	    }

	  //Case3: sys1 is negative and sys2 is positive
	  if(tmp1<0.0 && tmp2>0.0)
	    {
	      sysR[i][k] = tmp2;
	      sysL[i][k] = tmp1;
	    }
	  
	  fSys[i][2*k+fNData].add = sysR[i][k];
	  fSys[i][2*k+fNData].mult = fSys[i][2*k+fNData].add*1e2/fData[i];
	  fSys[i][2*k+fNData].type = MULT;
	  ostringstream sysnameR;
	  sysnameR << "ATLASWJ" << 2*k;
	  fSys[i][2*k+fNData].name = sysnameR.str();
	  
	  fSys[i][2*k+1+fNData].add = sysL[i][k];
	  fSys[i][2*k+1+fNData].mult = fSys[i][2*k+1+fNData].add*1e2/fData[i];
	  fSys[i][2*k+1+fNData].type = MULT;
	  ostringstream sysnameL;
	  sysnameL << "ATLASWJ" << 2*k+1;
	  fSys[i][2*k+1+fNData].name = sysnameL.str();

	}

      //Luminosity uncertainty
      fSys[i][2*50+fNData].add  = sysR[i][50];
      fSys[i][2*50+fNData].mult = fSys[i][2*50+fNData].add*1e2/fData[i];
      fSys[i][2*50+fNData].type = MULT;
      fSys[i][2*50+fNData].name = "ATLASLUMI12";

      //Uncorrelated unfolding uncertainties
      for(int k=51; k<nrealsys; k++)
	{
	  sysR[i][k] /= sqrt(2.);
	  sysL[i][k] /= sqrt(2.);
	  
	  double tmp1, tmp2;
	  tmp1=sysR[i][k];
	  tmp2=sysL[i][k];
	  
	  //Case 1: sysL and sysR are both negative
	  if(tmp1<0 && tmp2<0)
	    {
	      if(tmp2<tmp1 || tmp2==tmp1)
		{
		  sysR[i][k] = 0.0;
		  sysL[i][k] = tmp2;
		}
	      if(tmp2>tmp1)
		{
		  sysR[i][k] = 0.0;
		  sysL[i][k] = tmp1;
		}
	    }

	  //Case 2: sysL and sysR are both positive
	  if(tmp1>0.0 && tmp2>0.0)
	    {
	      if(tmp1>tmp2 || tmp1==tmp2)
		{
		  sysR[i][k] = tmp1;
		  sysL[i][k] = 0.0;
		}
	      if(tmp1<tmp2)
		{
		  sysR[i][k] = tmp2;
		  sysL[i][k] = 0.0;
		}
	    }
	  
	  //Case3: sys1 is negative and sys2 is positive
	  if(tmp1<0.0 && tmp2>0.0)
	    {
	      sysR[i][k] = tmp2;
	      sysL[i][k] = tmp1;
	    }
	  	  
	  fSys[i][2*k+fNData-1].add  = sysR[i][k];
	  fSys[i][2*k+fNData-1].mult = fSys[i][2*k+fNData-1].add*1e2/fData[i];
	  fSys[i][2*k+fNData-1].type = MULT;
	  fSys[i][2*k+fNData-1].name = "UNCORR";
	  
	  fSys[i][2*k+fNData].add  = sysL[i][k];
	  fSys[i][2*k+fNData].mult = fSys[i][2*k+fNData].add*1e2/fData[i];
	  fSys[i][2*k+fNData].type = MULT;
	  fSys[i][2*k+fNData].name = "UNCORR";
	}
      
      //Non-perturbative corrections
      fSys[i][123].add  = fData[i]*(1. - npcorr[i]);
      fSys[i][123].mult = fSys[i][123].add*1e2/fData[i];
      fSys[i][123].type = MULT;
      fSys[i][123].name = "SKIP";

    }

 // Clean-up
  for (int i=0; i<fNData; i++)
    delete[] syscor[i];
  
  delete[] syscor;

  for(int i=0; i<fNData; i++)
    delete[] corrmat[i];

  delete[] corrmat;

  for(int i=0; i<fNData; i++)
    delete[] sysR[i];

  delete[] sysR;

  for(int i=0; i<fNData; i++)
    delete[] sysL[i];

  delete[] sysL;
  
  f1.close();
  f2.close();

} 

//2)W- distribution differential in the transverse momentum of the W boson

void ATLAS_WM_JET_8TEV_PTFilter::ReadData()
{
  fstream f1;
  fstream f2;
  fstream f3;

  //Full breakdown of statistical uncertainties
  stringstream datafile_sys("");
  datafile_sys << dataPath()
	   << "rawdata/ATLAS_WJET_8TEV/HEPData-ins1635273-v1-Table_13.csv";
  f1.open(datafile_sys.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile_sys.str() << endl;
      exit(-1);
    }

  //Correlation matrix of statistical ucnertainties
  stringstream datafile_stat("");
  datafile_stat << dataPath()
		<< "rawdata/ATLAS_WJET_8TEV/HEPData-ins1635273-v1-Table_15.csv";

  f2.open(datafile_stat.str().c_str(), ios::in);

  if (f2.fail())
    {
      cerr << "Error opening data file " << datafile_stat.str() << endl;
      exit(-1);
    }

  //Non-perturbative corrections
  stringstream datafile_np("");
  datafile_np << dataPath()
		<< "rawdata/ATLAS_WJET_8TEV/HEPData-ins1635273-v1-Table_67.csv";

  f3.open(datafile_np.str().c_str(), ios::in);

  if (f3.fail())
    {
      cerr << "Error opening data file " << datafile_np.str() << endl;
      exit(-1);
    }

  //Read central value
  string line;
  for(int i=0; i<41; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<16; i++)
    {
      getline(f2,line);
    }

  for(int i=0; i<42; i++)
    {
      getline(f3,line);
    }

  double** corrmat = new double*[fNData];
  double** syscor  = new double*[fNData];
  double** sysR    = new double*[fNData];
  double** sysL    = new double*[fNData];
  double npcorr[fNData];
  const int nrealsys = 54;
  double ddum;
  char comma;

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      lstream >> fKin1[i] >> comma
	      >> ddum     >> comma
	      >> ddum     >> comma
	      >> fData[i] >> comma
	      >> fStat[i] >> comma
	      >> ddum     >> comma;

      fKin2[i] = 0;
      fKin3[i] = 8000; //GeV

      sysR[i]    = new double[nrealsys];
      sysL[i]    = new double[nrealsys];

      for(int k=0; k<nrealsys; k++)
	{ 
	  lstream >> sysR[i][k] >> comma >> sysL[i][k] >> comma;
	}

      corrmat[i] = new double[fNData];
      syscor[i]  = new double[fNData];

      for(int j=0; j<fNData; j++)
	{
	  getline(f2,line);
	  istringstream istream(line);
	  istream >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> corrmat[i][j];
	}

      getline(f3,line);
      istringstream kstream(line);
      kstream >> ddum >> comma
	      >> ddum >> comma
	      >> ddum >> comma
	      >> npcorr[i];

    }

  //Generate artificial systematics
  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  corrmat[i][j] *= fStat[i]*fStat[j];
	}
    }
  
  if(!genArtSys(fNData,corrmat,syscor))
    {
      throw runtime_error("Couldn't generate artificial systematics for " + fSetName);
    }

  for(int i=0; i<fNData; i++)
    {
      fStat[i]=0.;
      for(int j=0; j<fNData; j++)
	{
	  fSys[i][j].add  = syscor[i][j];
	  fSys[i][j].mult = fSys[i][j].add*1e2/fData[i];
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = "CORR";
	}
    }

  //Real systematics (correlated across bins and W=/W- distributions)
  for(int i=0; i<fNData; i++)
    {
      
      for(int k=0; k<50; k++)
	{
	  sysR[i][k] /= sqrt(2.);
	  sysL[i][k] /= sqrt(2.);

	  double tmp1, tmp2;
	  tmp1=sysR[i][k];
	  tmp2=sysL[i][k];
	  
	  //Case 1: sysL and sysR are both negative
	  if(tmp1<0 && tmp2<0)
	    {
	      if(tmp2<tmp1 || tmp2==tmp1)
		{
		  sysR[i][k] = 0.0;
		  sysL[i][k] = tmp2;
		}
	      if(tmp2>tmp1)
		{
		  sysR[i][k] = 0.0;
		  sysL[i][k] = tmp1;
		}
	    }

	  //Case 2: sysL and sysR are both positive
	  if(tmp1>0.0 && tmp2>0.0)
	    {
	      if(tmp1>tmp2 || tmp1==tmp2)
		{
		  sysR[i][k] = tmp1;
		  sysL[i][k] = 0.0;
		}
	      if(tmp1<tmp2)
		{
		  sysR[i][k] = tmp2;
		  sysL[i][k] = 0.0;
		}
	    }

	  //Case3: sys1 is negative and sys2 is positive
	  if(tmp1<0.0 && tmp2>0.0)
	    {
	      sysR[i][k] = tmp2;
	      sysL[i][k] = tmp1;
	    }
	  
	  fSys[i][2*k+fNData].add = sysR[i][k];
	  fSys[i][2*k+fNData].mult = fSys[i][2*k+fNData].add*1e2/fData[i];
	  fSys[i][2*k+fNData].type = MULT;
	  ostringstream sysnameR;
	  sysnameR << "ATLASWJ" << 2*k;
	  fSys[i][2*k+fNData].name = sysnameR.str();
	  
	  fSys[i][2*k+1+fNData].add = sysL[i][k];
	  fSys[i][2*k+1+fNData].mult = fSys[i][2*k+1+fNData].add*1e2/fData[i];
	  fSys[i][2*k+1+fNData].type = MULT;
	  ostringstream sysnameL;
	  sysnameL << "ATLASWJ" << 2*k+1;
	  fSys[i][2*k+1+fNData].name = sysnameL.str();

	}

      //Luminosity uncertainty
      fSys[i][2*50+fNData].add  = sysR[i][50];
      fSys[i][2*50+fNData].mult = fSys[i][2*50+fNData].add*1e2/fData[i];
      fSys[i][2*50+fNData].type = MULT;
      fSys[i][2*50+fNData].name = "ATLASLUMI12";

      //Uncorrelated unfolding uncertainties
      for(int k=51; k<nrealsys; k++)
	{
	  sysR[i][k] /= sqrt(2.);
	  sysL[i][k] /= sqrt(2.);
	  
	  double tmp1, tmp2;
	  tmp1=sysR[i][k];
	  tmp2=sysL[i][k];
	  
	  //Case 1: sysL and sysR are both negative
	  if(tmp1<0 && tmp2<0)
	    {
	      if(tmp2<tmp1 || tmp2==tmp1)
		{
		  sysR[i][k] = 0.0;
		  sysL[i][k] = tmp2;
		}
	      if(tmp2>tmp1)
		{
		  sysR[i][k] = 0.0;
		  sysL[i][k] = tmp1;
		}
	    }

	  //Case 2: sysL and sysR are both positive
	  if(tmp1>0.0 && tmp2>0.0)
	    {
	      if(tmp1>tmp2 || tmp1==tmp2)
		{
		  sysR[i][k] = tmp1;
		  sysL[i][k] = 0.0;
		}
	      if(tmp1<tmp2)
		{
		  sysR[i][k] = tmp2;
		  sysL[i][k] = 0.0;
		}
	    }
	  
	  //Case3: sys1 is negative and sys2 is positive
	  if(tmp1<0.0 && tmp2>0.0)
	    {
	      sysR[i][k] = tmp2;
	      sysL[i][k] = tmp1;
	    }
	  	  
	  fSys[i][2*k+fNData-1].add  = sysR[i][k];
	  fSys[i][2*k+fNData-1].mult = fSys[i][2*k+fNData-1].add*1e2/fData[i];
	  fSys[i][2*k+fNData-1].type = MULT;
	  fSys[i][2*k+fNData-1].name = "UNCORR";
	  
	  fSys[i][2*k+fNData].add  = sysL[i][k];
	  fSys[i][2*k+fNData].mult = fSys[i][2*k+fNData].add*1e2/fData[i];
	  fSys[i][2*k+fNData].type = MULT;
	  fSys[i][2*k+fNData].name = "UNCORR";
	}
      
      //Non-perturbative corrections
      fSys[i][123].add  = fData[i]*(1. - npcorr[i]);
      fSys[i][123].mult = fSys[i][123].add*1e2/fData[i];
      fSys[i][123].type = MULT;
      fSys[i][123].name = "SKIP";

    }

 // Clean-up
  for (int i=0; i<fNData; i++)
    delete[] syscor[i];
  
  delete[] syscor;

  for(int i=0; i<fNData; i++)
    delete[] corrmat[i];

  delete[] corrmat;

  for(int i=0; i<fNData; i++)
    delete[] sysR[i];

  delete[] sysR;

  for(int i=0; i<fNData; i++)
    delete[] sysL[i];

  delete[] sysL;
  
  f1.close();
  f2.close();

} 

//3)W+ distribution differential in the transverse momentum of the leading jet

void ATLAS_WP_JET_8TEV_PTJFilter::ReadData()
{
  fstream f1;
  fstream f2;
  fstream f3;

  //Full breakdown of statistical uncertainties
  stringstream datafile_sys("");
  datafile_sys << dataPath()
	   << "rawdata/ATLAS_WJET_8TEV/HEPData-ins1635273-v1-Table_18.csv";
  f1.open(datafile_sys.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile_sys.str() << endl;
      exit(-1);
    }

  //Correlation matrix of statistical ucnertainties
  stringstream datafile_stat("");
  datafile_stat << dataPath()
		<< "rawdata/ATLAS_WJET_8TEV/HEPData-ins1635273-v1-Table_19.csv";

  f2.open(datafile_stat.str().c_str(), ios::in);

  if (f2.fail())
    {
      cerr << "Error opening data file " << datafile_stat.str() << endl;
      exit(-1);
    }

  //Non-perturbative corrections
  stringstream datafile_np("");
  datafile_np << dataPath()
	      << "rawdata/ATLAS_WJET_8TEV/HEPData-ins1635273-v1-Table_69.csv";

  f3.open(datafile_np.str().c_str(), ios::in);

  if (f3.fail())
    {
      cerr << "Error opening data file " << datafile_np.str() << endl;
      exit(-1);
    }

  //Read central value
  string line;
  for(int i=0; i<16; i++)
    {
      getline(f1,line);
      getline(f2,line);
      getline(f3,line);
    }

  double** corrmat = new double*[fNData];
  double** syscor  = new double*[fNData];
  double** sysR    = new double*[fNData];
  double** sysL    = new double*[fNData];
  double npcorr[fNData];
  const int nrealsys = 54;
  double ddum;
  char comma;

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      lstream >> fKin1[i] >> comma
	      >> ddum     >> comma
	      >> ddum     >> comma
	      >> fData[i] >> comma
	      >> fStat[i] >> comma
	      >> ddum     >> comma;

      fKin2[i] = 0;
      fKin3[i] = 8000; //GeV

      sysR[i]    = new double[nrealsys];
      sysL[i]    = new double[nrealsys];

      for(int k=0; k<nrealsys; k++)
	{ 
	  lstream >> sysR[i][k] >> comma >> sysL[i][k] >> comma;
	}

      corrmat[i] = new double[fNData];
      syscor[i]  = new double[fNData];

      for(int j=0; j<fNData; j++)
	{
	  getline(f2,line);
	  istringstream istream(line);
	  istream >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> corrmat[i][j];
	}

      getline(f3,line);
      istringstream kstream(line);
      kstream >> ddum >> comma
	      >> ddum >> comma
	      >> ddum >> comma
	      >> npcorr[i];

    }

  //Generate artificial systematics
  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  corrmat[i][j] *= fStat[i]*fStat[j];
	}
    }
  
  if(!genArtSys(fNData,corrmat,syscor))
    {
      throw runtime_error("Couldn't generate artificial systematics for " + fSetName);
    }

  for(int i=0; i<fNData; i++)
    {
      fStat[i]=0.;
      for(int j=0; j<fNData; j++)
	{
	  fSys[i][j].add  = syscor[i][j];
	  fSys[i][j].mult = fSys[i][j].add*1e2/fData[i];
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = "CORR";
	}
    }

  //Real systematics (correlated across bins and W=/W- distributions)
  for(int i=0; i<fNData; i++)
    {
      
      for(int k=0; k<50; k++)
	{
	  sysR[i][k] /= sqrt(2.);
	  sysL[i][k] /= sqrt(2.);

	  double tmp1, tmp2;
	  tmp1=sysR[i][k];
	  tmp2=sysL[i][k];
	  
	  //Case 1: sysL and sysR are both negative
	  if(tmp1<0 && tmp2<0)
	    {
	      if(tmp2<tmp1 || tmp2==tmp1)
		{
		  sysR[i][k] = 0.0;
		  sysL[i][k] = tmp2;
		}
	      if(tmp2>tmp1)
		{
		  sysR[i][k] = 0.0;
		  sysL[i][k] = tmp1;
		}
	    }

	  //Case 2: sysL and sysR are both positive
	  if(tmp1>0.0 && tmp2>0.0)
	    {
	      if(tmp1>tmp2 || tmp1==tmp2)
		{
		  sysR[i][k] = tmp1;
		  sysL[i][k] = 0.0;
		}
	      if(tmp1<tmp2)
		{
		  sysR[i][k] = tmp2;
		  sysL[i][k] = 0.0;
		}
	    }

	  //Case3: sys1 is negative and sys2 is positive
	  if(tmp1<0.0 && tmp2>0.0)
	    {
	      sysR[i][k] = tmp2;
	      sysL[i][k] = tmp1;
	    }
	  
	  fSys[i][2*k+fNData].add = sysR[i][k];
	  fSys[i][2*k+fNData].mult = fSys[i][2*k+fNData].add*1e2/fData[i];
	  fSys[i][2*k+fNData].type = MULT;
	  ostringstream sysnameR;
	  sysnameR << "ATLASWJ" << 2*k;
	  fSys[i][2*k+fNData].name = sysnameR.str();
	  
	  fSys[i][2*k+1+fNData].add = sysL[i][k];
	  fSys[i][2*k+1+fNData].mult = fSys[i][2*k+1+fNData].add*1e2/fData[i];
	  fSys[i][2*k+1+fNData].type = MULT;
	  ostringstream sysnameL;
	  sysnameL << "ATLASWJ" << 2*k+1;
	  fSys[i][2*k+1+fNData].name = sysnameL.str();

	}

      //Luminosity uncertainty
      fSys[i][2*50+fNData].add  = sysR[i][50];
      fSys[i][2*50+fNData].mult = fSys[i][2*50+fNData].add*1e2/fData[i];
      fSys[i][2*50+fNData].type = MULT;
      fSys[i][2*50+fNData].name = "ATLASLUMI12";

      //Uncorrelated unfolding uncertainties
      for(int k=51; k<nrealsys; k++)
	{
	  sysR[i][k] /= sqrt(2.);
	  sysL[i][k] /= sqrt(2.);
	  
	  double tmp1, tmp2;
	  tmp1=sysR[i][k];
	  tmp2=sysL[i][k];
	  
	  //Case 1: sysL and sysR are both negative
	  if(tmp1<0 && tmp2<0)
	    {
	      if(tmp2<tmp1 || tmp2==tmp1)
		{
		  sysR[i][k] = 0.0;
		  sysL[i][k] = tmp2;
		}
	      if(tmp2>tmp1)
		{
		  sysR[i][k] = 0.0;
		  sysL[i][k] = tmp1;
		}
	    }

	  //Case 2: sysL and sysR are both positive
	  if(tmp1>0.0 && tmp2>0.0)
	    {
	      if(tmp1>tmp2 || tmp1==tmp2)
		{
		  sysR[i][k] = tmp1;
		  sysL[i][k] = 0.0;
		}
	      if(tmp1<tmp2)
		{
		  sysR[i][k] = tmp2;
		  sysL[i][k] = 0.0;
		}
	    }
	  
	  //Case3: sys1 is negative and sys2 is positive
	  if(tmp1<0.0 && tmp2>0.0)
	    {
	      sysR[i][k] = tmp2;
	      sysL[i][k] = tmp1;
	    }
	  	  
	  fSys[i][2*k+fNData-1].add  = sysR[i][k];
	  fSys[i][2*k+fNData-1].mult = fSys[i][2*k+fNData-1].add*1e2/fData[i];
	  fSys[i][2*k+fNData-1].type = MULT;
	  fSys[i][2*k+fNData-1].name = "UNCORR";
	  
	  fSys[i][2*k+fNData].add  = sysL[i][k];
	  fSys[i][2*k+fNData].mult = fSys[i][2*k+fNData].add*1e2/fData[i];
	  fSys[i][2*k+fNData].type = MULT;
	  fSys[i][2*k+fNData].name = "UNCORR";
	}
      
      //Non-perturbative corrections
      fSys[i][129].add  = fData[i]*(1. - npcorr[i]);
      fSys[i][129].mult = fSys[i][123].add*1e2/fData[i];
      fSys[i][129].type = MULT;
      fSys[i][129].name = "SKIP";

    }

 // Clean-up
  for (int i=0; i<fNData; i++)
    delete[] syscor[i];
  
  delete[] syscor;

  for(int i=0; i<fNData; i++)
    delete[] corrmat[i];

  delete[] corrmat;

  for(int i=0; i<fNData; i++)
    delete[] sysR[i];

  delete[] sysR;

  for(int i=0; i<fNData; i++)
    delete[] sysL[i];

  delete[] sysL;
  
  f1.close();
  f2.close();

} 

//4)W- distribution differential in the transverse momentum of the leading jet

void ATLAS_WM_JET_8TEV_PTJFilter::ReadData()
{
  fstream f1;
  fstream f2;
  fstream f3;

  //Full breakdown of statistical uncertainties
  stringstream datafile_sys("");
  datafile_sys << dataPath()
	   << "rawdata/ATLAS_WJET_8TEV/HEPData-ins1635273-v1-Table_18.csv";
  f1.open(datafile_sys.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile_sys.str() << endl;
      exit(-1);
    }

  //Correlation matrix of statistical ucnertainties
  stringstream datafile_stat("");
  datafile_stat << dataPath()
		<< "rawdata/ATLAS_WJET_8TEV/HEPData-ins1635273-v1-Table_20.csv";

  f2.open(datafile_stat.str().c_str(), ios::in);

  if (f2.fail())
    {
      cerr << "Error opening data file " << datafile_stat.str() << endl;
      exit(-1);
    }

  //Non-perturbative corrections
  stringstream datafile_np("");
  datafile_np << dataPath()
		<< "rawdata/ATLAS_WJET_8TEV/HEPData-ins1635273-v1-Table_69.csv";

  f3.open(datafile_np.str().c_str(), ios::in);

  if (f3.fail())
    {
      cerr << "Error opening data file " << datafile_np.str() << endl;
      exit(-1);
    }

  //Read central value
  string line;
  for(int i=0; i<47; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<16; i++)
    {
      getline(f2,line);
    }

  for(int i=0; i<48; i++)
    {
      getline(f3,line);
    }

  double** corrmat = new double*[fNData];
  double** syscor  = new double*[fNData];
  double** sysR    = new double*[fNData];
  double** sysL    = new double*[fNData];
  double npcorr[fNData];
  const int nrealsys = 54;
  double ddum;
  char comma;

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      lstream >> fKin1[i] >> comma
	      >> ddum     >> comma
	      >> ddum     >> comma
	      >> fData[i] >> comma
	      >> fStat[i] >> comma
	      >> ddum     >> comma;

      fKin2[i] = 0;
      fKin3[i] = 8000; //GeV

      sysR[i]    = new double[nrealsys];
      sysL[i]    = new double[nrealsys];

      for(int k=0; k<nrealsys; k++)
	{ 
	  lstream >> sysR[i][k] >> comma >> sysL[i][k] >> comma;
	}

      corrmat[i] = new double[fNData];
      syscor[i]  = new double[fNData];

      for(int j=0; j<fNData; j++)
	{
	  getline(f2,line);
	  istringstream istream(line);
	  istream >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> corrmat[i][j];
	}

      getline(f3,line);
      istringstream kstream(line);
      kstream >> ddum >> comma
	      >> ddum >> comma
	      >> ddum >> comma
	      >> npcorr[i];

    }

  //Generate artificial systematics
  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  corrmat[i][j] *= fStat[i]*fStat[j];
	}
    }
  
  if(!genArtSys(fNData,corrmat,syscor))
    {
      throw runtime_error("Couldn't generate artificial systematics for " + fSetName);
    }

  for(int i=0; i<fNData; i++)
    {
      fStat[i]=0.;
      for(int j=0; j<fNData; j++)
	{
	  fSys[i][j].add  = syscor[i][j];
	  fSys[i][j].mult = fSys[i][j].add*1e2/fData[i];
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = "CORR";
	}
    }

  //Real systematics (correlated across bins and W=/W- distributions)
  for(int i=0; i<fNData; i++)
    {
      
      for(int k=0; k<50; k++)
	{
	  sysR[i][k] /= sqrt(2.);
	  sysL[i][k] /= sqrt(2.);

	  double tmp1, tmp2;
	  tmp1=sysR[i][k];
	  tmp2=sysL[i][k];
	  
	  //Case 1: sysL and sysR are both negative
	  if(tmp1<0 && tmp2<0)
	    {
	      if(tmp2<tmp1 || tmp2==tmp1)
		{
		  sysR[i][k] = 0.0;
		  sysL[i][k] = tmp2;
		}
	      if(tmp2>tmp1)
		{
		  sysR[i][k] = 0.0;
		  sysL[i][k] = tmp1;
		}
	    }

	  //Case 2: sysL and sysR are both positive
	  if(tmp1>0.0 && tmp2>0.0)
	    {
	      if(tmp1>tmp2 || tmp1==tmp2)
		{
		  sysR[i][k] = tmp1;
		  sysL[i][k] = 0.0;
		}
	      if(tmp1<tmp2)
		{
		  sysR[i][k] = tmp2;
		  sysL[i][k] = 0.0;
		}
	    }

	  //Case3: sys1 is negative and sys2 is positive
	  if(tmp1<0.0 && tmp2>0.0)
	    {
	      sysR[i][k] = tmp2;
	      sysL[i][k] = tmp1;
	    }
	  
	  fSys[i][2*k+fNData].add = sysR[i][k];
	  fSys[i][2*k+fNData].mult = fSys[i][2*k+fNData].add*1e2/fData[i];
	  fSys[i][2*k+fNData].type = MULT;
	  ostringstream sysnameR;
	  sysnameR << "ATLASWJ" << 2*k;
	  fSys[i][2*k+fNData].name = sysnameR.str();
	  
	  fSys[i][2*k+1+fNData].add = sysL[i][k];
	  fSys[i][2*k+1+fNData].mult = fSys[i][2*k+1+fNData].add*1e2/fData[i];
	  fSys[i][2*k+1+fNData].type = MULT;
	  ostringstream sysnameL;
	  sysnameL << "ATLASWJ" << 2*k+1;
	  fSys[i][2*k+1+fNData].name = sysnameL.str();

	}

      //Luminosity uncertainty
      fSys[i][2*50+fNData].add  = sysR[i][50];
      fSys[i][2*50+fNData].mult = fSys[i][2*50+fNData].add*1e2/fData[i];
      fSys[i][2*50+fNData].type = MULT;
      fSys[i][2*50+fNData].name = "ATLASLUMI12";

      //Uncorrelated unfolding uncertainties
      for(int k=51; k<nrealsys; k++)
	{
	  sysR[i][k] /= sqrt(2.);
	  sysL[i][k] /= sqrt(2.);
	  
	  double tmp1, tmp2;
	  tmp1=sysR[i][k];
	  tmp2=sysL[i][k];
	  
	  //Case 1: sysL and sysR are both negative
	  if(tmp1<0 && tmp2<0)
	    {
	      if(tmp2<tmp1 || tmp2==tmp1)
		{
		  sysR[i][k] = 0.0;
		  sysL[i][k] = tmp2;
		}
	      if(tmp2>tmp1)
		{
		  sysR[i][k] = 0.0;
		  sysL[i][k] = tmp1;
		}
	    }

	  //Case 2: sysL and sysR are both positive
	  if(tmp1>0.0 && tmp2>0.0)
	    {
	      if(tmp1>tmp2 || tmp1==tmp2)
		{
		  sysR[i][k] = tmp1;
		  sysL[i][k] = 0.0;
		}
	      if(tmp1<tmp2)
		{
		  sysR[i][k] = tmp2;
		  sysL[i][k] = 0.0;
		}
	    }
	  
	  //Case3: sys1 is negative and sys2 is positive
	  if(tmp1<0.0 && tmp2>0.0)
	    {
	      sysR[i][k] = tmp2;
	      sysL[i][k] = tmp1;
	    }
	  	  
	  fSys[i][2*k+fNData-1].add  = sysR[i][k];
	  fSys[i][2*k+fNData-1].mult = fSys[i][2*k+fNData-1].add*1e2/fData[i];
	  fSys[i][2*k+fNData-1].type = MULT;
	  fSys[i][2*k+fNData-1].name = "UNCORR";
	  
	  fSys[i][2*k+fNData].add  = sysL[i][k];
	  fSys[i][2*k+fNData].mult = fSys[i][2*k+fNData].add*1e2/fData[i];
	  fSys[i][2*k+fNData].type = MULT;
	  fSys[i][2*k+fNData].name = "UNCORR";
	}
      
      //Non-perturbative corrections
      fSys[i][129].add  = fData[i]*(1. - npcorr[i]);
      fSys[i][129].mult = fSys[i][123].add*1e2/fData[i];
      fSys[i][129].type = MULT;
      fSys[i][129].name = "SKIP";

    }

 // Clean-up
  for (int i=0; i<fNData; i++)
    delete[] syscor[i];
  
  delete[] syscor;

  for(int i=0; i<fNData; i++)
    delete[] corrmat[i];

  delete[] corrmat;

  for(int i=0; i<fNData; i++)
    delete[] sysR[i];

  delete[] sysR;

  for(int i=0; i<fNData; i++)
    delete[] sysL[i];

  delete[] sysL;
  
  f1.close();
  f2.close();

} 




















































