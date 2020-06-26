/*
Reference:
   arXiv:     [1905.04242]
   hepdata:   https://www.hepdata.net/record/ins1734263
   published: Eur.Phys.J.C 79 (2019) 10, 884
Description:
   Measurements of fiducial and differential cross-sections for W+ W- 
   production in proton–proton collisions at a centre-of-mass energy of 13 TeV 
  with the ATLAS experiment at the Large Hadron Collider using data 
  corresponding to an integrated luminosity of 36.1 fb-1. Events with one 
  electron and one muon are selected. The fducial cross section and four
  distributions differential in the invariant mass of the lepton pair, the pT
  of the lepton pair, the leading pT and the rapidity of the lepton pair are
  implemented. The information is taken from Tables 7-8, 10-11, 4-5 and 13-14
  from the hepdata entry.
*/

#include "ATLAS_WW_13TEV.h"

//ATLAS_WW_13TEV_memu: combined mass spectrum
void ATLAS_WW_13TEV_memuFilter::ReadData()
{
  fstream f1;
  fstream f2;

  //Central values and breakdown of ucnertainties
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/ATLAS_WW_13TEV/Table7.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Statistical correlation matrix
  stringstream corrfile("");
  corrfile << dataPath()
	   << "rawdata/ATLAS_WW_13TEV/Table8.csv";
  f2.open(corrfile.str().c_str(), ios::in);

  if (f2.fail())
    {
      cerr << "Error opening data file " << corrfile.str() << endl;
      exit(-1);
    }
 
  //Read central values and uncertainties
  string line;
  const int realsys=11;

  double* Stat = new double[fNData];
  double** corrmat = new double*[fNData];
  double** syscor  = new double*[fNData];

  for(int i=0; i<21; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<21; i++)
    {
      getline(f2,line);
    }

  for(int i=0; i<fNData; i++)
    {
      double ddum;
      char comma;
      getline(f1,line);
      istringstream lstream(line);

      fKin2[i] = 0.;
      fKin3[i] = 13000.;    //c.m. energy
      lstream >> fKin1[i]        >> comma
	      >> ddum            >> comma
	      >> ddum            >> comma
	      >> fData[i]        >> comma
	      >> Stat[i];

      fStat[i] = 0.;

      for(int j=0; j<realsys; j++)
	{
	  lstream >> comma >> fSys[i][j].add
		  >> comma >> ddum;
	  fSys[i][j].mult = fSys[i][j].add/fData[i]*100.;
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = "CORR";
	}

      //Exception
      fSys[i][5].type = ADD;           //statistical background uncertainty
      fSys[i][7].name = "ATLASLUMI13"; //luminosity uncertainty

      corrmat[i] = new double[fNData];
      syscor[i]  = new double[fNData];

      for(int j=0; j<fNData; j++)
	{
	  getline(f2,line);
	  istringstream kstream(line);
	  kstream >> ddum >> comma 
		  >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> corrmat[i][j];
	}
    }

  //Generate covariance matrix from correlation matrix
    for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  corrmat[i][j] = corrmat[i][j]*Stat[i]*Stat[j];
	}
    }
    
    //Generate artificial systematics from covariance matrix
    if(!genArtSys(fNData,corrmat,syscor))
      {
	throw runtime_error("Couldn't generate artificial systematics for " + fSetName);
      }
    
    for(int i=0; i<fNData; i++)
      {
	for(int j=realsys; j<fNSys; j++)
	  {
	    fSys[i][j].add  = syscor[i][j];
	    fSys[i][j].mult = fSys[i][j].add*1e2/fData[i];
	    fSys[i][j].type = ADD;
	    fSys[i][j].name = "CORR";
	  }
      } 
    
    f1.close();
    f2.close();

}

//ATLAS_WW_13TEV_pTemu: combined pT spectrum
void ATLAS_WW_13TEV_pTemuFilter::ReadData()
{
  fstream f1;
  fstream f2;

  //Central values and breakdown of ucnertainties
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/ATLAS_WW_13TEV/Table10.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Statistical correlation matrix
  stringstream corrfile("");
  corrfile << dataPath()
	   << "rawdata/ATLAS_WW_13TEV/Table11.csv";
  f2.open(corrfile.str().c_str(), ios::in);

  if (f2.fail())
    {
      cerr << "Error opening data file " << corrfile.str() << endl;
      exit(-1);
    }
 
  //Read central values and uncertainties
  string line;
  const int realsys=11;

  double* Stat = new double[fNData];
  double** corrmat = new double*[fNData];
  double** syscor  = new double*[fNData];

  for(int i=0; i<21; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<21; i++)
    {
      getline(f2,line);
    }

  for(int i=0; i<fNData; i++)
    {
      double ddum;
      char comma;
      getline(f1,line);
      istringstream lstream(line);

      fKin2[i] = 0.;
      fKin3[i] = 13000.;    //c.m. energy
      lstream >> fKin1[i]        >> comma
	      >> ddum            >> comma
	      >> ddum            >> comma
	      >> fData[i]        >> comma
	      >> Stat[i];

      fStat[i] = 0.;

      for(int j=0; j<realsys; j++)
	{
	  lstream >> comma >> fSys[i][j].add
		  >> comma >> ddum;
	  fSys[i][j].mult = fSys[i][j].add/fData[i]*100.;
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = "CORR";
	}

      //Exception
      fSys[i][5].type = ADD;           //statistical background uncertainty
      fSys[i][7].name = "ATLASLUMI13"; //luminosity uncertainty

      corrmat[i] = new double[fNData];
      syscor[i]  = new double[fNData];

      for(int j=0; j<fNData; j++)
	{
	  getline(f2,line);
	  istringstream kstream(line);
	  kstream >> ddum >> comma 
		  >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> corrmat[i][j];
	}
    }

  //Generate covariance matrix from correlation matrix
    for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  corrmat[i][j] = corrmat[i][j]*Stat[i]*Stat[j];
	}
    }
    
    //Generate artificial systematics from covariance matrix
    if(!genArtSys(fNData,corrmat,syscor))
      {
	throw runtime_error("Couldn't generate artificial systematics for " + fSetName);
      }
    
    for(int i=0; i<fNData; i++)
      {
	for(int j=realsys; j<fNSys; j++)
	  {
	    fSys[i][j].add  = syscor[i][j];
	    fSys[i][j].mult = fSys[i][j].add*1e2/fData[i];
	    fSys[i][j].type = ADD;
	    fSys[i][j].name = "CORR";
	  }
      } 
    
    f1.close();
    f2.close();

}

//ATLAS_WW_13TEV_pTleademu: combined pT leading spectrum
void ATLAS_WW_13TEV_pTleadFilter::ReadData()
{
  fstream f1;
  fstream f2;

  //Central values and breakdown of ucnertainties
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/ATLAS_WW_13TEV/Table4.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Statistical correlation matrix
  stringstream corrfile("");
  corrfile << dataPath()
	   << "rawdata/ATLAS_WW_13TEV/Table5.csv";
  f2.open(corrfile.str().c_str(), ios::in);

  if (f2.fail())
    {
      cerr << "Error opening data file " << corrfile.str() << endl;
      exit(-1);
    }
 
  //Read central values and uncertainties
  string line;
  const int realsys=11;

  double* Stat = new double[fNData];
  double** corrmat = new double*[fNData];
  double** syscor  = new double*[fNData];

  for(int i=0; i<21; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<21; i++)
    {
      getline(f2,line);
    }

  for(int i=0; i<fNData; i++)
    {
      double ddum;
      char comma;
      getline(f1,line);
      istringstream lstream(line);

      fKin2[i] = 0.;
      fKin3[i] = 13000.;    //c.m. energy
      lstream >> fKin1[i]        >> comma
	      >> ddum            >> comma
	      >> ddum            >> comma
	      >> fData[i]        >> comma
	      >> Stat[i];

      fStat[i] = 0.;

      for(int j=0; j<realsys; j++)
	{
	  lstream >> comma >> fSys[i][j].add
		  >> comma >> ddum;
	  fSys[i][j].mult = fSys[i][j].add/fData[i]*100.;
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = "CORR";
	}

      //Exception
      fSys[i][5].type = ADD;           //statistical background uncertainty
      fSys[i][7].name = "ATLASLUMI13"; //luminosity uncertainty

      corrmat[i] = new double[fNData];
      syscor[i]  = new double[fNData];

      for(int j=0; j<fNData; j++)
	{
	  getline(f2,line);
	  istringstream kstream(line);
	  kstream >> ddum >> comma 
		  >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> corrmat[i][j];
	}
    }

  //Generate covariance matrix from correlation matrix
    for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  corrmat[i][j] = corrmat[i][j]*Stat[i]*Stat[j];
	}
    }
    
    //Generate artificial systematics from covariance matrix
    if(!genArtSys(fNData,corrmat,syscor))
      {
	throw runtime_error("Couldn't generate artificial systematics for " + fSetName);
      }
    
    for(int i=0; i<fNData; i++)
      {
	for(int j=realsys; j<fNSys; j++)
	  {
	    fSys[i][j].add  = syscor[i][j];
	    fSys[i][j].mult = fSys[i][j].add*1e2/fData[i];
	    fSys[i][j].type = ADD;
	    fSys[i][j].name = "CORR";
	  }
      } 
    
    f1.close();
    f2.close();

}

//ATLAS_WW_13TEV_yemu: combined rapidity spectrum
void ATLAS_WW_13TEV_yemuFilter::ReadData()
{
  fstream f1;
  fstream f2;

  //Central values and breakdown of ucnertainties
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/ATLAS_WW_13TEV/Table13.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Statistical correlation matrix
  stringstream corrfile("");
  corrfile << dataPath()
	   << "rawdata/ATLAS_WW_13TEV/Table14.csv";
  f2.open(corrfile.str().c_str(), ios::in);

  if (f2.fail())
    {
      cerr << "Error opening data file " << corrfile.str() << endl;
      exit(-1);
    }
 
  //Read central values and uncertainties
  string line;
  const int realsys=11;

  double* Stat = new double[fNData];
  double** corrmat = new double*[fNData];
  double** syscor  = new double*[fNData];

  for(int i=0; i<21; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<21; i++)
    {
      getline(f2,line);
    }

  for(int i=0; i<fNData; i++)
    {
      double ddum;
      char comma;
      getline(f1,line);
      istringstream lstream(line);

      fKin2[i] = 0.;
      fKin3[i] = 13000.;    //c.m. energy
      lstream >> fKin1[i]        >> comma
	      >> ddum            >> comma
	      >> ddum            >> comma
	      >> fData[i]        >> comma
	      >> Stat[i];

      fStat[i] = 0.;

      for(int j=0; j<realsys; j++)
	{
	  lstream >> comma >> fSys[i][j].add
		  >> comma >> ddum;
	  fSys[i][j].mult = fSys[i][j].add/fData[i]*100.;
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = "CORR";
	}

      //Exception
      fSys[i][5].type = ADD;           //statistical background uncertainty
      fSys[i][7].name = "ATLASLUMI13"; //luminosity uncertainty

      corrmat[i] = new double[fNData];
      syscor[i]  = new double[fNData];

      for(int j=0; j<fNData; j++)
	{
	  getline(f2,line);
	  istringstream kstream(line);
	  kstream >> ddum >> comma 
		  >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> ddum >> comma
		  >> corrmat[i][j];
	}
    }

  //Generate covariance matrix from correlation matrix
    for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  corrmat[i][j] = corrmat[i][j]*Stat[i]*Stat[j];
	}
    }
    
    //Generate artificial systematics from covariance matrix
    if(!genArtSys(fNData,corrmat,syscor))
      {
	throw runtime_error("Couldn't generate artificial systematics for " + fSetName);
      }
    
    for(int i=0; i<fNData; i++)
      {
	for(int j=realsys; j<fNSys; j++)
	  {
	    fSys[i][j].add  = syscor[i][j];
	    fSys[i][j].mult = fSys[i][j].add*1e2/fData[i];
	    fSys[i][j].type = ADD;
	    fSys[i][j].name = "CORR";
	  }
      } 
    
    f1.close();
    f2.close();

}

//ATLAS_WW_13TEV_totWW: combined mass spectrum
void ATLAS_WW_13TEV_totWWFilter::ReadData()
{
  fKin1[0] = 0.0;
  fKin2[0] = 0.0;
  fKin3[0] = 13000.;
  fData[0] = 379.1;        //fb
  fStat[0] = 5.0;          //statistical uncertainty
  fSys[0][0].add = 25.4;   //total systematic uncertainty
  fSys[0][0].mult = fSys[0][0].add/fData[0] * 100.;
  fSys[0][0].type = MULT;
  fSys[0][0].name = "CORR";
  fSys[0][1].add = 8.0;    //luminosity uncertainty
  fSys[0][1].mult = fSys[0][0].add/fData[0] * 100.;
  fSys[0][1].type = MULT;
  fSys[0][1].name = "ATLASLUMI13";
}
