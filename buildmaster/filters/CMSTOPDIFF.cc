/*
Experiment: CERN-LHC-CMS (CMS)
Preprinted as CERN-PH-EP-2015-117
Preprinted as CMS-TOP-12-028
Archived as: ARXIV:1505.04480
Published in Eur.Phys.J. C75 (2015) 11, 542

Record in: INSPIRE
Record in: CERN Document Server
Record in: HEPData 

Description of the measurement 
The normalized differential cross section for top quark pair (ttbar) 
production is measured in pp collisions at a centre-of-mass energy 
of 8 TeV at the CERN LHC using the CMS detector in data corresponding 
to an integrated luminosity of 19.7 inverse femtobarns. 
The measurements are performed in the lepton+jets (e/mu+jets) and in the 
dilepton (e e, mu mu, and e mu) decay channels. 
The ttbar cross section is measured as a function of the kinematic 
properties of the charged leptons, the jets associated to b quarks, 
the top quarks, and the ttbar system. 

Description of the buildmaster implementation
Normalized cross sections for the distributions (lepton+jets channel)
differential in the following variables are implemented:
1) top quark transverse momentum;       
2) top quark pair transverse momentum;  
3) top quark rapidity;                  
4) top quark pair rapidity;             
5) top quark pair invariant mass.       

Raw data and full breakdown of systematic uncertainties are from HepData:
http://hepdata.cedar.ac.uk/view/ins1370682
1) TABS 15-16-17 HepData; TAB A.6 in the preprint
2) TABS 33-34-35 HepData; TAB A.9 in the preprint
3) TABS 21-22-23 HepData; TAB A.7 in the preprint
4) TABS 36-37-38 Hepdata; TAB A.9 in the preprint
5) TABS 39-40-41 HepData; TAB A.9 in the preprint

Notes:
1) The data is available separately for the lepton+jets and for the 
   dilepton channel.
2) The statistical covariance matrix is also available, so far we take 
   into account only the principal diagonal values.
3) All systematic uncertainties are assumed to be multiplicative.
4) Custom uncertainty descriptions are assumed to allow for cross-correlations
   among the five differential distributions. 
5) Unnormalized differential distributions are obtained by multiplying the
   normalized distributions by the total inclusive ttbar cross section   
   available in CMS PAS TOP-13-004. 
   Statistical uncertainties (on the total cross section and on the 
   normalised distributions) are added in quadrature.
   Systematic uncertainties (of the normalised distributions) factorise.
   Two addidional sources of systematics coming from the total cross 
   section (total sys and lumi) are considered on top of the full 
   breakdown of the systematics on normalised distributions.
*/
 
#include "CMSTOPDIFF.h"

//Define custom uncertainty descriptions to allow for cross-correlations
const std::vector<std::string> sysdescr = {
  /*
  "CMSTOPDIFFLEP",
  "CMSTOPDIFFJES",
  "CMSTOPDIFFJER",
  "CMSTOPDIFFBG",
  "CMSTOPDIFFBtag",
  "CMSTOPDIFFPU",
  "CMSTOPDIFFTopScale",
  "CMSTOPDIFFTopMatch",
  "CMSTOPDIFFHadronization",
  "CMSTOPDIFFTopMass",
  "CMSTOPDIFFPDF",
  "CMSTOPDIFFTotXSec",
  "CMSTOPDIFFLumi",
  */  
  "CORR",
  "CORR",
  "CORR",
  "CORR",
  "CORR",
  "CORR",
  "CORR",
  "CORR",
  "CORR",
  "CORR",
  "CORR",
  "CORR",
  "CMSLUMI12",
};

//A - NORMALISED distributions

//1) Distribution differential in top quark transverse momentum
void  CMSTOPDIFF8TEVTPTNORMFilter::ReadData()
{
  //Opening files
  fstream f1, f2, f3;

  //Central values and statistical uncertainties
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/CMSTOPDIFF8TEVTPT/CMSTOPDIFF8TEVTPT.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Statistical covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/CMSTOPDIFF8TEVTPT/CMSTOPDIFF8TEVTPT.cov";
  f3.open(covfile.str().c_str(), ios::in);

  if (f3.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Full breakdown of systematic uncertainties
  stringstream sysfile("");
  sysfile << dataPath()  
	  << "rawdata/CMSTOPDIFF8TEVTPT/CMSTOPDIFF8TEVTPT.sys";
  f2.open(sysfile.str().c_str(), ios::in);

  if (f2.fail()) 
    {
      cerr << "Error opening data file " << sysfile.str() << endl;
      exit(-1);
    }

  //Read central values and statistical uncertainty
  string line;
  for(int i=0; i<10; i++)
    {
      getline(f1,line);
    }

  int idum;
  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      lstream >> idum;
      double adum;
      for(int j=0; j<5; j++)
      	{
	  lstream >> adum;
	}
      double pt_top;
      lstream >> pt_top;
      
      fKin1[i] = pt_top;   //P_T^(top)
      fKin2[i] = Mt*Mt;       
      fKin3[i] = 8000;     //sqrt(s)
      
      lstream >> fData[i]; //normalized differential distribution
      lstream >> fStat[i]; 
      fStat[i]=0.;

      for(int j=0; j<3; j++)
	{
	  lstream >> adum;
	}
    }

  //Read statistical covariance matrix
  for(int i=0; i<8; i++)
    {
      getline(f3,line);
    }

  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
  {
    covmat[i] = new double[fNData];
    getline(f3,line);
    istringstream lstream(line);
    double adum;
    lstream >> adum >> adum >> adum;
    for(int j=0; j<fNData; j++)
    {
      lstream >> covmat[i][j];
    }
  }

  //Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
   {
     cerr << " in " << fSetName << endl;
     exit(-1);
   }

  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  fSys[i][j].add  = syscor[i][j];
	  fSys[i][j].mult = fSys[i][j].add*100/fData[i];
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = "CORR";
	}
    }

  //Read full breakdown of systematics
  for(int i=0; i<8; i++)
    {
      getline(f2,line);
    }

  for(int j=fNData; j<fNSys; j++)
    {
      string sdum;
      getline(f2,line);
      istringstream lstream(line);
      lstream >> sdum;
      for(int i=0; i<fNData; i++)
	{
	  lstream >> fSys[i][j].mult;
	  fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = sysdescr[j-fNData];
	}
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
  f3.close();

}

//==================================================================

//2) Distribution differential in top quark pair transverse momentum
void  CMSTOPDIFF8TEVTTPTNORMFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;

  //Central values and statistical uncertainties
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/CMSTOPDIFF8TEVTTPT/CMSTOPDIFF8TEVTTPT.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Statistical covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/CMSTOPDIFF8TEVTTPT/CMSTOPDIFF8TEVTTPT.cov";
  f3.open(covfile.str().c_str(), ios::in);

  if (f3.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Full breakdown of systematic uncertainties
  stringstream sysfile("");
  sysfile << dataPath() 
	  <<  "rawdata/CMSTOPDIFF8TEVTTPT/CMSTOPDIFF8TEVTTPT.sys";
  f2.open(sysfile.str().c_str(), ios::in);

  if (f2.fail()) 
    {
      cerr << "Error opening data file " << sysfile.str() << endl;
      exit(-1);
    }

  //Read central values and statistical uncertainties
  string line;
  for(int i=0; i<10; i++)
    {
      getline(f1,line);
    }

  int idum;
  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      lstream >> idum;
      double adum;
      for(int j=0; j<5; j++)
	{
	  lstream >> adum;
	}
      double pt_top;
      lstream >> pt_top;
      
      fKin1[i] = pt_top;   //P_T^(top)
      fKin2[i] = Mt*Mt;    
      fKin3[i] = 8000;     //sqrt(s)

      lstream >> fData[i]; //normalized differential distribution
      lstream >> fStat[i]; //assume stat errors uncorrelated so far
      fStat[i] =0.;

      for(int j=0; j<3; j++)
	{
	  lstream >> adum;
	}
    }

  //Read statistical covariance matrix
  for(int i=0; i<8; i++)
    {
      getline(f3,line);
    }

  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
  {
    covmat[i] = new double[fNData];
    getline(f3,line);
    istringstream lstream(line);
    double adum;
    lstream >> adum >> adum >> adum;
    for(int j=0; j<fNData; j++)
    {
      lstream >> covmat[i][j];
    }
  }

  //Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
   {
     cerr << " in " << fSetName << endl;
     exit(-1);
   }

  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  fSys[i][j].add  = syscor[i][j];
	  fSys[i][j].mult = fSys[i][j].add*100/fData[i];
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = "CORR";
	}
    }

  //Read full breakdown of systematics
  for(int i=0; i<8; i++)
    {
      getline(f2,line);
    }

  for(int j=fNData; j<fNSys; j++)
    {
      string sdum;
      getline(f2,line);
      istringstream lstream(line);
      lstream >> sdum;
      for(int i=0; i<fNData; i++)
	{
	  lstream >> fSys[i][j].mult;
	  fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = sysdescr[j-fNData];
	}
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
  f3.close();

}

//==============================================================

//3) Distribution differential in top quark rapidity
void  CMSTOPDIFF8TEVTRAPNORMFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;

  //Central values and statistical uncertainties
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/CMSTOPDIFF8TEVTRAP/CMSTOPDIFF8TEVTRAP.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Statistical covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/CMSTOPDIFF8TEVTRAP/CMSTOPDIFF8TEVTRAP.cov";
  f3.open(covfile.str().c_str(), ios::in);

  if (f3.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Full breakdown of systematic uncertainties
  stringstream sysfile("");
  sysfile << dataPath() 
	  << "rawdata/CMSTOPDIFF8TEVTRAP/CMSTOPDIFF8TEVTRAP.sys";
  f2.open(sysfile.str().c_str(), ios::in);

  if (f2.fail()) 
    {
      cerr << "Error opening data file " << sysfile.str() << endl;
      exit(-1);
    }

  //Read central values and statistical uncertainties
  string line;
  for(int i=0; i<10; i++)
    {
      getline(f1,line);
    }

  int idum;
  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      lstream >> idum;
      double adum;
      for(int j=0; j<5; j++)
	{
	  lstream >> adum;
	}
      double pt_top;
      lstream >> pt_top;
      
      fKin1[i] = pt_top;   //P_T^(top)
      fKin2[i] = Mt*Mt;    
      fKin3[i] = 8000;     //sqrt(s)

      lstream >> fData[i]; //normalized differential distribution
      lstream >> fStat[i]; //assume stat errors uncorrelated so far
      fStat[i] = 0.;

      for(int j=0; j<3; j++)
	{
	  lstream >> adum;
	}
    }

  //Read statistical covariance matrix
  for(int i=0; i<8; i++)
    {
      getline(f3,line);
    }

  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
  {
    covmat[i] = new double[fNData];
    getline(f3,line);
    istringstream lstream(line);
    double adum;
    lstream >> adum >> adum >> adum;
    for(int j=0; j<fNData; j++)
    {
      lstream >> covmat[i][j];
    }
  }

  //Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
   {
     cerr << " in " << fSetName << endl;
     exit(-1);
   }

  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  fSys[i][j].add  = syscor[i][j];
	  fSys[i][j].mult = fSys[i][j].add*100/fData[i];
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = "CORR";
	}
    }

  //Read full breakdown of systematic uncertainties
  for(int i=0; i<8; i++)
    {
      getline(f2,line);
    }

  for(int j=fNData; j<fNSys; j++)
    {
      string sdum;
      getline(f2,line);
      istringstream lstream(line);
      lstream >> sdum;
      
      for(int i=0; i<fNData; i++)
	{
	  lstream >> fSys[i][j].mult;	  
	  fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = sysdescr[j-fNData];
	  
	}
      
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
  f3.close();

}

//=================================================================

//4) Distribution differential in top quark pair rapidity
void  CMSTOPDIFF8TEVTTRAPNORMFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;

  //Central values and statistical uncertainties
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/CMSTOPDIFF8TEVTTRAP/CMSTOPDIFF8TEVTTRAP.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Statistical covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/CMSTOPDIFF8TEVTTRAP/CMSTOPDIFF8TEVTTRAP.cov";
  f3.open(covfile.str().c_str(), ios::in);

  if (f3.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Full breakdown of systematic uncertainties
  stringstream sysfile("");
  sysfile << dataPath() 
	  << "rawdata/CMSTOPDIFF8TEVTTRAP/CMSTOPDIFF8TEVTTRAP.sys";
  f2.open(sysfile.str().c_str(), ios::in);

  if (f2.fail()) 
    {
      cerr << "Error opening data file " << sysfile.str() << endl;
      exit(-1);
    }

  //Read central values and statistical uncertainties
  string line;
  for(int i=0; i<10; i++)
    {
      getline(f1,line);
    }

  int idum;
  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      lstream >> idum;
      double adum;
      for(int j=0; j<5; j++)
	{
	  lstream >> adum;
	}
      double pt_top;
      lstream >> pt_top;
      
      fKin1[i] = pt_top;   //P_T^(top)
      fKin2[i] = Mt*Mt;   
      fKin3[i] = 8000;     //sqrt(s)

      lstream >> fData[i]; //normalized differential distribution
      lstream >> fStat[i]; //assume stat errors uncorrelated so far
      fStat[i] = 0.;

      for(int j=0; j<3; j++)
	{
	  lstream >> adum;
	}
    }

  //Read statistical covariance matrix
  for(int i=0; i<8; i++)
    {
      getline(f3,line);
    }

  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
  {
    covmat[i] = new double[fNData];
    getline(f3,line);
    istringstream lstream(line);
    double adum;
    lstream >> adum >> adum >> adum;
    for(int j=0; j<fNData; j++)
    {
      lstream >> covmat[i][j];
    }
  }

  //Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
   {
     cerr << " in " << fSetName << endl;
     exit(-1);
   }

  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  fSys[i][j].add  = syscor[i][j];
	  fSys[i][j].mult = fSys[i][j].add*100/fData[i];
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = "CORR";
	}
    }

  //Read full breakdown of systematic uncertainties
  for(int i=0; i<8; i++)
    {
      getline(f2,line);
    }

  for(int j=fNData; j<fNSys; j++)
    {
      string sdum;
      getline(f2,line);
      istringstream lstream(line);
      lstream >> sdum;
      for(int i=0; i<fNData; i++)
	{
	  lstream >> fSys[i][j].mult;
	  fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = sysdescr[j-fNData];
	}
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
  f3.close();

}

//=================================================================

//5) Distribution differential in top quark pair invariant mass
void  CMSTOPDIFF8TEVTTMNORMFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;

  //Central values and statistical uncertainties
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/CMSTOPDIFF8TEVTTM/CMSTOPDIFF8TEVTTM.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Statistical covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/CMSTOPDIFF8TEVTTM/CMSTOPDIFF8TEVTTM.cov";
  f3.open(covfile.str().c_str(), ios::in);

  if (f3.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Full breakdown of systematic uncertainties
  stringstream sysfile("");
  sysfile << dataPath() 
	  << "rawdata/CMSTOPDIFF8TEVTTM/CMSTOPDIFF8TEVTTM.sys";
  f2.open(sysfile.str().c_str(), ios::in);

  if (f2.fail()) 
    {
      cerr << "Error opening data file " << sysfile.str() << endl;
      exit(-1);
    }

  //Read central values and statistical uncertainties
  string line;
  for(int i=0; i<10; i++)
    {
      getline(f1,line);
    }

  int idum;
  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      lstream >> idum;
      double adum;
      for(int j=0; j<5; j++)
	{
	  lstream >> adum;
	}
      double pt_top;
      lstream >> pt_top;
      
      fKin1[i] = pt_top;   //P_T^(top)
      fKin2[i] = Mt*Mt;    
      fKin3[i] = 8000;     //sqrt(s)

      lstream >> fData[i]; //normalized differential distribution
      lstream >> fStat[i]; //assume stat errors uncorrelated so far
      fStat[i] = 0.;

      for(int j=0; j<3; j++)
	{
	  lstream >> adum;
	}
    }

  //Read statistical covariance matrix
  for(int i=0; i<8; i++)
    {
      getline(f3,line);
    }

  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
  {
    covmat[i] = new double[fNData];
    getline(f3,line);
    istringstream lstream(line);
    double adum;
    lstream >> adum >> adum >> adum;
    for(int j=0; j<fNData; j++)
    {
      lstream >> covmat[i][j];
    }
  }

 //Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
   {
     cerr << " in " << fSetName << endl;
     exit(-1);
   }

  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  fSys[i][j].add  = syscor[i][j];
	  fSys[i][j].mult = fSys[i][j].add*100/fData[i];
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = "CORR";
	}
    }

  //Read full breakdown of systematic uncertainties
  for(int i=0; i<8; i++)
    {
      getline(f2,line);
    }

  for(int j=fNData; j<fNSys; j++)
    {
      string sdum;
      getline(f2,line);
      istringstream lstream(line);
      lstream >> sdum;
      for(int i=0; i<fNData; i++)
	{
	  lstream >> fSys[i][j].mult;
	  fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = sysdescr[j-fNData];
	}
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
  f3.close();

}

/*
========================================================================
*/

//B - UNNORMALISED distributions

//1) Distribution differential in top quark transverse momentum
void CMSTOPDIFF8TEVTPTFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3, f4;

  //Central values and statistical uncertainties
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/CMSTOPDIFF8TEVTPT/CMSTOPDIFF8TEVTPT.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Statistical covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/CMSTOPDIFF8TEVTPT/CMSTOPDIFF8TEVTPT.cov";
  f4.open(covfile.str().c_str(), ios::in);
  
  if (f4.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }


  //Full breakdown of systematic uncertainties
  stringstream sysfile("");
  sysfile << dataPath() 
	  << "rawdata/CMSTOPDIFF8TEVTPT/CMSTOPDIFF8TEVTPT.sys";
  f2.open(sysfile.str().c_str(), ios::in);

  if (f2.fail()) 
    {
      cerr << "Error opening data file " << sysfile.str() << endl;
      exit(-1);
    }

  //Inclusive total cross section
  stringstream sigmatotfile("");
  sigmatotfile << dataPath() 
	       << "rawdata/CMSTOPDIFF8TEVTOT/CMSTOPDIFF8TEVTOT.data";
  f3.open(sigmatotfile.str().c_str(), ios::in);

  if (f3.fail()) 
    {
      cerr << "Error opening data file " << sysfile.str() << endl;
      exit(-1);
    }

  //Starting filter
  string line;

  for(int i=0; i<3; i++)
    {
      getline(f3,line);
    }

  double xscv, xsstat, xssystpl, xssystmi, xslumi;
  double dtmp, stmp;

  f3 >> xscv >> xsstat >> xssystpl >> xssystmi >> xslumi;

  for(int i=0; i<10; i++)
    {
      getline(f1,line);
    }

  int idum;

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      lstream >> idum;
      double adum;
      for(int j=0; j<5; j++)
	{
      	  lstream >> adum;
      	}
      double pt_top;
      lstream >> pt_top;
      
      fKin1[i] = pt_top;   //P_T^(top)
      fKin2[i] = Mt*Mt;    
      fKin3[i] = 8000;     //sqrt(s)

      lstream >> fData[i]; //normalized differential distribution
      lstream >> fStat[i]; //assume stat errors uncorrelated so far
      fStat[i] = 0.;

      //differential distributions are unnormalized
      fData[i] = fData[i]*xscv;
      //statistical uncertainties are added in quadrature
      fStat[i] = fData[i]*pow(pow(fStat[i]/(fData[i]/(xscv)),2)+pow(xsstat/(xscv),2),0.5); 

      for(int j=0; j<3; j++)
	{
	  lstream >> adum;
	}
    }

  //Read statistical covariance matrix
  for(int i=0; i<8; i++)
    {
      getline(f4,line);
    }

  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
  {
    covmat[i] = new double[fNData];
    getline(f4,line);
    istringstream lstream(line);
    double adum;
    lstream >> adum >> adum >> adum;
    for(int j=0; j<fNData; j++)
    {
      lstream >> covmat[i][j];
    }
  }

  //Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];
  
  if(!genArtSys(fNData,covmat,syscor))
   {
     cerr << " in " << fSetName << endl;
     exit(-1);
   }

  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  fSys[i][j].add  = syscor[i][j]*xscv;
	  fSys[i][j].mult = fSys[i][j].add*100/fData[i];
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = "CORR";
	}
    }
  
  //Read full breakdown of systematic uncertainties
  for(int i=0; i<8; i++)
    {
      getline(f2,line);
    }

  for(int j=fNData; j<fNSys-2; j++)
    {
      string sdum;
      getline(f2,line);
      istringstream lstream(line);
      lstream >> sdum;
      for(int i=0; i<fNData; i++)
	{
	  //Relative systematic uncertainties are the same for both normalised and unnormalised distributions
	  lstream >> fSys[i][j].mult; 
	  fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = sysdescr[j-fNData];
	}
    }
 
  //Two additional sources of systematic uncertainties affect the unnormalised distributions
  for(int i=0; i<fNData; i++)
    {
      xssystpl = xssystpl/xscv*100;
      xssystmi = xssystmi/xscv*100;

      symmetriseErrors(xssystpl, xssystmi, &stmp, &dtmp); //symmetrise systematics
      fSys[i][19].add  = stmp*fData[i]/100;
      fSys[i][19].mult = stmp; 
      fSys[i][19].type = MULT;
      fSys[i][19].name = sysdescr[19-fNData];

      fSys[i][20].add  = fData[i]/xscv*xslumi;
      fSys[i][20].mult = fSys[i][20].add/fData[i]*100;
      fSys[i][20].type = MULT;
      fSys[i][20].name = sysdescr[20-fNData];

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
  f3.close();

}

//========================================================================

//2) Distribution differential in top quark pair transverse momentum
void  CMSTOPDIFF8TEVTTPTFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3, f4;

  //Central values and statistical uncertainties
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/CMSTOPDIFF8TEVTTPT/CMSTOPDIFF8TEVTTPT.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Statistical covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/CMSTOPDIFF8TEVTTPT/CMSTOPDIFF8TEVTTPT.cov";
  f4.open(covfile.str().c_str(), ios::in);

  if (f4.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Full breakdown of systematic uncertainties
  stringstream sysfile("");
  sysfile << dataPath() 
	  << "rawdata/CMSTOPDIFF8TEVTTPT/CMSTOPDIFF8TEVTTPT.sys";
  f2.open(sysfile.str().c_str(), ios::in);

  if (f2.fail()) 
    {
      cerr << "Error opening data file " << sysfile.str() << endl;
      exit(-1);
    }

  //Inclusive total cross section
  stringstream sigmatotfile("");
  sigmatotfile << dataPath() 
	       << "rawdata/CMSTOPDIFF8TEVTOT/CMSTOPDIFF8TEVTOT.data";
  f3.open(sigmatotfile.str().c_str(), ios::in);

  if (f3.fail()) 
    {
      cerr << "Error opening data file " << sysfile.str() << endl;
      exit(-1);
    }

  //Read central values and statistical uncertainties
  string line;

  for(int i=0; i<3; i++)
    {
      getline(f3,line);
    }

  double xscv, xsstat, xssystpl, xssystmi, xslumi;
  double dtmp, stmp;

  f3 >> xscv >> xsstat >> xssystpl >> xssystmi >> xslumi;

  for(int i=0; i<10; i++)
    {
      getline(f1,line);
    }

  int idum;

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      lstream >> idum;
      double adum;
      for(int j=0; j<5; j++)
	{
	  lstream >> adum;
	}
      double pt_top;
      lstream >> pt_top;
      
      fKin1[i] = pt_top;   //P_T^(top)
      fKin2[i] = Mt*Mt;    
      fKin3[i] = 8000;     //sqrt(s)

      lstream >> fData[i]; //normalized differential distribution
      lstream >> fStat[i]; //assume stat errors uncorrelated so far
      fStat[i] = 0.;

      //differential distributions are unnormalized
      fData[i] = fData[i]*xscv;
      //statistical uncertainties are added in quadrature
      fStat[i] = fData[i]*pow(pow(fStat[i]/(fData[i]/xscv),2)+pow(xsstat/xscv,2),0.5); 

      for(int j=0; j<3; j++)
	{
	  lstream >> adum;
	}
    }

  //Read statistical covariance matrix
  for(int i=0; i<8; i++)
    {
      getline(f4,line);
    }

  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
  {
    covmat[i] = new double[fNData];
    getline(f4,line);
    istringstream lstream(line);
    double adum;
    lstream >> adum >> adum >> adum;
    for(int j=0; j<fNData; j++)
    {
      lstream >> covmat[i][j];
    }
  }

  //Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
   {
     cerr << " in " << fSetName << endl;
     exit(-1);
   }

  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  fSys[i][j].add  = syscor[i][j]*xscv;
	  fSys[i][j].mult = fSys[i][j].add*100/fData[i];
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = "CORR";
	}
    }

  //Read full breakdown of systematic uncertainties
  for(int i=0; i<8; i++)
    {
      getline(f2,line);
    }

  for(int j=fNData; j<fNSys-2; j++)
    {
      string sdum;
      getline(f2,line);
      istringstream lstream(line);
      lstream >> sdum;
      for(int i=0; i<fNData; i++)
	{
	  //Relative systematic uncertainties are the same for both normalised and unnormalised distributions
	  lstream >> fSys[i][j].mult; 
	  fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = sysdescr[j-fNData];
	}
    }
 
  //Two additional sources of systematic uncertainties affect the unnormalised distributions
  for(int i=0; i<fNData; i++)
    { 
      xssystpl = xssystpl/xscv*100;
      xssystmi = xssystmi/xscv*100;

      symmetriseErrors(xssystpl, xssystmi, &stmp, &dtmp); //symmetrise systematics
      fSys[i][17].add  = fData[i]/100*stmp;
      fSys[i][17].mult = stmp; 
      fSys[i][17].type = MULT;
      fSys[i][17].name = sysdescr[17-fNData];

      fSys[i][18].add  = fData[i]/xscv*xslumi;
      fSys[i][18].mult = fSys[i][18].add/fData[i]*100;
      fSys[i][18].type = MULT;
      fSys[i][18].name = sysdescr[18-fNData];
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
  f3.close();
  f4.close();

}

//===================================================================

//3) Distribution differential in top quark rapidity
void  CMSTOPDIFF8TEVTRAPFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3, f4;

  //Central values and statistical uncertainties
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/CMSTOPDIFF8TEVTRAP/CMSTOPDIFF8TEVTRAP.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Statistical covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/CMSTOPDIFF8TEVTRAP/CMSTOPDIFF8TEVTRAP.cov";
  f4.open(covfile.str().c_str(), ios::in);

  if (f4.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Full breakdown of systematic uncertainties
  stringstream sysfile("");
  sysfile << dataPath() 
	  << "rawdata/CMSTOPDIFF8TEVTRAP/CMSTOPDIFF8TEVTRAP.sys";
  f2.open(sysfile.str().c_str(), ios::in);

  if (f2.fail()) 
    {
      cerr << "Error opening data file " << sysfile.str() << endl;
      exit(-1);
    }

  //Inclusive total cross section
  stringstream sigmatotfile("");
  sigmatotfile << dataPath() 
	       << "rawdata/CMSTOPDIFF8TEVTOT/CMSTOPDIFF8TEVTOT.data";
  f3.open(sigmatotfile.str().c_str(), ios::in);

  if (f3.fail()) 
    {
      cerr << "Error opening data file " << sysfile.str() << endl;
      exit(-1);
    }

  //Read central values and statistical uncertainties
  string line;

  for(int i=0; i<3; i++)
    {
      getline(f3,line);
    }

  double xscv, xsstat, xssystpl, xssystmi, xslumi;
  double dtmp, stmp;

  f3 >> xscv >> xsstat >> xssystpl >> xssystmi >> xslumi;

  for(int i=0; i<10; i++)
    {
      getline(f1,line);
    }

  int idum;

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      lstream >> idum;
      double adum;
      for(int j=0; j<5; j++)
	{
	  lstream >> adum;
	}
      double pt_top;
      lstream >> pt_top;
      
      fKin1[i] = pt_top;   //P_T^(top)
      fKin2[i] = Mt*Mt;    
      fKin3[i] = 8000;     //sqrt(s)

      lstream >> fData[i]; //normalized differential distribution
      lstream >> fStat[i]; //assume stat errors uncorrelated so far
      fStat[i] =0.;

      //differential distributions are unnormalized
      fData[i] = fData[i]*xscv;
      //statistical uncertainties are added in quadrature
      fStat[i] = fData[i]*pow(pow(fStat[i]/(fData[i]/xscv),2)+pow(xsstat/xscv,2),0.5); 

      for(int j=0; j<3; j++)
	{
	  lstream >> adum;
	}
    }

  //Read statistical covariance matrix
  for(int i=0; i<8; i++)
    {
      getline(f4,line);
    }

  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
  {
    covmat[i] = new double[fNData];
    getline(f4,line);
    istringstream lstream(line);
    double adum;
    lstream >> adum >> adum >> adum;
    for(int j=0; j<fNData; j++)
    {
      lstream >> covmat[i][j];
    }
  }

  //Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
   {
     cerr << " in " << fSetName << endl;
     exit(-1);
   }

  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  fSys[i][j].add  = syscor[i][j]*xscv;
	  fSys[i][j].mult = fSys[i][j].add*100/fData[i];
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = "CORR";
	}
    }

  //Read full breakdown of systematic uncertainties
  for(int i=0; i<8; i++)
    {
      getline(f2,line);
    }

  for(int j=fNData; j<fNSys-2; j++)
    {
      string sdum;
      getline(f2,line);
      istringstream lstream(line);
      lstream >> sdum;
      for(int i=0; i<fNData; i++)
	{
	  //Relative systematic uncertainties are the same for both normalised and unnormalised distributions
	  lstream >> fSys[i][j].mult; 
	  fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = sysdescr[j-fNData];
	}
    }

  //Two additional sources of systematic uncertainties affect the unnormalised distributions
  for(int i=0; i<fNData; i++)
    { 
      xssystpl = xssystpl/xscv*100;
      xssystmi = xssystmi/xscv*100;

      symmetriseErrors(xssystpl, xssystmi, &stmp, &dtmp); //symmetrise systematics
      fSys[i][21].add  = fData[i]/100*stmp;
      fSys[i][21].mult = stmp; 
      fSys[i][21].type = MULT;
      fSys[i][21].name = sysdescr[21-fNData];

      fSys[i][22].add  = fData[i]/xscv*xslumi;
      fSys[i][22].mult = fSys[i][22].add/fData[i]*100;
      fSys[i][22].type = MULT;
      fSys[i][22].name = sysdescr[22-fNData];
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
  f3.close();
  f4.close();

}

//==================================================================

//4) Distribution differential in top quark pairrapidity
void  CMSTOPDIFF8TEVTTRAPFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3, f4;

  //Central values and statistical uncertainties
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/CMSTOPDIFF8TEVTTRAP/CMSTOPDIFF8TEVTTRAP.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Statistical covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/CMSTOPDIFF8TEVTTRAP/CMSTOPDIFF8TEVTTRAP.cov";
  f4.open(covfile.str().c_str(), ios::in);

  if (f4.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Full breakdown of systematic uncertainties
  stringstream sysfile("");
  sysfile << dataPath() 
	  << "rawdata/CMSTOPDIFF8TEVTTRAP/CMSTOPDIFF8TEVTTRAP.sys";
  f2.open(sysfile.str().c_str(), ios::in);

  if (f2.fail()) 
    {
      cerr << "Error opening data file " << sysfile.str() << endl;
      exit(-1);
    }

  //Inclusive total cross section
  stringstream sigmatotfile("");
  sigmatotfile << dataPath() 
	       << "rawdata/CMSTOPDIFF8TEVTOT/CMSTOPDIFF8TEVTOT.data";
  f3.open(sigmatotfile.str().c_str(), ios::in);

  if (f3.fail()) 
    {
      cerr << "Error opening data file " << sysfile.str() << endl;
      exit(-1);
    }

  //Read central values and statistical uncertainties
  string line;

  for(int i=0; i<3; i++)
    {
      getline(f3,line);
    }

  double xscv, xsstat, xssystpl, xssystmi, xslumi;
  double dtmp, stmp;

  f3 >> xscv >> xsstat >> xssystpl >> xssystmi >> xslumi;

  for(int i=0; i<10; i++)
    {
      getline(f1,line);
    }

  int idum;

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      lstream >> idum;
      double adum;
      for(int j=0; j<5; j++)
	{
	  lstream >> adum;
	}
      double pt_top;
      lstream >> pt_top;
      
      fKin1[i] = pt_top;   //P_T^(top)
      fKin2[i] = Mt*Mt;    
      fKin3[i] = 8000;     //sqrt(s)

      lstream >> fData[i]; //normalized differential distribution
      lstream >> fStat[i]; //assume stat errors uncorrelated so far
      fStat[i] = 0.;

      //differential distributions are unnormalized
      fData[i] = fData[i]*xscv;
      //statistical uncertainties are added in quadrature
      fStat[i] = fData[i]*pow(pow(fStat[i]/(fData[i]/xscv),2)+pow(xsstat/xscv,2),0.5); 

      for(int j=0; j<3; j++)
	{
	  lstream >> adum;
	}
    }

 //Read statistical covariance matrix
  for(int i=0; i<8; i++)
    {
      getline(f4,line);
    }

  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
  {
    covmat[i] = new double[fNData];
    getline(f4,line);
    istringstream lstream(line);
    double adum;
    lstream >> adum >> adum >> adum;
    for(int j=0; j<fNData; j++)
    {
      lstream >> covmat[i][j];
    }
  }

  //Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
   {
     cerr << " in " << fSetName << endl;
     exit(-1);
   }

  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  fSys[i][j].add  = syscor[i][j]*xscv;
	  fSys[i][j].mult = fSys[i][j].add*100/fData[i];
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = "CORR";
	}
    }

  //Read full breakdown of systematic uncertainties
  for(int i=0; i<8; i++)
    {
      getline(f2,line);
    }

  for(int j=fNData; j<fNSys-2; j++)
    {
      string sdum;
      getline(f2,line);
      istringstream lstream(line);
      lstream >> sdum;
      for(int i=0; i<fNData; i++)
	{
	  //Relative systematic uncertainties are the same for both normalised and unnormalised distributions
	  lstream >> fSys[i][j].mult; 
	  fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = sysdescr[j-fNData];
	}
    }
 
  //Two additional sources of systematic uncertainties affect the unnormalised distributions
  for(int i=0; i<fNData; i++)
    {

      xssystpl = xssystpl/xscv*100;
      xssystmi = xssystmi/xscv*100;
      symmetriseErrors(xssystpl, xssystmi, &stmp, &dtmp); //symmetrise systematics

      fSys[i][21].add  = fData[i]/100*stmp;
      fSys[i][21].mult = stmp; 
      fSys[i][21].type = MULT;
      fSys[i][21].name = sysdescr[21-fNData];

      fSys[i][22].add  = fData[i]/xscv*xslumi;
      fSys[i][22].mult = fSys[i][22].add/fData[i]*100;
      fSys[i][22].type = MULT;
      fSys[i][22].name = sysdescr[22-fNData];
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
  f3.close();
  f4.close();

}

//====================================================================

//5) Distribution differential in top quark pair invariant mass
void  CMSTOPDIFF8TEVTTMFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3, f4;

  //Central values and statistical uncertainties
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/CMSTOPDIFF8TEVTTM/CMSTOPDIFF8TEVTTM.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Statistical covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/CMSTOPDIFF8TEVTTM/CMSTOPDIFF8TEVTTM.cov";
  f4.open(covfile.str().c_str(), ios::in);

  if (f4.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Full breakdown of systematic uncertainties
  stringstream sysfile("");
  sysfile << dataPath() 
	  << "rawdata/CMSTOPDIFF8TEVTTM/CMSTOPDIFF8TEVTTM.sys";
  f2.open(sysfile.str().c_str(), ios::in);

  if (f2.fail()) 
    {
      cerr << "Error opening data file " << sysfile.str() << endl;
      exit(-1);
    }

  //Inclusive total cross section
  stringstream sigmatotfile("");
  sigmatotfile << dataPath() 
	       << "rawdata/CMSTOPDIFF8TEVTOT/CMSTOPDIFF8TEVTOT.data";
  f3.open(sigmatotfile.str().c_str(), ios::in);

  if (f3.fail()) 
    {
      cerr << "Error opening data file " << sysfile.str() << endl;
      exit(-1);
    }

  //Read central values and statistical uncertainties
  string line;

  for(int i=0; i<3; i++)
    {
      getline(f3,line);
    }

  double xscv, xsstat, xssystpl, xssystmi, xslumi;
  double dtmp, stmp;

  f3 >> xscv >> xsstat >> xssystpl >> xssystmi >> xslumi;

  for(int i=0; i<10; i++)
    {
      getline(f1,line);
    }

  int idum;

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      lstream >> idum;
      double adum;
      for(int j=0; j<5; j++)
	{
	  lstream >> adum;
	}
      double pt_top;
      lstream >> pt_top;
      
      fKin1[i] = pt_top;   //P_T^(top)
      fKin2[i] = Mt*Mt;       
      fKin3[i] = 8000;     //sqrt(s)

      lstream >> fData[i]; //normalized differential distribution
      lstream >> fStat[i]; //assume stat errors uncorrelated so far
      fStat[i] = 0.;

      //differential distributions are unnormalized
      fData[i] = fData[i]*xscv;
      //statistical uncertainties are added in quadrature
      fStat[i] = fData[i]*pow(pow(fStat[i]/(fData[i]/xscv),2)+pow(xsstat/xscv,2),0.5); 

      for(int j=0; j<3; j++)
	{
	  lstream >> adum;
	}
    }

 //Read statistical covariance matrix
  for(int i=0; i<8; i++)
    {
      getline(f4,line);
    }

  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
  {
    covmat[i] = new double[fNData];
    getline(f4,line);
    istringstream lstream(line);
    double adum;
    lstream >> adum >> adum >> adum;
    for(int j=0; j<fNData; j++)
    {
      lstream >> covmat[i][j];
    }
  }

  //Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
   {
     cerr << " in " << fSetName << endl;
     exit(-1);
   }

  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  fSys[i][j].add  = syscor[i][j]*xscv;
	  fSys[i][j].mult = fSys[i][j].add*100/fData[i];
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = "CORR";
	}
    }

  //Read full breakdown of systematic uncertainties
  for(int i=0; i<8; i++)
    {
      getline(f2,line);
    }

  for(int j=fNData; j<fNSys-2; j++)
    {
      string sdum;
      getline(f2,line);
      istringstream lstream(line);
      lstream >> sdum;
      for(int i=0; i<fNData; i++)
	{
	  //Relative systematic uncertainties are the same for both normalised and unnormalised distributions
	  lstream >> fSys[i][j].mult; 
	  fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = sysdescr[j-fNData];
	}
    }
 
  //Two additional sources of systematic uncertainties affect the unnormalised distributions
  for(int i=0; i<fNData; i++)
    {
      
      xssystpl = xssystpl/xscv*100;
      xssystmi = xssystmi/xscv*100;
      symmetriseErrors(xssystpl, xssystmi, &stmp, &dtmp); //symmetrise systematics
      
      fSys[i][18].add  = fData[i]/100*stmp;
      fSys[i][18].mult = stmp; 
      fSys[i][18].type = MULT;
      fSys[i][18].name = sysdescr[18-fNData];

      fSys[i][19].add  = fData[i]/xscv*xslumi;
      fSys[i][19].mult = fSys[i][19].add/fData[i]*100;
      fSys[i][19].type = MULT;
      fSys[i][19].name = sysdescr[19-fNData];
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
  f3.close();
  f4.close();

}
