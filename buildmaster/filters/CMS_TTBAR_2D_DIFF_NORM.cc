/**
 * Experiment: CERN-LHC-CMS (CMS)
 * Preprinted as CMS-TOP-14-013
 * Archived as: ARXIV:1703.01630
 * Published in Eur.Phys.J. C77 (2017), 459
 *
 * Record in: INSPIRE
 * Record in: CERN Document Server
 * Record in: HEPData
 *
 * Normalized double-differential cross sections for top quark pair
 * (ttbar) production are measured in pp collisions at a
 * centre-of-mass energy of 8TeV with the CMS experiment at the
 * LHC. The analyzed data correspond to an integrated luminosity of
 * 19.7fb^-1. The measurement is performed in the dilepton e \pm
 * \mu^\mp final state. The ttbar cross section is determined as a
 * function of various pairs of observables characterizing the
 * kinematics of the top quark and ttbar system.
 *
 * Description of the buildmaster implementation
 * Normalized double differential cross sections for the distributions 
 * (lepton+jets channel) differential in the following variables are implemented:
 * 1) top quark transverse momentum and top quark rapidity;       
 * 2) top quark pair invariant mass and top quark rapidity;  
 * 3) top quark pair invariant mass and top quark pair rapidity;                   
 *
 * Raw data and full breakdown of systematic uncertainties are from HepData:
 * https://www.hepdata.net/record/ins1516191
 * 1) TABS 5-6-7 HepData; 
 * 2) TABS 8-9-10 HepData; 
 * 3) TABS 11-12-13 HepData; 
 *
 * Notes:
 * 1) Custom uncertainty descriptions are assumed to allow for cross-correlations
 *   among the three differential distributions. 
 * 2) Careful treatment is needed for the systematics in the .sys files. Where sys1
 *   and sys2 have opposite signs, + corresponds to a right (R) error and - to a 
 *   left (L) error.
 *   Where they have opposite signs, the largest
 *   value is the corresponding R or L error, and the other error is set to 0.
 *   Where there is only one value (in the case of Hadronization and Hard 
 *   scattering), the value (irrespective of sign) applies to both R and L 
 *   errors; L and R errors are treated as separate nuisance parameters and
 *   are normalised by sqrt(2) in order to be consistent with Eq.(6) in 
 *   arXiv:1703.01630.
 *
 **/
#include "CMS_TTBAR_2D_DIFF_NORM.h"

//Define custom uncertainty descriptions to allow for cross-correlations
const std::vector<std::string> sysdescr = { 
  "CORR",       // "CMSTTBAR2DDIFFJES+"
  "CORR",       // "CMSTTBAR2DDIFFJER+"
  "CORR",       // "CMSTTBAR2DDIFFKINR+"
  "CORR",       // "CMSTTBAR2DDIFFPU+"
  "CORR",       // "CMSTTBAR2DDIFFTRIG+" 
  "CORR",       // "CMSTTBAR2DDIFFBG1+"
  "CORR",       // "CMSTTBAR2DDIFFBG2+"
  "CORR",       // "CMSTTBAR2DDIFFBtag+"
  "CORR",       // "CMSTTBAR2DDIFFLumi+"
  "CORR",       // "CMSTTBAR2DDIFFTopMass+"
  "CORR",       // "CMSTTBAR2DDIFFTopScale+"
  "CORR",       // "CMSTTBAR2DDIFFTopMatch+"
  "CORR",       // "CMSTTBAR2DDIFFPDF+"
  "CORR",       // "CMSTTBAR2DDIFFJES-"
  "CORR",       // "CMSTTBAR2DDIFFJER-"
  "CORR",       // "CMSTTBAR2DDIFFKINR-"
  "CORR",       // "CMSTTBAR2DDIFFPU-"
  "CORR",       // "CMSTTBAR2DDIFFTRIG-" 
  "CORR",       // "CMSTTBAR2DDIFFBG1-"
  "CORR",       // "CMSTTBAR2DDIFFBG2-"
  "CORR",       // "CMSTTBAR2DDIFFBtag-"
  "CORR",       // "CMSTTBAR2DDIFFLumi-"
  "CORR",       // "CMSTTBAR2DDIFFTopMass-"
  "CORR",       // "CMSTTBAR2DDIFFTopScale-"
  "CORR",       // "CMSTTBAR2DDIFFTopMatch-"
  "CORR",       // "CMSTTBAR2DDIFFPDF-"
  "CORR",       // "CMSTTBAR2DDIFFHadronization"
  "CORR",       // "CMSTTBAR2DDIFFHardScat"
};

// Defining function to carry out routine for each distribution

//void process_data(string, string, string);

//1) Distribution differential 
//   in top quark transverse momentum and top quark rapidity
void  CMS_TTBAR_2D_DIFF_PT_TRAP_NORMFilter::ReadData()
{
  string data_filename = "rawdata/CMSTTBAR2DDIFF8TEVTPTTRAP/CMSTTBAR2DDIFF8TEVTPTTRAP.data";
  string cov_filename  = "rawdata/CMSTTBAR2DDIFF8TEVTPTTRAP/CMSTTBAR2DDIFF8TEVTPTTRAP.cov";
  string sys_filename  = "rawdata/CMSTTBAR2DDIFF8TEVTPTTRAP/CMSTTBAR2DDIFF8TEVTPTTRAP.sys";

  //Opening files
  fstream f1, f2, f3;
  
  //Central values and statistical uncertainties
  stringstream datafile("");
  datafile << dataPath() 
	   << data_filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Statistical correlation matrix
  stringstream covfile("");
  covfile << dataPath()
	  << cov_filename;
  f3.open(covfile.str().c_str(), ios::in);
  
  if (f3.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Full breakdown of systematic uncertainties
  stringstream sysfile("");
  sysfile << dataPath()  
	  << sys_filename;
  f2.open(sysfile.str().c_str(), ios::in);
  
  if (f2.fail()) 
    {
      cerr << "Error opening data file " << sysfile.str() << endl;
      exit(-1);
    }
  
  //Read central values and statistical uncertainty
  string line;
  // Skip comments
  for(int i=0; i<9; i++)
    {
      getline(f1,line);
      if (line[0] != '#' && line[0] != '\'')
	throw std::runtime_error("Error: line doesn't start with comment symbol");
    }
  
  vector<double> stat(fNData);
  
  for(int i=0; i<fNData; i++)
    {
      double var_1, var_2, ddum;
      char comma;
      
      getline(f1,line);
      istringstream lstream(line);
      lstream >> var_1  >> comma >> ddum >> comma >> ddum >> comma 
              >> var_2 >> comma >> ddum >> comma >> ddum >> comma
              >> fData[i] >> comma 
              >> stat[i] >> comma >> ddum >> comma
              >> ddum >> comma >> ddum;
      
      fKin1[i] = var_1;          
      fKin2[i] = pow(var_2,2);    
      fKin3[i] = 8000;           //sqrt(s) 
      fStat[i] = 0.;
    }
  
  //Read statistical correlation matrix
  // Skip comments
  for(int i=0; i<9; i++)
    {
      getline(f3,line);
      if (line[0] != '#' && line[0] != '\'')
	throw std::runtime_error("Error: line doesn't start with comment symbol");
    }

  //Create covmat of correct dimensions
  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
    {
      covmat[i] = new double[fNData];
      
      for (int j=i; j<fNData; j++)
	{
	  int row, col;
	  char comma;
	  
	  getline(f3,line);
	  istringstream lstream(line);
	  lstream >> row >> comma >> col >> comma >> covmat[i][j];
	  covmat[i][j] = stat[i]*fData[i]/100*stat[j]*fData[j]/100*covmat[i][j]/100;
	}
    }
  
  //Symmetrise covariance matrix
  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  if(i!=j)
	    covmat[j][i]=covmat[i][j];
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

  for(int i = 0; i < fNData; i++) 
    {
      delete[] syscor[i];
      delete[] covmat[i];
    }
  
  delete [] covmat; 
  delete [] syscor;
  
  //Read full breakdown of systematics
  // Skip comments
  for(int i=0; i<9; i++)
    {
      getline(f2,line);
      if (line[0] != '#' && line[0] != '\'')
	throw std::runtime_error("Error: line doesn't start with comment symbol");
    }

  //Create matrices to hold asymmmetric systematics  
  double** sys1 = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    sys1[i] = new double[15];
  double** sys2 = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    sys2[i] = new double[15];
  
  //Load sys1
  for(int idat=0; idat<fNData; idat++)
    {
      for(int isys=0; isys<15; isys++)
        {
          string bin;
          string sys_tag;
          string value;
	  
          getline(f2, line);
          stringstream ss(line);
          getline(ss,bin,',');
          getline(ss,sys_tag,',');
	  //Removing additional comma separated entry for TopScale
          if(isys == 10){
	    getline(ss,sys_tag,',');
          }
          getline(ss,value,',');       
          sys1[idat][isys] = atof(value.c_str());	  
	}
    }
  
  // Discard 4 intermediate lines
  // Skip comments
  for(int i=0; i<4; i++)
    {
      getline(f2,line);
      cout << line << endl;
      if (line[0] != '#' && line[0] != '\'' && line[0] != '\0')
	throw std::runtime_error("Error: line doesn't start with comment symbol");
    }

  //Load sys2
  for(int idat=0; idat<fNData; idat++)
    {
      for(int isys=0; isys<15; isys++)
        {
          string bin;
          string sys_tag;
          string value;

          getline(f2, line);
          stringstream ss(line);
          getline(ss,bin,',');
          getline(ss,sys_tag,',');
	  //Removing additional comma separated entry for TopScale
          if(isys == 10){
	    getline(ss,sys_tag,',');
          }
          getline(ss,value,',');     
	  sys2[idat][isys] = atof(value.c_str());
	}
    }
  
  //Sort out sys1 and sys2 if sys1 and sys2 have different signs arXiv:1703.01630
  for(int idat=0; idat<fNData; idat++)
    {
      for(int isys=0; isys<13; isys++) 
        {
	  double tmp1 = sys1[idat][isys];
	  double tmp2 = sys2[idat][isys];

	  //Case 1: sys1 and sys2 are both negative
	  if(tmp1<0.0 && tmp2<0.0)
	    {
	      if(tmp2<tmp1)
		{
		  sys1[idat][isys] = 0.0;
		  sys2[idat][isys] = tmp2;
		}
	      if(tmp2>tmp1)
		{
		  sys1[idat][isys] = 0.0;
		  sys2[idat][isys] = tmp1;
		}
	    }

	  //Case 2: sys1 and sys2 are both positive
	  if(tmp1>0.0 && tmp2>0.0)
	    {
	      if(tmp1>tmp2)
		{
		  sys1[idat][isys] = tmp1;
		  sys2[idat][isys] = 0.0;
		}
	      if(tmp1<tmp2)
		{
		  sys1[idat][isys] = tmp2;
		  sys2[idat][isys] = 0.0;
		}
	    }

	  //Case3: sys1 is negative and sys2 is positive
	  if(tmp1<0.0 && tmp2>0.0)
	    {
	      sys1[idat][isys] = tmp2;
	      sys2[idat][isys] = tmp1;
	    }

	  sys1[idat][isys] /= sqrt(2.);
	  sys2[idat][isys] /= sqrt(2.);
	}
    }
  
  //Array with full systematics
  const int  new_systematics = fNSys-fNData;
  double **sys = new double*[fNData];
  for(int i=0; i<fNData; i++)
    {
      sys[i] = new double[new_systematics];
      for(int isys=0; isys<(new_systematics-2)/2; isys++)
	sys[i][isys] = sys1[i][isys];
      
      for(int isys=(new_systematics-2)/2; isys<new_systematics-2; isys++)
	sys[i][isys] = sys2[i][isys-(new_systematics-2)/2];

      for(int isys=new_systematics-2; isys<new_systematics; isys++)
	sys[i][isys] = sys1[i][isys-(new_systematics-2)/2];
    }

  //Write systematics on file  
  for(int i=0; i<fNData; i++)
    {
      for(int j=fNData; j<fNSys; j++)
	{
	  fSys[i][j].mult = sys[i][j-fNData];
	  fSys[i][j].add  = fSys[i][j].mult*fData[i]/100;
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = sysdescr[j-fNData];
	}
      delete[] sys[i];
      delete[] sys1[i];
      delete[] sys2[i];
    }

  delete[] sys;
  delete[] sys1;
  delete[] sys2; 
  
  f1.close();
  f2.close();
  f3.close();

}

//2)Distribution differential in top quark pair invariant mass and top quark rapidity
void  CMS_TTBAR_2D_DIFF_MTT_TRAP_NORMFilter::ReadData()
{
  string data_filename = "rawdata/CMSTTBAR2DDIFF8TEVTTMTRAP/CMSTTBAR2DDIFF8TEVTPTTMTRAP.data";
  string cov_filename  = "rawdata/CMSTTBAR2DDIFF8TEVTTMTRAP/CMSTTBAR2DDIFF8TEVTPTTMTRAP.cov";
  string sys_filename  = "rawdata/CMSTTBAR2DDIFF8TEVTTMTRAP/CMSTTBAR2DDIFF8TEVTPTTMTRAP.sys";

  //Opening files
  fstream f1, f2, f3;
  
  //Central values and statistical uncertainties
  stringstream datafile("");
  datafile << dataPath() 
	   << data_filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Statistical correlation matrix
  stringstream covfile("");
  covfile << dataPath()
	  << cov_filename;
  f3.open(covfile.str().c_str(), ios::in);
  
  if (f3.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Full breakdown of systematic uncertainties
  stringstream sysfile("");
  sysfile << dataPath()  
	  << sys_filename;
  f2.open(sysfile.str().c_str(), ios::in);
  
  if (f2.fail()) 
    {
      cerr << "Error opening data file " << sysfile.str() << endl;
      exit(-1);
    }
  
  //Read central values and statistical uncertainty
  string line;
  // Skip comments
  for(int i=0; i<9; i++)
    {
      getline(f1,line);
      if (line[0] != '#' && line[0] != '\'')
	throw std::runtime_error("Error: line doesn't start with comment symbol");
    }
  
  vector<double> stat(fNData);
  
  for(int i=0; i<fNData; i++)
    {
      double var_1, var_2, ddum;
      char comma;
      
      getline(f1,line);
      istringstream lstream(line);
      lstream >> var_1  >> comma >> ddum >> comma >> ddum >> comma 
              >> var_2 >> comma >> ddum >> comma >> ddum >> comma
              >> fData[i] >> comma 
              >> stat[i] >> comma >> ddum >> comma
              >> ddum >> comma >> ddum;
      
      fKin1[i] = var_2;          
      fKin2[i] = pow(var_1,2);    
      fKin3[i] = 8000;           //sqrt(s) 
      fStat[i] = 0.;
    }
  
  //Read statistical correlation matrix
  // Skip comments
  for(int i=0; i<9; i++)
    {
      getline(f3,line);
      if (line[0] != '#' && line[0] != '\'')
	throw std::runtime_error("Error: line doesn't start with comment symbol");
    }

  //Create covmat of correct dimensions
  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
    {
      covmat[i] = new double[fNData];
      
      for (int j=i; j<fNData; j++)
	{
	  int row, col;
	  char comma;
	  
	  getline(f3,line);
	  istringstream lstream(line);
	  lstream >> row >> comma >> col >> comma >> covmat[i][j];
	  covmat[i][j] = stat[i]*fData[i]/100*stat[j]*fData[j]/100*covmat[i][j]/100;
	}
    }
  
  //Symmetrise covariance matrix
  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  if(i!=j)
	    covmat[j][i]=covmat[i][j];
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
  
  for(int i=0; i<fNData; i++)
   {
     delete [] covmat[i]; 
     delete [] syscor[i];
   }

  delete [] covmat; 
  delete [] syscor;
  
  //Read full breakdown of systematics
  // Skip comments
  for(int i=0; i<9; i++)
    {
      getline(f2,line);
      if (line[0] != '#' && line[0] != '\'')
	throw std::runtime_error("Error: line doesn't start with comment symbol");
    }

  //Create matrices to hold asymmmetric systematics  
  double** sys1 = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    sys1[i] = new double[15];
  double** sys2 = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    sys2[i] = new double[15];
  
  //Load sys1
  for(int idat=0; idat<fNData; idat++)
    {
      for(int isys=0; isys<15; isys++)
        {
          string bin;
          string sys_tag;
          string value;
	  
          getline(f2, line);
          stringstream ss(line);
          getline(ss,bin,',');
          getline(ss,sys_tag,',');
	  //Removing additional comma separated entry for TopScale
          if(isys == 10){
	    getline(ss,sys_tag,',');
          }
          getline(ss,value,',');       
          sys1[idat][isys] = atof(value.c_str());	  
	}
    }
  
  // Discard 4 intermediate lines
  // Skip comments
  for(int i=0; i<4; i++)
    {
      getline(f2,line);
      cout << line << endl;
      if (line[0] != '#' && line[0] != '\'' && line[0] != '\0')
	throw std::runtime_error("Error: line doesn't start with comment symbol");
    }

  //Load sys2
  for(int idat=0; idat<fNData; idat++)
    {
      for(int isys=0; isys<15; isys++)
        {
          string bin;
          string sys_tag;
          string value;

          getline(f2, line);
          stringstream ss(line);
          getline(ss,bin,',');
          getline(ss,sys_tag,',');
	  //Removing additional comma separated entry for TopScale
          if(isys == 10){
	    getline(ss,sys_tag,',');
          }
          getline(ss,value,',');     
	  sys2[idat][isys] = atof(value.c_str());
	}
    }
  
  //Sort out sys1 and sys2 if sys1 and sys2 have different signs arXiv:1703.01630
  for(int idat=0; idat<fNData; idat++)
    {
      for(int isys=0; isys<13; isys++) 
        {
	  double tmp1 = sys1[idat][isys];
	  double tmp2 = sys2[idat][isys];

	  //Case 1: sys1 and sys2 are both negative
	  if(tmp1<0.0 && tmp2<0.0)
	    {
	      if(tmp2<tmp1)
		{
		  sys1[idat][isys] = 0.0;
		  sys2[idat][isys] = tmp2;
		}
	      if(tmp2>tmp1)
		{
		  sys1[idat][isys] = 0.0;
		  sys2[idat][isys] = tmp1;
		}
	    }

	  //Case 2: sys1 and sys2 are both positive
	  if(tmp1>0.0 && tmp2>0.0)
	    {
	      if(tmp1>tmp2)
		{
		  sys1[idat][isys] = tmp1;
		  sys2[idat][isys] = 0.0;
		}
	      if(tmp1<tmp2)
		{
		  sys1[idat][isys] = tmp2;
		  sys2[idat][isys] = 0.0;
		}
	    }

	  //Case3: sys1 is negative and sys2 is positive
	  if(tmp1<0.0 && tmp2>0.0)
	    {
	      sys1[idat][isys] = tmp2;
	      sys2[idat][isys] = tmp1;
	    }

	  sys1[idat][isys] /= sqrt(2.);
	  sys2[idat][isys] /= sqrt(2.);
	}
    }
  
  //Array with full systematics
  const int  new_systematics = fNSys-fNData;
  double **sys = new double*[fNData];
  for(int i=0; i<fNData; i++)
    {
      sys[i] = new double[new_systematics];
      for(int isys=0; isys<(new_systematics-2)/2; isys++)
	sys[i][isys] = sys1[i][isys];
      
      for(int isys=(new_systematics-2)/2; isys<new_systematics-2; isys++)
	sys[i][isys] = sys2[i][isys-(new_systematics-2)/2];

      for(int isys=new_systematics-2; isys<new_systematics; isys++)
	sys[i][isys] = sys1[i][isys-(new_systematics-2)/2];
    }

  //Write systematics on file  
  for(int i=0; i<fNData; i++)
    {
      for(int j=fNData; j<fNSys; j++)
	{
	  fSys[i][j].mult = sys[i][j-fNData];
	  fSys[i][j].add  = fSys[i][j].mult*fData[i]/100;
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = sysdescr[j-fNData];
	}
      delete[] sys[i];
      delete[] sys1[i];
      delete[] sys2[i];
    }

  delete[] sys;
  delete[] sys1;
  delete[] sys2; 
  
  f1.close();
  f2.close();
  f3.close();

}

//3) Distribution differential in top quark pair invariant mass and top quark pair rapidity
void  CMS_TTBAR_2D_DIFF_MTT_TTRAP_NORMFilter::ReadData()
{
  string data_filename = "rawdata/CMSTTBAR2DDIFF8TEVTTMTTRAP/CMSTTBAR2DDIFF8TEVTPTTMTTRAP.data";
  string cov_filename  = "rawdata/CMSTTBAR2DDIFF8TEVTTMTTRAP/CMSTTBAR2DDIFF8TEVTPTTMTTRAP.cov";
  string sys_filename  = "rawdata/CMSTTBAR2DDIFF8TEVTTMTTRAP/CMSTTBAR2DDIFF8TEVTPTTMTTRAP.sys";

  //Opening files
  fstream f1, f2, f3;
  
  //Central values and statistical uncertainties
  stringstream datafile("");
  datafile << dataPath() 
	   << data_filename;
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Statistical correlation matrix
  stringstream covfile("");
  covfile << dataPath()
	  << cov_filename;
  f3.open(covfile.str().c_str(), ios::in);
  
  if (f3.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Full breakdown of systematic uncertainties
  stringstream sysfile("");
  sysfile << dataPath()  
	  << sys_filename;
  f2.open(sysfile.str().c_str(), ios::in);
  
  if (f2.fail()) 
    {
      cerr << "Error opening data file " << sysfile.str() << endl;
      exit(-1);
    }
  
  //Read central values and statistical uncertainty
  string line;
  // Skip comments
  for(int i=0; i<9; i++)
    {
      getline(f1,line);
      if (line[0] != '#' && line[0] != '\'')
	throw std::runtime_error("Error: line doesn't start with comment symbol");
    }
  
  vector<double> stat(fNData);
  
  for(int i=0; i<fNData; i++)
    {
      double var_1, var_2, ddum;
      char comma;
      
      getline(f1,line);
      istringstream lstream(line);
      lstream >> var_1  >> comma >> ddum >> comma >> ddum >> comma 
              >> var_2 >> comma >> ddum >> comma >> ddum >> comma
              >> fData[i] >> comma 
              >> stat[i] >> comma >> ddum >> comma
              >> ddum >> comma >> ddum;
      
      fKin1[i] = var_2;          
      fKin2[i] = pow(var_1,2);    
      fKin3[i] = 8000;           //sqrt(s) 
      fStat[i] = 0.;
    }
  
  //Read statistical correlation matrix
  // Skip comments
  for(int i=0; i<9; i++)
    {
      getline(f3,line);
      if (line[0] != '#' && line[0] != '\'')
	throw std::runtime_error("Error: line doesn't start with comment symbol");
    }

  //Create covmat of correct dimensions
  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
    {
      covmat[i] = new double[fNData];
      
      for (int j=i; j<fNData; j++)
	{
	  int row, col;
	  char comma;
	  
	  getline(f3,line);
	  istringstream lstream(line);
	  lstream >> row >> comma >> col >> comma >> covmat[i][j];
	  covmat[i][j] = stat[i]*fData[i]/100*stat[j]*fData[j]/100*covmat[i][j]/100;
	}
    }
  
  //Symmetrise covariance matrix
  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  if(i!=j)
	    covmat[j][i]=covmat[i][j];
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
  
  for(int i=0; i<fNData; i++)
   {
     delete [] covmat[i]; 
     delete [] syscor[i];
   }

  delete [] covmat; 
  delete [] syscor;
  
  //Read full breakdown of systematics
  // Skip comments
  for(int i=0; i<9; i++)
    {
      getline(f2,line);
      if (line[0] != '#' && line[0] != '\'')
	throw std::runtime_error("Error: line doesn't start with comment symbol");
    }

  //Create matrices to hold asymmmetric systematics  
  double** sys1 = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    sys1[i] = new double[15];
  double** sys2 = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    sys2[i] = new double[15];
  
  //Load sys1
  for(int idat=0; idat<fNData; idat++)
    {
      for(int isys=0; isys<15; isys++)
        {
          string bin;
          string sys_tag;
          string value;
	  
          getline(f2, line);
          stringstream ss(line);
          getline(ss,bin,',');
          getline(ss,sys_tag,',');
	  //Removing additional comma separated entry for TopScale
          if(isys == 10){
	    getline(ss,sys_tag,',');
          }
          getline(ss,value,',');       
          sys1[idat][isys] = atof(value.c_str());	  
	}
    }
  
  // Discard 4 intermediate lines
  // Skip comments
  for(int i=0; i<4; i++)
    {
      getline(f2,line);
      cout << line << endl;
      if (line[0] != '#' && line[0] != '\'' && line[0] != '\0')
	throw std::runtime_error("Error: line doesn't start with comment symbol");
    }

  //Load sys2
  for(int idat=0; idat<fNData; idat++)
    {
      for(int isys=0; isys<15; isys++)
        {
          string bin;
          string sys_tag;
          string value;

          getline(f2, line);
          stringstream ss(line);
          getline(ss,bin,',');
          getline(ss,sys_tag,',');
	  //Removing additional comma separated entry for TopScale
          if(isys == 10){
	    getline(ss,sys_tag,',');
          }
          getline(ss,value,',');     
	  sys2[idat][isys] = atof(value.c_str());
	}
    }
  
  //Sort out sys1 and sys2 if sys1 and sys2 have different signs arXiv:1703.01630
  for(int idat=0; idat<fNData; idat++)
    {
      for(int isys=0; isys<13; isys++) 
        {
	  double tmp1 = sys1[idat][isys];
	  double tmp2 = sys2[idat][isys];

	  //Case 1: sys1 and sys2 are both negative
	  if(tmp1<0.0 && tmp2<0.0)
	    {
	      if(tmp2<tmp1)
		{
		  sys1[idat][isys] = 0.0;
		  sys2[idat][isys] = tmp2;
		}
	      if(tmp2>tmp1)
		{
		  sys1[idat][isys] = 0.0;
		  sys2[idat][isys] = tmp1;
		}
	    }

	  //Case 2: sys1 and sys2 are both positive
	  if(tmp1>0.0 && tmp2>0.0)
	    {
	      if(tmp1>tmp2)
		{
		  sys1[idat][isys] = tmp1;
		  sys2[idat][isys] = 0.0;
		}
	      if(tmp1<tmp2)
		{
		  sys1[idat][isys] = tmp2;
		  sys2[idat][isys] = 0.0;
		}
	    }

	  //Case3: sys1 is negative and sys2 is positive
	  if(tmp1<0.0 && tmp2>0.0)
	    {
	      sys1[idat][isys] = tmp2;
	      sys2[idat][isys] = tmp1;
	    }

	  sys1[idat][isys] /= sqrt(2.);
	  sys2[idat][isys] /= sqrt(2.);
	}
    }
  
  //Array with full systematics
  const int  new_systematics = fNSys-fNData;
  double **sys = new double*[fNData];
  for(int i=0; i<fNData; i++)
    {
      sys[i] = new double[new_systematics];
      for(int isys=0; isys<(new_systematics-2)/2; isys++)
	sys[i][isys] = sys1[i][isys];
      
      for(int isys=(new_systematics-2)/2; isys<new_systematics-2; isys++)
	sys[i][isys] = sys2[i][isys-(new_systematics-2)/2];

      for(int isys=new_systematics-2; isys<new_systematics; isys++)
	sys[i][isys] = sys1[i][isys-(new_systematics-2)/2];
    }

  //Write systematics on file  
  for(int i=0; i<fNData; i++)
    {
      for(int j=fNData; j<fNSys; j++)
	{
	  fSys[i][j].mult = sys[i][j-fNData];
	  fSys[i][j].add  = fSys[i][j].mult*fData[i]/100;
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = sysdescr[j-fNData];
	}
      delete[] sys[i];
      delete[] sys1[i];
      delete[] sys2[i];
    }

  delete[] sys;
  delete[] sys1;
  delete[] sys2; 
  
  f1.close();
  f2.close();
  f3.close();

}

