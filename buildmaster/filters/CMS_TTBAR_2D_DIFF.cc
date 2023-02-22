
/*
Experiment: CERN-LHC-CMS (CMS)
Preprinted as CMS-TOP-14-013
Archived as: ARXIV:1703.01630
Published in Eur.Phys.J. C77 (2017), 459

Record in: INSPIRE
Record in: CERN Document Server
Record in: HEPData 

Normalized double-differential cross sections for top quark pair (tt⎯⎯) 
production are measured in pp collisions at a centre-of-mass energy of 8TeV 
with the CMS experiment at the LHC. The analyzed data correspond to an 
integrated luminosity of 19.7fb−1. The measurement is performed in the 
dilepton e±μ∓ final state. The tt⎯⎯ cross section is determined as a function 
of various pairs of observables characterizing the kinematics of the top quark 
and tt⎯⎯ system. The data are compared to calculations using perturbative 
quantum chromodynamics at next-to-leading and approximate 
next-to-next-to-leading orders. They are also compared to predictions of 
Monte Carlo event generators that complement fixed-order computations with 
parton showers, hadronization, and multiple-parton interactions. 
Overall agreement is observed with the predictions, which is improved when 
the latest global sets of proton parton distribution functions are used. 
The inclusion of the measured tt⎯⎯ cross sections in a fit of parametrized 
parton distribution functions is shown to have significant impact on the gluon 
distribution.

Description of the buildmaster implementation
Normalized double differential cross sections for the distributions 
(lepton+jets channel) differential in the following variables are implemented:
1) top quark transverse momentum and top quark rapidity;       
2) top quark pair invariant mass and top quark rapidity;  
3) top quark pair invariant mass and top quark pair rapidity;                  
4) top quark pair invariant mass and top quark pair transverse momentum;   

Raw data and full breakdown of systematic uncertainties are from HepData:
https://www.hepdata.net/record/ins1516191
1) TABS 5-6-7 HepData; 
2) TABS 8-9-10 HepData; 
3) TABS 11-12-13 HepData; 
4) TABS 17-18-19 Hepdata;  


Notes:
1) The statistical covariance matrix is also available, so far we take 
   into account only the principal diagonal values.
2) All systematic uncertainties are assumed to be multiplicative.
3) Custom uncertainty descriptions are assumed to allow for cross-correlations
   among the five differential distributions. 
4) Unnormalized differential distributions are obtained by multiplying the
   normalized distributions by the total inclusive ttbar cross section   
   available in CMS PAS TOP-13-004. 
   Statistical uncertainties (on the total cross section and on the 
   normalised distributions) are added in quadrature.
   Systematic uncertainties (of the normalised distributions) factorise.
   Two addidional sources of systematics coming from the total cross 
   section (total sys and lumi) are considered on top of the full 
   breakdown of the systematics on normalised distributions.
5) Careful treatment is needed for the systematics in the .sys files. Where sys1
   and sys2 have opposite signs, + corresponds to a R error and - to a L error.
   Where they have opposite signs, the largest
   value is the corresponding R or L error, and the other error is set to 0.
   Where there is only one value (in the case of Hadronization and Hard 
   scattering), the value (irrespective of sign) applies to both R and L 
   errors; L and R errors are treated as separate nuisance parameters and
   are normalised by sqrt(2) in order to be consistent with Eq.(6) in 
   arXiv:1703.01630.
*/

 
#include "CMS_TTBAR_2D_DIFF.h"

//Define custom uncertainty descriptions to allow for cross-correlations
const std::vector<std::string> sysdescr = { 
  "CORR",       // "CMSTTBAR2DDIFFUNNORMJES+"
  "CORR",       // "CMSTTBAR2DDIFFUNNORMJER+"
  "CORR",       // "CMSTTBAR2DDIFFUNNORMKINR+"
  "CORR",       // "CMSTTBAR2DDIFFUNNORMPU+"
  "CORR",       // "CMSTTBAR2DDIFFUNNORMTRIG+" 
  "CORR",       // "CMSTTBAR2DDIFFUNNORMBG1+"
  "CORR",       // "CMSTTBAR2DDIFFUNNORMBG2+"
  "CORR",       // "CMSTTBAR2DDIFFUNNORMBtag+"
  "CORR",       // "CMSTTBAR2DDIFFUNNORMLumi+"
  "CORR",       // "CMSTTBAR2DDIFFUNNORMTopMass+"
  "CORR",       // "CMSTTBAR2DDIFFUNNORMTopScale+"
  "CORR",       // "CMSTTBAR2DDIFFUNNORMTopMatch+"
  "CORR",       // "CMSTTBAR2DDIFFUNNORMPDF+"
  "CORR",       // "CMSTTBAR2DDIFFUNNORMJES-"
  "CORR",       // "CMSTTBAR2DDIFFUNNORMJER-"
  "CORR",       // "CMSTTBAR2DDIFFUNNORMKINR-"
  "CORR",       // "CMSTTBAR2DDIFFUNNORMPU-"
  "CORR",       // "CMSTTBAR2DDIFFUNNORMTRIG-" 
  "CORR",       // "CMSTTBAR2DDIFFUNNORMBG1-"
  "CORR",       // "CMSTTBAR2DDIFFUNNORMBG2-"
  "CORR",       // "CMSTTBAR2DDIFFUNNORMBtag-"
  "CORR",       // "CMSTTBAR2DDIFFUNNORMLumi-"
  "CORR",       // "CMSTTBAR2DDIFFUNNORMTopMass-"
  "CORR",       // "CMSTTBAR2DDIFFUNNORMTopScale-"
  "CORR",       // "CMSTTBAR2DDIFFUNNORMTopMatch-"
  "CORR",       // "CMSTTBAR2DDIFFUNNORMPDF-"
  "CORR",       // "CMSTTBAR2DDIFFUNNORMHadronization"
  "CORR",       // "CMSTTBAR2DDIFFUNNORMHardScat"
};

//1) Distribution differential 
//   in top quark transverse momentum and top quark rapidity
void  CMS_TTBAR_2D_DIFF_PT_TRAPFilter::ReadData()
{
  //Fiducial cross section
  const double xsec      = 244.9*1000; //[fb]
  const double xsec_stat = 1.4*1000;   //[fb]
  const double xsec_syst = 5.9*1000;   //[fb]
  const double xsec_lumi = 6.4*1000;   //[fb]

  //Opening files
  fstream f1, f2, f3;
  
  //Central values and statistical uncertainties
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/CMSTTBAR2DDIFF8TEVTPTTRAP/CMSTTBAR2DDIFF8TEVTPTTRAP.data";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Statistical correlation matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/CMSTTBAR2DDIFF8TEVTPTTRAP/CMSTTBAR2DDIFF8TEVTPTTRAP.cov";
  f3.open(covfile.str().c_str(), ios::in);
  
  if (f3.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Full breakdown of systematic uncertainties
  stringstream sysfile("");
  sysfile << dataPath()  
	  << "rawdata/CMSTTBAR2DDIFF8TEVTPTTRAP/CMSTTBAR2DDIFF8TEVTPTTRAP.sys";
  f2.open(sysfile.str().c_str(), ios::in);
  
  if (f2.fail()) 
    {
      cerr << "Error opening data file " << sysfile.str() << endl;
      exit(-1);
    }
  
  //Read central values and statistical uncertainty
  string line;
  for(int i=0; i<9; i++)
    {
      getline(f1,line);
    }
  
  double stat[fNData];
  
  for(int i=0; i<fNData; i++)
    {
      double y_top, pt_top, ddum;
      char comma;
      
      getline(f1,line);
      istringstream lstream(line);
      lstream >> y_top  >> comma >> ddum >> comma >> ddum >> comma 
              >> pt_top >> comma >> ddum >> comma >> ddum >> comma
              >> fData[i] >> comma 
              >> stat[i] >> comma >> ddum >> comma
              >> ddum >> comma >> ddum;
      
      fKin1[i] = y_top;          //yt
      fKin2[i] = pow(pt_top,2);  //pTt  
      fKin3[i] = 8000;           //sqrt(s) 
      fStat[i] = 0.;
    }
  
  //Read statistical correlation matrix
  for(int i=0; i<9; i++)
    {
      getline(f3,line);
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
      fData[i] *= xsec;

      for(int j=0; j<fNData; j++)
	{
	  fSys[i][j].add  = syscor[i][j]*xsec;
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
  for(int i=0; i<9; i++)
    {
      getline(f2,line);
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
  
  //Discard 4 intermediate lines
  for(int i=0; i<4; i++)
    {
      getline(f2,line);
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

	  sys1[idat][isys] = sys1[idat][isys]/sqrt(2.);
	  sys2[idat][isys] = sys2[idat][isys]/sqrt(2.);
	}
    }
  
  //Array with full systematics
  double sys[fNData][fNSys-3-fNData];
  for(int i=0; i<fNData; i++)
    {
      for(int isys=0; isys<(fNSys-3-fNData-2)/2; isys++)
	sys[i][isys] = sys1[i][isys];
      
      for(int isys=(fNSys-3-fNData-2)/2; isys<fNSys-3-fNData-2; isys++)
	sys[i][isys] = sys2[i][isys-(fNSys-3-fNData-2)/2];

      for(int isys=fNSys-3-fNData-2; isys<fNSys-3-fNData; isys++)
	sys[i][isys] = sys1[i][isys-(fNSys-3-fNData-2)/2];
    }

  //Write systematics on file  
  for(int i=0; i<fNData; i++)
    {
      for(int j=fNData; j<fNSys-3; j++)
	{
	  fSys[i][j].mult = sys[i][j-fNData];
	  fSys[i][j].add  = fSys[i][j].mult*fData[i]/100;
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = sysdescr[j-fNData];
	}
      
      fSys[i][fNSys-3].mult = xsec_stat/xsec ;
      fSys[i][fNSys-3].add  = fSys[i][fNSys-3].mult*fData[i]/100;
      fSys[i][fNSys-3].type = MULT;
      fSys[i][fNSys-3].name = "UNCORR";

      fSys[i][fNSys-2].mult = xsec_syst/xsec ;
      fSys[i][fNSys-2].add  = fSys[i][fNSys-2].mult*fData[i]/100;
      fSys[i][fNSys-2].type = MULT;
      fSys[i][fNSys-2].name = "CORR";

      fSys[i][fNSys-1].mult = xsec_lumi/xsec ;
      fSys[i][fNSys-1].add  = fSys[i][fNSys-1].mult*fData[i]/100;
      fSys[i][fNSys-1].type = MULT;
      fSys[i][fNSys-1].name = "CMSLUMI12";

    }

  for(int i = 0; i < fNData; i++) 
    {
      delete[] sys1[i];
      delete[] sys2[i];
    }
  
  delete [] sys1;
  delete [] sys2; 
  
  f1.close();
  f2.close();
  f3.close();
  
}

//2)Distribution differential in top quark pair invariant mass and top quark rapidity
void  CMS_TTBAR_2D_DIFF_MTT_TRAPFilter::ReadData()
{
  //Fiducial cross section
  const double xsec      = 244.9*1000; //[pb]
  const double xsec_stat = 1.4*1000;   //[pb]
  const double xsec_syst = 5.9*1000;   //[pb]
  const double xsec_lumi = 6.4*1000;   //[pb]

  //Opening files
  fstream f1, f2, f3;
  
  //Central values and statistical uncertainties
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/CMSTTBAR2DDIFF8TEVTTMTRAP/CMSTTBAR2DDIFF8TEVTPTTMTRAP.data";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Statistical correlation matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/CMSTTBAR2DDIFF8TEVTTMTRAP/CMSTTBAR2DDIFF8TEVTPTTMTRAP.cov";
  f3.open(covfile.str().c_str(), ios::in);
  
  if (f3.fail()) 
    {
      cerr << "Error opening data file " << covfile.str() << endl;
      exit(-1);
    }
  
  //Full breakdown of systematic uncertainties
  stringstream sysfile("");
  sysfile << dataPath()  
	  << "rawdata/CMSTTBAR2DDIFF8TEVTTMTRAP/CMSTTBAR2DDIFF8TEVTPTTMTRAP.sys";
  f2.open(sysfile.str().c_str(), ios::in);
  
  if (f2.fail()) 
    {
      cerr << "Error opening data file " << sysfile.str() << endl;
      exit(-1);
    }
  
  //Read central values and statistical uncertainty
  string line;
  for(int i=0; i<9; i++)
    {
      getline(f1,line);
    }
  
  double stat[fNData];
  
  for(int i=0; i<fNData; i++)
    {
      double m_tt, y_top, ddum;
      char comma;
      
      getline(f1,line);
      istringstream lstream(line);
      lstream >> m_tt  >> comma >> ddum >> comma >> ddum >> comma 
              >> y_top >> comma >> ddum >> comma >> ddum >> comma
              >> fData[i] >> comma 
              >> stat[i] >> comma >> ddum >> comma
              >> ddum >> comma >> ddum;

      fKin1[i] = y_top;        //yt   
      fKin2[i] = pow(m_tt,2);  //mtt  
      fKin3[i] = 8000;         //sqrt(s) 
      fStat[i] = 0.;
    }

  //Read statistical correlation matrix
  for(int i=0; i<9; i++)
    {
      getline(f3,line);
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
	    covmat[j][i]= covmat[i][j];
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
      fData[i] *= xsec;

      for(int j=0; j<fNData; j++)
	{
	  fSys[i][j].add  = syscor[i][j]*xsec;
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
  for(int i=0; i<9; i++)
    {
      getline(f2,line);
      
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
  
  //Discard 4 intermediate lines
  for(int i=0; i<4; i++)
    {
      getline(f2,line);
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

	  sys1[idat][isys] = sys1[idat][isys]/sqrt(2.);
	  sys2[idat][isys] = sys2[idat][isys]/sqrt(2.);
	}
    }
  
  //Array with full systematics
  double sys[fNData][fNSys-3-fNData];
  for(int i=0; i<fNData; i++)
    {
      for(int isys=0; isys<(fNSys-3-fNData-2)/2; isys++)
	sys[i][isys] = sys1[i][isys];
      
      for(int isys=(fNSys-3-fNData-2)/2; isys<fNSys-3-fNData-2; isys++)
	sys[i][isys] = sys2[i][isys-(fNSys-3-fNData-2)/2];

      for(int isys=fNSys-3-fNData-2; isys<fNSys-3-fNData; isys++)
	sys[i][isys] = sys1[i][isys-(fNSys-3-fNData-2)/2];
    }

  //Write systematics on file  
  for(int i=0; i<fNData; i++)
    {
      for(int j=fNData; j<fNSys-3; j++)
	{
	  fSys[i][j].mult = sys[i][j-fNData];
	  fSys[i][j].add  = fSys[i][j].mult*fData[i]/100;
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = sysdescr[j-fNData];
	}

      fSys[i][fNSys-3].mult = xsec_stat/xsec ;
      fSys[i][fNSys-3].add  = fSys[i][fNSys-3].mult*fData[i]/100;
      fSys[i][fNSys-3].type = MULT;
      fSys[i][fNSys-3].name = "UNCORR";

      fSys[i][fNSys-2].mult = xsec_syst/xsec ;
      fSys[i][fNSys-2].add  = fSys[i][fNSys-2].mult*fData[i]/100;
      fSys[i][fNSys-2].type = MULT;
      fSys[i][fNSys-2].name = "CORR";

      fSys[i][fNSys-1].mult = xsec_lumi/xsec ;
      fSys[i][fNSys-1].add  = fSys[i][fNSys-1].mult*fData[i]/100;
      fSys[i][fNSys-1].type = MULT;
      fSys[i][fNSys-1].name = "CMSLUMI12";

    }

  for(int i = 0; i < fNData; i++) 
    {
      delete[] sys1[i];
      delete[] sys2[i];
    }
  
  delete [] sys1;
  delete [] sys2;

  f1.close();
  f2.close();
  f3.close();

}

//3) Distribution differential in top quark pair invariant mass and top quark pair rapidity
void  CMS_TTBAR_2D_DIFF_MTT_TTRAPFilter::ReadData()
{
  //Fiducial cross section
  const double xsec      = 244.9*1000; //[pb]
  const double xsec_stat = 1.4*1000;   //[pb]
  const double xsec_syst = 5.9*1000;   //[pb]
  const double xsec_lumi = 6.4*1000;   //[pb]

  //Opening files
  fstream f1, f2, f3;
  
  //Central values and statistical uncertainties
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/CMSTTBAR2DDIFF8TEVTTMTTRAP/CMSTTBAR2DDIFF8TEVTPTTMTTRAP.data";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Statistical correlation matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/CMSTTBAR2DDIFF8TEVTTMTTRAP/CMSTTBAR2DDIFF8TEVTPTTMTTRAP.cov";
  f3.open(covfile.str().c_str(), ios::in);
  
  if (f3.fail()) 
    {
      cerr << "Error opening data file " << covfile.str() << endl;
      exit(-1);
    }
  
  //Full breakdown of systematic uncertainties
  stringstream sysfile("");
  sysfile << dataPath()  
	  << "rawdata/CMSTTBAR2DDIFF8TEVTTMTTRAP/CMSTTBAR2DDIFF8TEVTPTTMTTRAP.sys";
  f2.open(sysfile.str().c_str(), ios::in);
  
  if (f2.fail()) 
    {
      cerr << "Error opening data file " << sysfile.str() << endl;
      exit(-1);
    }
  
  //Read central values and statistical uncertainty
  string line;
  for(int i=0; i<9; i++)
    {
      getline(f1,line);
    }
  
  double stat[fNData];
  
  for(int i=0; i<fNData; i++)
    {
      double m_tt, y_tt, ddum;
      char comma;
      
      getline(f1,line);
      istringstream lstream(line);
      lstream >> m_tt  >> comma >> ddum >> comma >> ddum >> comma 
              >> y_tt >> comma >> ddum >> comma >> ddum >> comma
              >> fData[i] >> comma 
              >> stat[i] >> comma >> ddum >> comma
              >> ddum >> comma >> ddum;

      fKin1[i] = y_tt;        //y_tt
      fKin2[i] = pow(m_tt,2); //m_tt     
      fKin3[i] = 8000;        //sqrt(s) 
      fStat[i] = 0.;
    }
  
  //Read statistical correlation matrix
  for(int i=0; i<9; i++)
    {
      getline(f3,line);
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
	    covmat[j][i]= covmat[i][j];
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
      fData[i] *= xsec;

      for(int j=0; j<fNData; j++)
	{
	  fSys[i][j].add  = syscor[i][j]*xsec;
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
  for(int i=0; i<9; i++)
    {
      getline(f2,line);
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
  
  //Discard 4 intermediate lines
  for(int i=0; i<4; i++)
    {
      getline(f2,line);
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

	  sys1[idat][isys] = sys1[idat][isys]/sqrt(2.);
	  sys2[idat][isys] = sys2[idat][isys]/sqrt(2.);
	}
    }
  
  //Array with full systematics
  double sys[fNData][fNSys-3-fNData];
  for(int i=0; i<fNData; i++)
    {
      for(int isys=0; isys<(fNSys-3-fNData-2)/2; isys++)
	sys[i][isys] = sys1[i][isys];
      
      for(int isys=(fNSys-3-fNData-2)/2; isys<fNSys-3-fNData-2; isys++)
	sys[i][isys] = sys2[i][isys-(fNSys-3-fNData-2)/2];

      for(int isys=fNSys-3-fNData-2; isys<fNSys-3-fNData; isys++)
	sys[i][isys] = sys1[i][isys-(fNSys-3-fNData-2)/2];
    }

  //Write systematics on file  
  for(int i=0; i<fNData; i++)
    {
      for(int j=fNData; j<fNSys-3; j++)
	{
	  fSys[i][j].mult = sys[i][j-fNData];
	  fSys[i][j].add  = fSys[i][j].mult*fData[i]/100;
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = sysdescr[j-fNData];
	}

      fSys[i][fNSys-3].mult = xsec_stat/xsec ;
      fSys[i][fNSys-3].add  = fSys[i][fNSys-3].mult*fData[i]/100;
      fSys[i][fNSys-3].type = MULT;
      fSys[i][fNSys-3].name = "UNCORR";

      fSys[i][fNSys-2].mult = xsec_syst/xsec ;
      fSys[i][fNSys-2].add  = fSys[i][fNSys-2].mult*fData[i]/100;
      fSys[i][fNSys-2].type = MULT;
      fSys[i][fNSys-2].name = "CORR";

      fSys[i][fNSys-1].mult = xsec_lumi/xsec ;
      fSys[i][fNSys-1].add  = fSys[i][fNSys-1].mult*fData[i]/100;
      fSys[i][fNSys-1].type = MULT;
      fSys[i][fNSys-1].name = "CMSLUMI12";

    }

  for(int i = 0; i < fNData; i++) 
    {
      delete[] sys1[i];
      delete[] sys2[i];
    }
  
  delete [] sys1;
  delete [] sys2;

  f1.close();
  f2.close();
  f3.close();

}
