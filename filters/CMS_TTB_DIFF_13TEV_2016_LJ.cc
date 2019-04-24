/*
Experiment: CERN-LHC-CMS (CMS)
Preprinted as CERN-EP-2018-039
Preprinted as CMS-TOP-17-002
Archived as: ARXIV:1803.08856
Published in Phys. Rev. D97 (2018) no 11, 112003
Measurement of differential cross sections for the production of top quark pairs
and of additional jets in lepton+jets events from pp collisions at sqrt(s)= 13
TeV

differential in the following variables are implemented:

1U) unnormalised top quark transverse momentum; 
    
Raw data and covariance matrices are from HepData:
https://www.hepdata.net/record/ins1663958

1U) TABS 182-183 HepData; Fig 11, TAB 3 1803.08856

*/
 
#include "CMS_TTB_DIFF_13TEV_2016_LJ.h"

//U - UNNORMALISED distributions

//1U) Distribution differential in top quark transverse momentum
void  CMS_TTB_DIFF_13TEV_2016_LJ_TPTFilter::ReadData()
{
  //Opening files
  fstream f1, f2;

  //Data values
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/CMS_TTB_DIFF_13TEV_2016_LJ_TPT/HEPData-ins1663958-v2-Table_182.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/CMS_TTB_DIFF_13TEV_2016_LJ_TPT/HEPData-ins1663958-v2-Table_183.csv";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read central values
  string line;
  for(int i=0; i<11; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<fNData; i++)
    {
      double pt_top, dud;
      char comma;
      
      getline(f1,line);
      istringstream lstream(line);
      lstream >> pt_top >> comma 
	      >> dud >> comma 
	      >> dud >> comma  
	      >> fData[i] >> comma
	      >> dud >> comma 
	      >> dud >> comma 
	      >> dud >> comma 
	      >> dud;

      fKin1[i] = pt_top;  //pTt
      fKin2[i] = Mt*Mt;   
      fKin3[i] = 13000;  //sqrt(s) [GeV] 
      fStat[i] = 0.;
    }

  //Read covariance matrix
  for(int i=0; i<11; i++)
    {
      getline(f2,line);
    }

  //Create covmat of correct dimensions
  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
    {
      covmat[i] = new double[fNData];
      
      for(int j=0; j<fNData; j++)
	{
	  double row, col;
	  char comma;

	  getline(f2,line);
	  istringstream lstream(line);
	  lstream >> row >> comma >> col >> comma  >> covmat[i][j];
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
	  fSys[i][j].mult  = fSys[i][j].add*1e2/fData[i];
	  fSys[i][j].type = ADD; // Should this be MULT? e.g see CMSTTBARDIFF13TEV2
	  fSys[i][j].name = "CORR";
	}
    }
  
  delete [] covmat; 
  delete [] syscor;
  
  f1.close();
  f2.close();

}
