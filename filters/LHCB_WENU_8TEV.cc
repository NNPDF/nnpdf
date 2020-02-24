/*
Name_exp   : LHCB_WENU_8TEV
Reference  : Measurement of forward W → eν production in pp collisions 
             at √s = 8 TeV
ArXiv      : arXiv:1608.01484.
Published  : JHEP 1610 (2016) 030
Hepdata    : n/a
Description: A measurement of the cross-section for W → eν production in pp 
             collisions is presented using data corresponding to an integrated 
             luminosity of 2 fb^−1 collected by the LHCb experiment at a 
             centre-of-mass energy of √s = 8 TeV.

Implemented by TG February 2020. 
*/

#include "LHCb.h"
void  LHCB_WENU_8TEVFilter::ReadData()
{
  fstream fWp, fWm, fCorr;

  stringstream datafileWp("");
  datafileWp << dataPath() << "rawdata/" 
	     << fSetName << "/LHCB_WENU_8TEV_Wp.dat";
  fWp.open(datafileWp.str().c_str(), ios::in);

  if (fWp.fail()) {
    cerr << "Error opening data file " << datafileWp.str() << endl;
    exit(-1);
  }

  stringstream datafileWm("");
  datafileWm << dataPath() << "rawdata/" 
	     << fSetName << "/LHCB_WENU_8TEV_Wm.dat";
  fWm.open(datafileWm.str().c_str(), ios::in);

  if (fWm.fail()) {
    cerr << "Error opening data file " << datafileWm.str() << endl;
    exit(-1);
  }

  stringstream datafileCorr("");
  datafileCorr << dataPath() << "rawdata/" 
	       << fSetName << "/corr.dat";
  fCorr.open(datafileCorr.str().c_str(), ios::in);

  if (fCorr.fail()) {
    cerr << "Error opening data file " << datafileCorr.str() << endl;
    exit(-1);
  }

  string line;
  double MW2 = pow(MW,2.0);
  double s = 8;         //centre-of-mass energy     
  int ndata_Wp = 8;	//Wp points
  double etamin, etamax, sigma, fsr;
  std::vector<double> totsys(fNData);
  
  //Read Wp data
  for (int i = 0; i < 2; i++)
    getline(fWp,line);
  
  for (int i = 0; i < ndata_Wp; i++)
  {
    getline(fWp,line);                     
    istringstream lstream(line);         

    lstream >> etamin >> etamax;
    fKin1[i] = (etamin + etamax)/2.;         
    fKin2[i] = MW2;             
    fKin3[i] = s;               

    lstream >> sigma;        
    lstream >> fStat[i]; 
    lstream >> totsys[i];

    //Beam energy sys
    lstream >> fSys[i][fNSys-2].add;  
    fSys[i][fNSys-2].mult = fSys[i][fNSys-1].add/fData[i]*1e2;
    fSys[i][fNSys-2].type = MULT;       
    fSys[i][fNSys-2].name = "LHCBBEAM8TEV";	       

    //Lumi sys
    lstream >> fSys[i][fNSys-1].add;  
    fSys[i][fNSys-1].mult = fSys[i][fNSys-1].add/fData[i]*1e2;
    fSys[i][fNSys-1].type = MULT;		
    fSys[i][fNSys-1].name = "LHCBLUMI8TEV";     

    lstream >> fsr >> fsr;
  }

  //Read Wm data
  for (int i = 0; i < 2; i++)
    getline(fWm,line);
  
  for (int i = ndata_Wp; i < fNData; i++)
  {
    getline(fWm,line);                     
    istringstream lstream(line);         

    lstream >> etamin >> etamax;
    fKin1[i] = (etamin + etamax)/2.;         
    fKin2[i] = MW2;             
    fKin3[i] = s;               

    lstream >> sigma;        
    lstream >> fStat[i]; 
    lstream >> totsys[i];

    //Beam energy sys
    lstream >> fSys[i][fNSys-2].add;  
    fSys[i][fNSys-2].mult = fSys[i][fNSys-1].add/fData[i]*1e2;
    fSys[i][fNSys-2].type = MULT;
    fSys[i][fNSys-2].name = "LHCBBEAM8TEV";

    //Lumi sys
    lstream >> fSys[i][fNSys-1].add;  
    fSys[i][fNSys-1].mult = fSys[i][fNSys-1].add/fData[i]*1e2;
    fSys[i][fNSys-1].type = MULT;
    fSys[i][fNSys-1].name = "LHCBLUMI8TEV";

    lstream >> fsr >> fsr;
  }

  //Defining covariance matrix
  double** covmat = new double*[fNData];
  for (int i = 0; i < fNData; i++) 
    covmat[i] = new double[fNData];
 
  //Reading Covariance Matrix
  for (int i = 0; i < 2; i++) getline(fCorr,line);

  for (int i = 0; i < fNData; i++)
    { 
      getline(fCorr,line);
      for (int j = 0; j < i+1; j++) 
	{    
	  fCorr >> covmat[i][j];
	  //covmat[i][j] = covmat[i][j]*totsys[i]*totsys[j];  
	  covmat[j][i] = covmat[i][j];
	}
    }

  
  //Check
  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  cout << covmat[i][j] << "  " ;
	}
      cout << std::endl;
    }

  //Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
   {
     cerr << " in " << fSetName << " : cannot generate artificial systematics" << endl;
     exit(-1);
   }

  //Copy the artificial systematics in the fSys matrix
  for (int i = 0; i < fNData; i++)
    for (int l = 0; l < fNSys-3; l++)    
    {
      fSys[i][l].add  = syscor[i][l];
      fSys[i][l].mult = fSys[i][l].add/fData[i]*1e2;
      fSys[i][l].type = ADD;  
      fSys[i][l].name = "CORR";
    }

  fWp.close();
  fWm.close();
  fCorr.close();

  for(int i = 0; i < fNData; i++) 
    delete[] covmat[i];
  delete[] covmat;

  for(int i = 0; i < fNData; i++) 
    delete[] syscor[i];
  delete[] syscor;     
}





