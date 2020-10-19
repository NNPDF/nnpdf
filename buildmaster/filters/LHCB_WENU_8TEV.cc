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
Note:        The covaraince matrix provided in the paper is not positive
             semi-definite, therefore individual measurements of W+, W- 
	     rapidity distributions are not implemented.

Implemented by TG February 2020. 
*/

#include "LHCb.h"


//W+ to W− cross-section ratio in bins of electron pseudorapidity
void  LHCB_WENU_8TEV_RFilter::ReadData()
{
  fstream fR;

  stringstream datafileR("");
  datafileR << dataPath() << "rawdata/" 
	     << fSetName << "/LHCB_WENU_8TEV_R.dat";
  fR.open(datafileR.str().c_str(), ios::in);

  if (fR.fail()) {
    cerr << "Error opening data file " << datafileR.str() << endl;
    exit(-1);
  }

  string line;
  double MW2 = pow(MW,2.0);
  double s = 8000.;         //centre-of-mass energy [GeV]     
  double etamin, etamax;
  
  //Read R data
  for (int i = 0; i < 2; i++)
    getline(fR,line);
  
  for (int i = 0; i < fNData; i++)
  {
    getline(fR,line);                     
    istringstream lstream(line);         

    lstream >> etamin >> etamax;
    fKin1[i] = (etamin + etamax)/2.;         
    fKin2[i] = MW2;             
    fKin3[i] = s;               

    lstream >> fData[i];        
    lstream >> fStat[i]; 

    //Uncorrelated systematic
    lstream >> fSys[i][0].add;
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = MULT;       
    fSys[i][0].name = "UNCORR";

    //Beam energy sys
    lstream >> fSys[i][1].add;  
    fSys[i][1].mult = fSys[i][1].add/fData[i]*1e2;
    fSys[i][1].type = MULT;       
    fSys[i][1].name = "LHCBBEAM8TEV";	       

  }
}

//W boson production charge asymmetry in bins of electron pseudorapidity
void  LHCB_WENU_8TEV_AFilter::ReadData()
{
  fstream fA;

  stringstream datafileA("");
  datafileA << dataPath() << "rawdata/" 
	     << fSetName << "/LHCB_WENU_8TEV_A.dat";
  fA.open(datafileA.str().c_str(), ios::in);

  if (fA.fail()) {
    cerr << "Error opening data file " << datafileA.str() << endl;
    exit(-1);
  }

  string line;
  double MW2 = pow(MW,2.0);
  double s = 8000.;         //centre-of-mass energy [GeV]    
  double etamin, etamax;
  
  //Read A data
  for (int i = 0; i < 2; i++)
    getline(fA,line);
  
  for (int i = 0; i < fNData; i++)
  {
    getline(fA,line);                     
    istringstream lstream(line);         

    lstream >> etamin >> etamax;
    fKin1[i] = (etamin + etamax)/2.;         
    fKin2[i] = MW2;             
    fKin3[i] = s;               

    lstream >> fData[i];        
    lstream >> fStat[i]; 
    //Convert from percentage to absolute values
    fData[i] /= 100.;
    fStat[i] /= 100.;

    //Uncorrelated systematic
    lstream >> fSys[i][0].add;
    fSys[i][0].add /= 100.;
    fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
    fSys[i][0].type = MULT;       
    fSys[i][0].name = "UNCORR";

    //Beam energy sys
    lstream >> fSys[i][1].add;  
    fSys[i][1].add /= 100.;
    fSys[i][1].mult = fSys[i][1].add/fData[i]*1e2;
    fSys[i][1].type = MULT;       
    fSys[i][1].name = "LHCBBEAM8TEV";	       

  }
}





