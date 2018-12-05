/*
*   Experiment: ATLAS 
*   Process: proton + proton -> W + charm
*   Center of mass energy = 7 TeV
*   Intergrated luminosity  = 4.6 1/fb
*   Reference: https://arxiv.org/abs/1402.6263
*   HEP Data: http://hepdata.cedar.ac.uk/view/ins1282447
*
*   The raw data file is formatted as:
*   
*   The systematics are as follows:
*   
*/
#include "ATLAS_WCHARM_TOT_UNNORM_7TEV.h"
void ATLAS_WCHARM_TOT_UNNORM_7TEVFilter::ReadData()
{
  // Opening files
  fstream f1, f2;
  stringstream datafileWP("");
  datafileWP << dataPath() 
       << "rawdata/" << fSetName << "/" << fSetName << "_WP.data";
  f1.open(datafileWP.str().c_str(), ios::in);
  if (f1.fail()) 
  {
      cerr << "Error opening data file " << datafileWP.str() << endl;
      exit(-1);
  }

  stringstream datafileWM("");
  datafileWP << dataPath() 
       << "rawdata/" << fSetName << "/" << fSetName << "_WM.data";
  f2.open(datafileWM.str().c_str(), ios::in);
  if (f2.fail()) 
  {
      cerr << "Error opening data file " << datafileWM.str() << endl;
      exit(-1);
  }
  
  //Starting filter
  string line;
  double etamin, etamax;          //pT ranges
  double MW2 = pow(MW,2.0);     //W mass
  double s = 7000;             //LHC at 7TeV
  for(int i=0; i<fNData;i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    //Reading in an interpretation of each column
    lstream >> pTmin >> pTmax >> fData[i] >> fStat[i];
     fData[i] = fData[i]*1000; // changing pb to fb for APPLgrid
    fStat[i] = fStat[i]*1000; // changing pb to fb for APPLgrid
    //Defining the kinematic variables
    fKin1[i] = (pTmin + pTmax)*0.5;    // eta
    fKin2[i] = MZ2;                      // Mass W squared
    fKin3[i] = s;                        // sqrt(s)
     //Reading in the systematics
    for(int k=0; k<fNSys;k++) // Factor of 2 because of +sys -sys format of rawdata
    {
      lstream >> fSys[i][k].mult;
      fSys[i][k].type = MULT;
      if(k == 0){
        fSys[i][k].name = "CMSLUMI12"; //We treat luminosity as a special case
      }
      else{
        fSys[i][k].name = "CORR";
      }
       //Add the additive uncertainties
      fSys[i][k].add = fSys[i][k].mult*fData[i]*1e-2;
    }      
     
  }  
   f1.close();
} 