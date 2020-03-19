/*
*   Experiment: ATLAS 
*   Process: proton + proton -> W + charm
*   Center of mass energy = 7 TeV
*   Intergrated luminosity  = 4.6 1/fb
*   Reference: https://link.springer.com/content/pdf/10.1007/JHEP05(2014)068.pdf
*   HEP Data: https://www.hepdata.net/record/ins1282447
*             with data taken from Tables 1, 12, and 13
*
*   The two raw data file correspond to c and cbar jets.
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
  datafileWM << dataPath() 
       << "rawdata/" << fSetName << "/" << fSetName << "_WM.data";
  f2.open(datafileWM.str().c_str(), ios::in);
  if (f2.fail()) 
  {
    cerr << "Error opening data file " << datafileWM.str() << endl;
    exit(-1);
  }

  //Starting filter
  string line;
  double etamin, etamax;       //rapidity binning
  double MW2 = pow(MW, 2.0);   //W mass
  double s = 7000;             //LHC at 7TeV
  double fSystematics[fNSys];
  for(int i=0; i<fNData;i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    //Reading in an interpretation of each column
    lstream >> etamin >> etamax >> fData[i] >> fStat[i];
    for(int k=0; k<fNSys; k++)
    {
      lstream >> fSystematics[k];
    }      

    fData[i] = fData[i]*1000; // changing pb to fb for APPLgrid
    fStat[i] = fStat[i]*1000; // changing pb to fb for APPLgrid
    //Defining the kinematic variables
    fKin1[i] = (etamin + etamax)*0.5;    // eta
    fKin2[i] = MW2;                      // Mass W squared
    fKin3[i] = s;                        // sqrt(s)

    //Reading in the systematics
    for(int k=0; k<fNSys; k++)
    {
      fSys[i][k].mult = fSystematics[k];
      std::cout << fSys[i][k].mult << std::endl;
      fSys[i][k].type = MULT;
      if(k == 0){
        fSys[i][k].name = "UNCORR"; //We treat luminosity as a special case
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
