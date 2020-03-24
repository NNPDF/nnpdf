/*
*   Experiment: ATLAS *   Process: proton + proton -> W+ + charm
*   Center of mass energy = 7 TeV
*   Intergrated luminosity  = 4.6 1/fb
*   Reference: https://link.springer.com/content/pdf/10.1007/JHEP05(2014)068.pdf
*   HEP Data: https://www.hepdata.net/record/ins1282447
*             with data taken from Tables 1 and 13
*
*
*/
#include "ATLAS_WPCHARM_TOT_UNNORM_7TEV.h"
void ATLAS_WPCHARM_TOT_UNNORM_7TEVFilter::ReadData()
{
  // Opening files
  fstream f1;
  stringstream datafileWP("");
  datafileWP << dataPath()
       << "rawdata/" << fSetName << "/" << fSetName << ".data";
  f1.open(datafileWP.str().c_str(), ios::in);
  if (f1.fail())
  {
    cerr << "Error opening data file " << datafileWP.str() << endl;
    exit(-1);
  }

  //Starting filter
  string line;
  double etamin, etamax;       //rapidity binning
  double MW2 = pow(MW, 2.0);   //W mass
  double s = 7000;             //LHC at 7TeV

  for(int i=0; i<fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    //Reading in an interpretation of each column
    lstream >> etamin >> etamax >> fData[i] >> fStat[i];

    fData[i] = fData[i]*1000; // changing pb to fb for APPLgrid
    fStat[i] = fStat[i]*1000; // changing pb to fb for APPLgrid
    //Defining the kinematic variables
    fKin1[i] = (etamin + etamax)*0.5;    // eta
    fKin2[i] = MW2;                      // Mass W squared
    fKin3[i] = s;                        // sqrt(s)

    //Reading in the systematics
    for(int k=0; k<fNSys; k++)
    {
      lstream >> fSys[i][k].mult;
      fSys[i][k].type = MULT;
      fSys[i][k].add = fSys[i][k].mult*fData[i]/100;
      if(k == 0){
        fSys[i][k].name = "UNCORR";
      }
      else{
        fSys[i][k].name = "CORR";
      }
    }
  }
  f1.close();
}
