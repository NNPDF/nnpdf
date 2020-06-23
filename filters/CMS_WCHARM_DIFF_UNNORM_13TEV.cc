/*
*   Experiment: CMS 
*   Process: proton + proton -> W + charm
*   Center of mass energy = 13 TeV
*   Intergrated luminosity  = 35.7 1/fb
*   Reference: https://cds.cern.ch/record/2314570/files/SMP-17-014-pas.pdf
*
*   The raw data file is formatted as:
*   eta_min - eta_max - differential cross section - statistical uncertainty - 
*   positive systematic (%) - negative systematic (%)
*   The systematics are broken into positive and negative components to allow 
*   for asymmetric uncertainties 
*   The systematics are as follows:
*   Luminosity, Tracking, Branching, Muons, Nsel determination, 
*   D*(2010)± kinematics, Bg normalization, p^T_miss, Pile Up, PDF, 
*   Secondary vertex, Fragmenation, Monte Carlo statistics
*   Systematici uncertainties are supplemented by a theoretical uncertainty
*   that takes into account missing NNLO corrections in the matrix element.
*   This uncertainty was estimated by means as the asymmetric envelope of the 
*   3pt renormalisation scale variation (Eq. 4.16 in 1906.10698).
*/

#include "CMS_WCHARM_DIFF_UNNORM_13TEV.h"

void CMS_WCHARM_DIFF_UNNORM_13TEVFilter::ReadData()
{
  // Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() 
       << "rawdata/" << fSetName << "/" << fSetName << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Starting filter
  string line;
  double etamin, etamax;        //eta ranges
  double MW2 = pow(MW,2.0);     //W mass
  double s = 13000;             //LHC at 13TeV
  double stmp, dtmp;
  //Systematics
  double fSystematics[2*fNSys];

  for(int i=0; i<fNData;i++)
    {
      getline(f1,line);
      istringstream lstream(line);

      //Reading in an interpretation of each column
      lstream >> etamin >> etamax >> fData[i] >> fStat[i];
      //Reading in the systematics
      for(int k=0; k<2*(fNSys-2);k++) // Factor of 2 because of +sys -sys format of rawdata
      {
        lstream >> fSystematics[k];
      }      
      lstream >> fSystematics[fNSys-2] >> fSystematics[fNSys-1];
      
      fData[i] = fData[i]*1000; // changing pb to fb for APPLgrid
      fStat[i] = fStat[i]*1000; // changing pb to fb for APPLgrid

      //Defining the kinematic variables
      fKin1[i] = (etamax + etamin)*0.5;    // eta
      fKin2[i] = MW2;                      // Mass W squared
      fKin3[i] = s;                        // sqrt(s)

      //Symmetrising the systematics
      for(int m=0;m < 2*(fNSys-2); m = m + 2)
      {
        symmetriseErrors(fSystematics[m],fSystematics[m+1],&stmp,&dtmp);
        fSys[i][m/2].mult=stmp;
        fSys[i][m/2].type = MULT;
        if(m == 0){
          fSys[i][0].name = "CMSLUMI13"; //We treat luminosity as a special case
          }
        else{
          fSys[i][m/2].name = "CORR";
        }
        //Add the additive uncertainties
        fSys[i][m/2].add = fSys[i][m/2].mult*fData[i]*1e-2;
      }

      fSys[i][fNSys-2].mult = fSystematics[fNSys-2]/sqrt(2.);
      fSys[i][fNSys-2].add  = fSys[i][fNSys-2].mult*fData[i]*1e-2;
      fSys[i][fNSys-2].type = MULT;
      fSys[i][fNSys-2].name = "SKIP";

      fSys[i][fNSys-1].mult = fSystematics[fNSys-1]/sqrt(2.);
      fSys[i][fNSys-1].add  = fSys[i][fNSys-1].mult*fData[i]*1e-2;
      fSys[i][fNSys-1].type = MULT;
      fSys[i][fNSys-1].name = "SKIP";
      
    }  

  f1.close();

}
