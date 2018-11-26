/*
*   Experiment: CMS 
*   Process: proton + proton -> Z + charm
*   Center of mass energy = 8 TeV
*   Intergrated luminosity  = 19.7 1/fb
*   Reference: https://arxiv.org/pdf/1711.02143.pdf
*
*   NEED TO CHANGE THIS
*   The raw data file is formatted as:
*   eta_min - eta_max - differential cross section - statistical uncertainty - positive systematic (%) - negative systematic (%)
*   The systematics are broken into positive and negative components to allow for asymmetric uncertainties 
*   The systematics are as follows:
*   Luminosity, Tracking, Branching, Muons, Nsel determination, D*(2010)± kinematics,
*   Bg normalization, p^T_miss, Pile Up, PDF, Secondary vertex, Fragmenation, Monte Carlo statistics
*/
 #include "CMS_ZCHARM_DIFF_UNNORM_8TEV.h"
 void CMS_ZCHARM_DIFF_UNNORM_8TEVFilter::ReadData()
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
  double pTmin, pTmax;          //pT ranges
  double MZ2 = pow(MZ,2.0);     //Z mass
  double s = 13000;             //LHC at 8TeV

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