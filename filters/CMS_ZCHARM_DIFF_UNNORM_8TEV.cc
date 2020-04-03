/*
*   Experiment: CMS 
*   Process: proton + proton -> Z + charm
*   Center of mass energy = 8 TeV
*   Intergrated luminosity  = 19.7 1/fb
*   Reference: https://arxiv.org/pdf/1711.02143.pdf
*
*   The raw data file is formatted as:
*   pT_min - pT_max - differential cross section - statistical uncertainty - systematic (%)
*   The systematics are as follows:
*   Luminosity, Branching c->l, Branching b->l, D±, D*(2010)±, Pile up, Missing transverse energy
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
  double pTmin, pTmax;         //pT ranges
  double MZ2 = pow(MZ,2.0);    //Z mass
  double s = 8000;             //LHC at 8TeV

  for(int i=0; i<fNData;i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    //Reading in an interpretation of each column
    lstream >> pTmin >> pTmax >> fData[i] >> fStat[i];
    
    fKin1[i] = 0.5*(pTmin + pTmax);	 // pT
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
  
 }
