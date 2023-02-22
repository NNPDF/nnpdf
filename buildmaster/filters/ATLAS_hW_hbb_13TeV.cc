/*
Measurement of WH, H->bbar production as a function of the vector-boson
transverse momentum in 13 TeV ppp collisions with the ATLAS detector
JHEP 1905 (2019) 141
[arXiv:1903.04618]
 */

#include "ATLAS_hW_hbb_13TeV.h"

void ATLAS_hW_hbb_13TeVFilter::ReadData()
{
  // Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/" << fSetName << "/" << fSetName << ".txt";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Starting filter
  for(int i=0; i<fNData;i++)
    {
      string line;
      getline(f1,line);
      istringstream lstream(line);
      double ddum;
      
      fKin1[i] = 0.;
      fKin2[i] = MH*MH;          //Higgs mass squared
      fKin3[i] = 13000;          //sqrt(s) [TeV]

      lstream >> fData[i];       //central value
      lstream >> ddum;           //total uncertainty
      lstream >> fStat[i];       //statistical uncertainty
      lstream >> fSys[i][0].add; //systematic uncertainty
      lstream >> fSys[i][1].add; //luminosity uncertainty
      lstream >> fSys[i][2].add; //beam energy uncertainty
      
      fSys[i][0].mult = fSys[i][0].add/fData[i]*100;
      fSys[i][0].type = MULT;
      fSys[i][0].name = "CORR";      

      fSys[i][1].mult = fSys[i][1].add/fData[i]*100;
      fSys[i][1].type = MULT;
      fSys[i][1].name = "CORR";
      
      fSys[i][2].mult = fSys[i][2].add/fData[i]*100;
      fSys[i][2].type = MULT;
      fSys[i][2].name = "CORR";
      
    }  

  f1.close();

}
