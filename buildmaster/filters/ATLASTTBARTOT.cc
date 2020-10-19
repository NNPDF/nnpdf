/*
Inclusive total cross section for ttbar production @LHC ATLAS
There are three points

LHC-ATLAS 7 TeV
---------------

emu events with b-tagged jets (L=4.6 1/fb)
Eur.Phys.J. C74 (2014) no.10, 3109
sigma ttbar = 182.9 ± 3.1 (stat) ± 4.2 (syst) ± 3.6 (lumi) ± 3.3 (beam) pb
 
LHC-ATLAS 8 TeV
---------------

emu events with b-tagged jets (L=20.3 1/fb)
Eur.Phys.J. C74 (2014) no.10, 3109
sigma ttbar = 242.4 ± 1.7 (stat) ± 5.5 (syst) ± 7.5 (lumi) ± 4.2 (beam) pb

LHC-ATLAS 13 TeV
----------------

emu events with b-tagged jets (L=3.2 1/fb)
[1606.02699]
sigma ttbar = 818 ± 8 (stat) ± 27 (syst) ± 19 (lumi) ± 12 (beam) pb 
 */

#include "ATLASTTBARTOT.h"

void ATLASTTBARTOTFilter::ReadData()
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
  for(int i=0; i<fNData;i++)
    {
      string line;
      int idum;
      double cme;
      getline(f1,line);
      istringstream lstream(line);
      lstream >> idum >> cme;
      
      fKin1[i] = 0.;
      fKin2[i] = Mt*Mt;             //top mass
      fKin3[i] = cme*1000;       //sqrt(s)

      lstream >> fData[i];       //central value
      lstream >> fStat[i];       //statistical uncertainty
      lstream >> fSys[i][0].add; //systematic uncertainty
      lstream >> fSys[i][1].add; //luminosity uncertainty
      lstream >> fSys[i][2].add; //beam energy uncertainty
      
      fSys[i][0].mult = fSys[i][0].add/fData[i]*100;
      fSys[i][0].type = MULT;
      fSys[i][0].name = "UNCORR";      

      fSys[i][1].mult = fSys[i][1].add/fData[i]*100;
      fSys[i][1].type = MULT;
      fSys[i][1].name = "UNCORR";
      
      fSys[i][2].mult = fSys[i][2].add/fData[i]*100;
      fSys[i][2].type = MULT;
      fSys[i][2].name = "UNCORR";
      
    }  

  f1.close();

}
