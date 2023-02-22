/*
Inclusive total cross section for ttbar production @LHC CMS 5 TeV
There is a single point

LHC-CMS 5 TeV
---------------

emu events with b-tagged jets (L=27.4 1/pb)
[1711.03143]
sigma ttbar = 69.5 ± 6.1 (stat) +5.6 - 5.6 (syst) ± 1.6 (lumi) pb
 */

#include "CMSTTBARTOT5TEV.h"
#include <NNPDF/exceptions.h>

void CMSTTBARTOT5TEVFilter::ReadData()
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
      double sys1, sys2;

      getline(f1,line);
      istringstream lstream(line);
      lstream >> idum >> cme;

      double mt=172.5;

      fKin1[i] = 0.;
      fKin2[i] = mt*mt;          //top mass
      fKin3[i] = cme*1000;       //sqrt(s)

      lstream >> fData[i];       //central value
      lstream >> fStat[i];       //statistical uncertainty
      lstream >> sys1 >> sys2;   //Symmetric systematic uncertainty
     
      if (sys2 != -sys1)
          throw NNPDF::LogicException("CMSTTBARTOT5TEVFilter::ReadData", "systematics are not symmetric.");

      sys1 = sys1/fData[i]*100;

      fSys[i][0].mult = sys1;    //Symmetric systematic uncertainty
      lstream >> fSys[i][1].add; //Luminosity uncertainty

      fSys[i][0].add = fSys[i][0].mult*fData[i]/100;
      fSys[i][0].type = MULT;
      fSys[i][0].name = "UNCORR";

      fSys[i][1].mult = fSys[i][1].add/fData[i]*100;
      fSys[i][1].type = MULT;
      fSys[i][1].name = "CMSLUMI5";


    }

  f1.close();

}
