/*
Measurement of the inclusive single top production cross section in the t-channel @LHC CMS 7 TeV
There is a single point

LHC-CMS 7 TeV
---------------

Data corresponds to the muon (L=1.17 1/fb) and electron (L=1.56 1/fb) final states collected in 2011
[1209.4533]
sigma t-ch = 67.2 ± 3.7 (stat.) ± 3.0 (syst.) ± 3.5 (theor.) ± 1.5 (lumi.) pb

For more details of uncertainties taken into account in the above measurement, see Table 2 in the paper.
Note that the theoretical (theor.) uncertainty is based on the theory uncertainty on the background subtraction
and on the signal templates, and includes scale, matching, signal generator and PDF uncertainties.
*/

#include "CMS_SINGLETOP_TCH_TOT_7TEV.h"

void CMS_SINGLETOP_TCH_TOT_7TEVFilter::ReadData()
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
      int idum;
      double cme;

      getline(f1,line);
      istringstream lstream(line);
      lstream >> idum >> cme;
      
      fKin1[0] = 0.;
      fKin2[0] = Mt*Mt;          //top mass
      fKin3[0] = cme*1000;       //sqrt(s)

      lstream >> fData[0];       //central value
      lstream >> fStat[0];       //statistical uncertainty

      lstream >> fSys[0][0].add; //systematic uncertainty
      lstream >> fSys[0][1].add; //theoretical uncertainty
      lstream >> fSys[0][2].add; //luminosity uncertainty

      fSys[0][0].mult = fSys[0][0].add/fData[0]*100;
      fSys[0][0].type = MULT;
      fSys[0][0].name = "UNCORR";

      fSys[0][1].mult = fSys[0][1].add/fData[0]*100;
      fSys[0][1].type = MULT;
      fSys[0][1].name = "UNCORR";

      fSys[0][2].mult = fSys[0][2].add/fData[0]*100;
      fSys[0][2].type = MULT;
      fSys[0][2].name = "CMSLUMI11";

  f1.close();

}
