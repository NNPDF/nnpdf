/* 
   This file implements the Z production subset of the ATLASWZRAP36PB data set.
   This is required to separate NC DY from CC DY.
   Implemented by ERN June 2023.
*/

#include "ATLASZRAP36PB.h"

void ATLASZRAP36PBFilter::ReadData()
{
  // Opening files
  fstream fWZ;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/ATLASWZRAP36PB/ATLAS-36pb-Zrap.data";
  fWZ.open(datafile.str().c_str(), ios::in);

  if (fWZ.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  const double lcorr = 1.0187; // correction factor due to luminosity upgrade
  const int ndataWZ = 8;  //Number of data for W+, W- and Z respectively
  const double convfac = lcorr*1000.; // Must multiply from pb to fb
  const double MWZ2 = pow(MZ,2.0);   //Mass squared of W (+ and -) and Z
  
  string line;
  double etamin,etamax,tmp;
  
  // rapidity
  getline(fWZ,line);
  istringstream lstream(line);
  for (int i = 0; i < ndataWZ; i++)
    {
      lstream >> etamin >> etamax;
      fKin1[i] = etamin + (etamax-etamin)*0.5;
    }

  //Z
  for (int i = 0; i < ndataWZ; i++)
    fKin2[i] = MWZ2;
  
  // sqrt(s)
  for (int i = 0; i < ndataWZ; i++)
    fKin3[i] = 7000;
  
  // obs
  getline(fWZ,line);
  istringstream lstream2(line);
  for (int i = 0; i < ndataWZ; i++)
    {
      lstream2 >> fData[i];
      fData[i] *= convfac;
    }

  // stat (%, converted later)
  getline(fWZ,line);
  istringstream lstream3(line);
  for (int i = 0; i < ndataWZ; i++)
    lstream3 >> fStat[i];
  
  // uncorrelated sys
  getline(fWZ,line);
  istringstream lstream4(line);
  for (int i = 0; i < ndataWZ; i++)
    {
      lstream4 >> fSys[i][0].mult;
      fSys[i][0].name = "UNCORR";
    }
  
  // total correlated sys (unused)
  getline(fWZ,line);
  
  // total uncertainty (unused)
  getline(fWZ,line);
  
  // correlated systematics
  for (int isys = 2; isys < fNSys; isys++)  //2 to skip uncorr and lumi
    {
      getline(fWZ,line);
      istringstream lstream(line);
      lstream >> tmp;
      for (int i = 0; i < ndataWZ; i++)
	{
	  lstream >> fSys[i][isys].mult;
	  ostringstream sysname;
	  sysname << "ATLASWZRAP36PB_" << isys-2;
	  fSys[i][isys].name = sysname.str();
	}
    }
  
  // luminosity: 3.4%
  for (int i = 0; i < ndataWZ; i++)
    {
      fSys[i][1].mult = 3.5;
      fSys[i][1].name = "ATLASLUMI10";
    }
  
  // Convert additive uncertainties to absolute form
  for (int i = 0; i < fNData; i++)
    {
      fStat[i] *= fData[i]*1e-2;
      for(int l = 0; l < fNSys; l++)
	{
	  fSys[i][l].type = MULT; // All systematics multiplicative
	  fSys[i][l].add = fSys[i][l].mult*fData[i]*1e-2;
	}
    }
  
  fWZ;
}
