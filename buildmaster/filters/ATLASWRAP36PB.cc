/* 
   This file implements the W production subset of the ATLASWZRAP36PB data set.
   This is required to separate CC DY from NC DY.
   Implemented by ERN June 2023.
*/

#include "ATLASWRAP36PB.h"

void ATLASWRAP36PBFilter::ReadData()
{
  // Opening files
  fstream fWZ[2];

  stringstream datafile("");
  datafile << dataPath() << "rawdata/ATLASWZRAP36PB/ATLAS-36pb-Wplrap.data";
  fWZ[0].open(datafile.str().c_str(), ios::in);

  if (fWZ[0].fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream datafile2("");
  datafile2 << dataPath() << "rawdata/ATLASWZRAP36PB/ATLAS-36pb-Wmlrap.data";
  fWZ[1].open(datafile2.str().c_str(), ios::in);

  if (fWZ[1].fail()) {
    cerr << "Error opening data file " << datafile2.str() << endl;
    exit(-1);
  }

  // Starting filter
  const double lcorr = 1.0187; // correction factor due to luminosity upgrade
  const int ndataWZ[2] = {11,11};  //Number of data for W+, W-
  const double convfac = lcorr*1000.; // Must multiply from pb to fb
  const double MWZ2[2] = {pow(MW,2.0), pow(MW,2.0)};   //Mass squared of W (+ and -) and Z

  string line;
  int idat = 0;
  double etamin,etamax,tmp;

  for (int iWZ = 0; iWZ < 2; iWZ++)
  {
    // rapidity
    getline(fWZ[iWZ],line);
    istringstream lstream(line);
    for (int i = 0; i < ndataWZ[iWZ]; i++)
    {
      lstream >> etamin >> etamax;
      fKin1[idat+i] = etamin + (etamax-etamin)*0.5;
    }

    // M_W
    for (int i = 0; i < ndataWZ[iWZ]; i++)
      fKin2[idat+i] = MWZ2[iWZ];

    // sqrt(s)
    for (int i = 0; i < ndataWZ[iWZ]; i++)
      fKin3[idat+i] = 7000;

    // obs
    getline(fWZ[iWZ],line);
    istringstream lstream2(line);
    for (int i = 0; i < ndataWZ[iWZ]; i++)
      {
        lstream2 >> fData[idat+i];
        fData[idat+i] *= convfac;
      }

    // stat (%, converted later)
    getline(fWZ[iWZ],line);
    istringstream lstream3(line);
    for (int i = 0; i < ndataWZ[iWZ]; i++)
      lstream3 >> fStat[idat+i];

    // uncorrelated sys
    getline(fWZ[iWZ],line);
    istringstream lstream4(line);
    for (int i = 0; i < ndataWZ[iWZ]; i++)
    {
      lstream4 >> fSys[idat+i][0].mult;
      fSys[idat+i][0].name = "UNCORR";
    }

    // total correlated sys (unused)
    getline(fWZ[iWZ],line);

    // total uncertainty (unused)
    getline(fWZ[iWZ],line);

    // correlated systematics
    for (int isys = 2; isys < fNSys; isys++)  //2 to skip uncorr and lumi
    {
      getline(fWZ[iWZ],line);
      istringstream lstream(line);
      lstream >> tmp;
      for (int i = 0; i < ndataWZ[iWZ]; i++)
      {
        lstream >> fSys[idat+i][isys].mult;
	ostringstream sysname;
	sysname << "ATLASWZRAP36PB_" << isys-2;
	fSys[idat+1][isys].name = sysname.str();
      }
    }

    // luminosity: 3.4%
    for (int i = 0; i < ndataWZ[iWZ]; i++)
    {
      fSys[idat+i][1].mult = 3.5;
      fSys[idat+i][1].name = "ATLASLUMI10";
    }

    idat+=ndataWZ[iWZ];
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


  fWZ[0].close();
  fWZ[1].close();
}
