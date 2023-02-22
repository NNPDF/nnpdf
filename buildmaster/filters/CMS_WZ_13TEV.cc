/*
Reference:
   arXiv:     [1901.03428]
   hepdata:   n/a
   published:     JHEP 04 (2019) 122
Description:
   The WZ production cross section is measured in proton-proton collisions at a 
   centre-of-mass energy of 13 TeV using data collected with the CMS detector, 
   corresponding to an integrated luminosity of 35.9 fb-1. Normalised 
   differential cross section measurements are also presented with respect to 
   three variables: the Z boson transverse momentum pT, the leading jet pT, and 
   the M (W,Z) variable, defined as the invariant mass of the system composed 
   of the three leptons and the missing transverse momentum. The implementation 
   is based on Tabs. 9, 11 and 12 in the paper.
*/

#include "CMS_WZ_13TEV.h"

//CMS_WZ_13TEV_pTZ: combined pT spectrum
void CMS_WZ_13TEV_pTZFilter::ReadData()
{
  fstream f1;

  //Central values and breakdown of uncertainties
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/CMS_WZ_13TEV/pT.dat";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  for(int i=0; i<fNData; i++)
    {
      double ddum;
      string line;
      getline(f1,line);
      istringstream lstream(line);
      
      fKin1[i] = 0.;
      fKin2[i] = 0.;
      fKin3[i] = 13000.;    //c.m. energy
      lstream >> ddum 
	      >> ddum
	      >> fData[i] 
	      >> fStat[i]
	      >> fSys[i][0].add
	      >> fSys[i][1].add;
      for(int j=0; j<fNSys; j++)
	{
	  fSys[i][j].mult = fSys[i][j].add/fData[i]*100.;
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = "CORR";
	}

    }

  f1.close();

}

//CMS_WZ_13TEV_mTZ: combined mass spectrum
void CMS_WZ_13TEV_mTZFilter::ReadData()
{
  fstream f1;

  //Central values and breakdown of uncertainties
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/CMS_WZ_13TEV/mT.dat";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  for(int i=0; i<fNData; i++)
    {
      double ddum;
      string line;
      getline(f1,line);
      istringstream lstream(line);
      
      fKin1[i] = 0.;
      fKin2[i] = 0.;
      fKin3[i] = 13000.;    //c.m. energy
      lstream >> ddum 
	      >> ddum
	      >> fData[i] 
	      >> fStat[i]
	      >> fSys[i][0].add
	      >> fSys[i][1].add;
      for(int j=0; j<fNSys; j++)
	{
	  fSys[i][j].mult = fSys[i][j].add/fData[i]*100.;
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = "CORR";
	}

    }

  f1.close();

}

//CMS_WZ_13TEV_pTlead: combined pT leading spectrum
void CMS_WZ_13TEV_pTleadFilter::ReadData()
{
  fstream f1;

  //Central values and breakdown of uncertainties
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/CMS_WZ_13TEV/pTlead.dat";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  for(int i=0; i<fNData; i++)
    {
      double ddum;
      string line;
      getline(f1,line);
      istringstream lstream(line);
      
      fKin1[i] = 0.;
      fKin2[i] = 0.;
      fKin3[i] = 13000.;    //c.m. energy
      lstream >> ddum 
	      >> ddum
	      >> fData[i] 
	      >> fStat[i]
	      >> fSys[i][0].add
	      >> fSys[i][1].add;
      for(int j=0; j<fNSys; j++)
	{
	  fSys[i][j].mult = fSys[i][j].add/fData[i]*100.;
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = "CORR";
	}

    }

  f1.close();

}
