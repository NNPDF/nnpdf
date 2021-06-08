/*
Name_exp   : ATLAS_DY_2D_8TEV
Reference  : Measurement of the Drell–Yan triple-differentialcross section 
             in pp collisions at √s=8 TeV
ArXiv      : 1710.05167
Published  : JHEP 12 (2017) 059
Hepdata    : https://www.hepdata.net/record/ins1630886
Description: 
Measurement of the triple-differential cross section for Drell-Yan electron/muon
production. The measurement is performed for invariant masses of the lepton 
pairs, between 46 and 200 GeV using a sample of 20.2  fb−1 of pp collisions  
at a centre-of-mass energy of √s=8 TeV. The data are presented in bins of 
invariant mass, absolute dilepton rapidity and the angular variable cosθ∗
between the outgoing lepton and the incoming quark in the Collins–Soper frame.
Measurements considered here are integrated over the angular variable and are
taken from Table 5 of the Hepdata entry.
*/


#include "ATLAS_DY_2D_8TEV_LOWMASS.h"

void ATLAS_DY_2D_8TEV_LOWMASSFilter::ReadData()
{
  fstream f1;

  //Full breakdown of systematic ucnertainties
  stringstream datafile("");
  datafile << dataPath()
	  << "rawdata/ATLAS_DY_2D_8TEV_LOWMASS/HEPData-ins1630886-v3-Table_5.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read central value and uncertainties
  string line;
  for(int i=0; i<15; i++)
    {
      getline(f1,line);     
    }

  double ddum;
  char comma;

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      lstream >> ddum >> comma
	      >> ddum >> comma
	      >> ddum >> comma
	      >> fKin1[i] >> comma
	      >> ddum >> comma
	      >> ddum >> comma
	      >> fKin2[i] >> comma
	      >> ddum >> comma
	      >> ddum >> comma
	      >> fData[i] >> comma  //[pb/GeV]
	      >> fStat[i] >> comma
	      >> ddum;

      fKin2[i] = fKin2[i]*fKin2[i]; //[GeV2]
      fKin3[i] = 8000.;             //[GeV]

      //Convert to [fb]
      fData[i] *= 1000.;
      fStat[i] *= 1000.;

      //Systematic uncertainties
      for(int j=0; j<fNSys-1; j++)
	{
	  lstream >> comma >> fSys[i][j].add
		  >> comma >> ddum;

	  fSys[i][j].add *= 1000.; //Convert to [fb]
	  
	  fSys[i][j].mult = fSys[i][j].add/fData[i]*100.;
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = "CORR";
	}

      //Uncorrelated uncertainty
      fSys[i][276].mult = 1.9; //%
      fSys[i][276].add = fSys[i][276].mult/100. * fData[i];
      fSys[i][276].type = MULT;
      fSys[i][276].name = "UNCORR";

      //Uncorrelated uncertainties
      fSys[i][275].name = "UNCORR";      
    }

  f1.close();
}



