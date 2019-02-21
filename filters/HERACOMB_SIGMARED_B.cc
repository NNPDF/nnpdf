/*
Name_exp  : HERACOMB_SIGMA_C
Reference : Combination and QCD analysis of charm and beauty
            production cross-section measurements in deep inelastic ep
            scattering at HERA
ArXiv     : arXiv:1804.01019
Published : Eur.Phys.J. C78 (2018) no.6, 473

Combination of H1 and ZEUS measurements of open bottom production reduced 
cross sections in deep-inelastic ep scattering at HERA. The data files 
contain information on the statistical error, the uncorrelated error, and
the full breakdown of 167 bin-by-bin correlated errors. The statistical and 
the uncorrelated errors are added in quadrature. Results are given for a 
centre-of-mass energy sqrts=318 GeV.

This data set supersedes H1HERAF2B and ZEUSHERAF2B

Implemented by ERN January 2019.
*/

#include "HERACOMB_SIGMARED_B.h"

void HERACOMB_SIGMARED_BFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
	   << fSetName << "/d18-037.tableBeauty.txt";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  //Starting filter
  string line;
  //Header
  for (int i=0; i<36; i++) getline(f1,line);
  //Data
  for (int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      lstream >> fKin2[i];                                //Q2 [GeV2]
      lstream >> fKin1[i];                                //x
      double sqrts = 318.0;                               //centre-of-mass energy [GeV]
      fKin3[i] = fKin2[i] / fKin1[i] / (sqrts*sqrts);     //inelasticity
      lstream >> fData[i];                                //observable
      double stat, uncorr;
      lstream >> stat >> uncorr;
      stat   *= fData[i]*1e-2;
      uncorr *= fData[i]*1e-2;
      fStat[i] = pow(stat*stat+uncorr*uncorr,0.5);        //total uncorrelated error
      for(int l=0; l<fNSys; l++)
	{
	  lstream >> fSys[i][l].mult;                     //percentage systematic
	  fSys[i][l].add = fSys[i][l].mult*fData[i]*1e-2; //absolute systematic
	  fSys[i][l].type = MULT;
	  fSys[i][l].name = "CORR";
	}
    }
  f1.close();
}
