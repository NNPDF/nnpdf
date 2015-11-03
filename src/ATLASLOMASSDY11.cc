/*
 *
 * ATLASLOMASSDY11 - ATLAS low mass Drell-Yan 2011
 *
 * ATLAS differential cross-section for Z/gamma* -> ll production,
 * as function of the dilepton invariant mass [26 GeV-66 GeV ].
 * LHC 7 TeV, 1.6 fb^-1
 * Reference: JHEP1406,112(2014), [arXiv:1404.1212]
 * HepData: http://hepdata.cedar.ac.uk/view/ins1288706
 *
 */

#include "ATLASLOMASSDY11.h"

/*
 * ATLASLOMASSDY11 dataset
 *
 * d sigma/dM_{LL} [pb / GeV], |eta^{gamma}|
 */

void ATLASLOMASSDY11Filter::ReadData()
{
  // Opening files
  fstream f1;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/ATLASLOMASSDY11/ATLASLOMASSDY11.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double mbin[fNData+1];//, shift[fNData], stmp, dtmp;

  // Read data
  string line;

  //Skip two description lines
  getline(f1, line);
  getline(f1, line);

  // Filter data
  for (int idat = 0; idat < fNData; idat++) {

    getline(f1, line);

    double fStatP, fStatM;
    double fUncorrP, fUncorrM;
    double fCorr[fNSys-2];
    double fATLAS2011Luminosity;

    istringstream lstream(line);
    lstream >> mbin[idat] >> mbin[idat+1] >> fData[idat]
	    >> fStatP >> fStatM
	    >> fUncorrP >> fUncorrM
      >> fCorr[1] >> fCorr[2] >> fCorr[3] >> fCorr[4] >> fCorr[5]
	    >> fCorr[6] >> fCorr[7] >> fCorr[8] >> fCorr[9] >> fCorr[10]
	    >> fCorr[11] >> fCorr[12] >> fCorr[13]
 	    >> fATLAS2011Luminosity;

    // Convert from pb (paper) to fb (APPLrid)

    fData[idat] *= 1000.;
    fStatP      *= 1000.;
    fStatM      *= 1000.;
    fUncorrP    *= 1000.;
    fUncorrM    *= 1000.;

    // Statisltical Uncenrtainty (absolute value in the data file)
    fStat[idat] = 0.5 * (fStatP - fStatM);

    // Total Uncorrelated systematics (absolute value in the data file)
    fSys[idat][0].add  = 0.5 * (fUncorrP - fUncorrM);
    fSys[idat][0].mult = fSys[idat][0].add/fData[idat]*1e2;
    fSys[idat][0].type = ADD;
    fSys[idat][0].name = "UNCORR";

    // Correlated systematics (given in % in the data file)
    for ( int k = 1; k < fNSys-1; k++ )
      {
	      fSys[idat][k].mult = fCorr[k];
	      fSys[idat][k].add  = fSys[idat][k].mult*fData[idat]*1e-2;
	      fSys[idat][k].type = ADD;
	      fSys[idat][k].name = "CORR";
      }

    // ATLAS 2011 Luminosity (given in % in the data file)
    fSys[idat][fNSys-1].mult = fATLAS2011Luminosity;
    fSys[idat][fNSys-1].add  = fSys[idat][fNSys-1].mult*fData[idat]*1e-2;
    fSys[idat][fNSys-1].type = MULT;
    fSys[idat][fNSys-1].name   = "ATLASLUMI11";

    // Kinematic variables
    fKin1[idat] = 0.0;                                // Dummy
    fKin2[idat] = (mbin[idat] + mbin[idat+1]) * 0.5;  // Avg. M_ll of each bin
    fKin3[idat] = 7E3;                                // LHC 7 TeV

  }

  f1.close();

}
