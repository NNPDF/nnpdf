/**
 * ATLAS_WP_7TEV
 * (see also ATLASWZRAP11.cc)
 * Implementation of ATLAS 2011 W plus rapidity data
 *
 * This filter contains the implementation of dataset: ATLAS_WP_7TEV,
 * Corresponding to the "central" selection cuts. The specific cuts 
 * can be found at the top of page eight of https://arxiv.org/pdf/1612.03016.pdf.
 *
 * Archived as: https://arxiv.org/abs/1612.03016
 * Published as: https://link.springer.com/article/10.1140%2Fepjc%2Fs10052-017-4911-9
 * HEPData entry: https://www.hepdata.net/record/ins1502620
 *
 */

#include "ATLAS.h"


// Central selection
void ATLAS_WP_7TEVFilter::ReadData()
{

  cout << "********** WARNING: Converting pb to fb to match ApplGrid output ********" << endl;

  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << "ATLAS_WP_7TEV/wplus.dat";

  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Dummy string  
  std::string dummy;

  // Read the names of systematics
  std::string *sysNames = new std::string[fNSys];
  sysNames[0] = "UNCORR";
  for (int i=0; i<5; i++) f1 >> dummy;
  for (int i=0; i<fNSys; i++) 
    f1 >> sysNames[i];
 
  const int nBins = 5;
  const int ndataWZ[nBins] = {11};  // Data thresholds for W+, W- and (Z_low, Z_peak, Z_high) respectively
  const double MWZ2[nBins]= {pow(MW,2.0), pow(MW,2.0), pow(56.0,2.0), pow(91.0,2.0), pow(133.0,2.0)};   //Mass squared of W (+ and -) and (Z_low, Z_peak, Z_high)

  int low_bin = 0;
  for (int b = 0; b < nBins; b++)
  {
    for (int i = low_bin; i < ndataWZ[b]; i++)
    {
      double etamin, etamax;

      // Kinematics
      f1 >> dummy; f1 >> etamin; f1 >> etamax;
      fKin1[i] = etamin + (etamax - etamin)/2.0;
      fKin2[i] = MWZ2[b];
      fKin3[i] = 7000;

      // Observable
      f1 >> fData[i];
      fData[i] *= 1000; // pb -> fb

      // Statistical errors - percentage with respect the observable
      f1 >> fStat[i];
      fStat[i] *= fData[i]*1e-2;

      // Correlated systematic errors
      for (int l = 0; l < fNSys; l++)
      {
        f1 >> fSys[i][l].mult;
        fSys[i][l].type = MULT;
        fSys[i][l].name = sysNames[l];
      }

      // Additive errors   
      for (int l = 0; l < fNSys; l++)
        fSys[i][l].add = fSys[i][l].mult*fData[i]*1e-2;
    }
    // Update lowest point in bin
    low_bin = ndataWZ[b];
  }

  delete[] sysNames;
  
  f1.close();
}
