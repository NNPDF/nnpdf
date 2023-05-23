/**
 * ATLAS_ZRAPPC_7TEV
 * (see also ATLASWZRAP11.cc)
 * Implementation of ATLAS 2011 Z peak rapidity data (central region)
 * The data is in the combined (electron + muon) channel
 *
 * This filter contains the implementation of the ATLAS_ZRAPPC_7TEV set,
 * cuts: Central selection. The specific cuts used to find this
 * set can be found at the top of page eight of https://arxiv.org/pdf/1612.03016.pdf.
 *
 * Archived as: https://arxiv.org/abs/1612.03016
 * Published as: https://link.springer.com/article/10.1140%2Fepjc%2Fs10052-017-4911-9
 * HEPData entry: https://www.hepdata.net/record/ins1502620?table=Table%2011
 *
 * Note:
 * It is believed that the second kinematic variable, i.e. the elements of fKin2, for the Z data were
 * originally set to incorrect values for these two data sets. This should not have led to any incorrect
 * results and this has now been fixed. See https://github.com/NNPDF/buildmaster/pull/154 for more details.

 * Paper Table 13 (top)
 * Differential cross section for the central Z/gamma -> ll processes in the invariant mass region 
 * 46<mll<66 GeV combined from l=e,mu decays extrapolated to the common fiducial region. 
 */

#include "ATLAS.h"

// Central selection
void ATLAS_ZRAPPC_7TEVFilter::ReadData()
{

  cout << "********** WARNING: Converting pb to fb to match ApplGrid output ********" << endl;

  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << "ATLAS_ZRAPPC_7TEV/zypeak.dat";

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
 
  const int nBins = 1;
  const int ndataWZ[nBins] = {12};  // Data thresholds for Z_peak
  const double MWZ2[nBins]= {pow(91.0,2.0)};   //Mass squared of Z_peak

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

