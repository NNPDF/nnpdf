/*
 *
 * ATLASPHT11 - ATLAS inclusive photons 2011
 *
 * ATLAS isolated photon production, LHC 7 TeV, 4.6 fb^-1
 * Reference: Phys. Rev. D 89, 052004 (2014), [arXiv:1311.1440] 
 * HepData: http://hepdata.cedar.ac.uk/view/ins1300647
 *
 */

#include "ATLASPHT11.h"

/*
 * ATLASPHT11ETCTR dataset
 * 
 * d sigma/dE_{T}^{gamma} [pb / GeV], |eta^{gamma}| < 1.37
 *
 * Correlation information:
 * 
 * EnergyScale+/- is to be treated fully correlated between bins in one 
 * region of pseudorapidity (|eta|<1.37 or 1.52<|eta|<2.37), but uncorrelated 
 * between different regions of pseudorapidity.
 *
 * AlternativeMC and EnergyResolution are to be treated fully uncorrelated 
 * between bins.
 *
 * All other sources are to be treated fully correlated between all bins.
 *
 */

void ATLASPHT11ETGCTRFilter::ReadData()
{
  // Opening files
  fstream f1;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/ATLASPHT11/ATLASPHT11ETGCTR.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double mbin[fNData+1], shift[fNData], stmp, dtmp;

  // Reading data

  string line;
  
  // Filtering data
  for (int idat = 0; idat < fNData; idat++) {

    getline(f1, line);

    double fIsolationCutP, fIsolationCutM;
    double fEnergyScaleP, fEnergyScaleM;
    double fAlternativeMC;
    double fBackgroundMethod, fBackgroundIsolation;
    double fEnergyResolution;
    double fFragmentationP, fFragmentationM;
    double fPhotonEfficiencyP, fPhotonEfficiencyM;
    double fATLAS2011Luminosity;
    
    istringstream lstream(line);
    lstream >> mbin[idat] >> mbin[idat+1] >> fData[idat] >> fStat[idat]
	    >> fIsolationCutP >> fIsolationCutM
	    >> fEnergyScaleP >> fEnergyScaleM
	    >> fAlternativeMC
	    >> fBackgroundMethod >> fBackgroundIsolation
	    >> fEnergyResolution
	    >> fFragmentationP >> fFragmentationM
	    >> fPhotonEfficiencyP >> fPhotonEfficiencyM
	    >> fATLAS2011Luminosity; 

    //    cout << mbin[idat] << "   " << mbin[idat+1] << endl;
    
    shift[idat] = 0;

    // Systematics are in %

    // Isolation Cut: asymmetric, fully correlated between all bins
    symmetriseErrors(fIsolationCutP,fIsolationCutM,&stmp,&dtmp);
    shift[idat]+=dtmp;
    fSys[idat][0].mult=stmp;
    fSys[idat][0].type = MULT;
    fSys[idat][0].name = "CORR";

    // Energy Scale: asymmetric, fully correlated between bins in one region
    //               of pseudorapidity (|eta|<1.37 or 1.52<|eta|<2.37), but
    //               uncorrelated between different regions of pseudorapidity
    symmetriseErrors(fEnergyScaleP,fEnergyScaleM,&stmp,&dtmp);
    shift[idat]+=dtmp;
    fSys[idat][1].mult=stmp;
    fSys[idat][1].type = MULT;
    fSys[idat][1].name = "CORR";

    // Alternative MC: symmetric, fully uncorrelated between bins
    fSys[idat][2].mult=fAlternativeMC;
    fSys[idat][2].type = MULT;
    fSys[idat][2].name = "UNCORR";

    // Background Method: symmetric, fully correlated between all bins
    fSys[idat][3].mult=fBackgroundMethod;
    fSys[idat][3].type = MULT;
    fSys[idat][3].name = "CORR";

    // Background Isolation: symmetric, fully correlated between all bins
    fSys[idat][4].mult=fBackgroundIsolation;
    fSys[idat][4].type = MULT;
    fSys[idat][4].name = "CORR";

    // Energy Resolution: symmetric, fully uncorrelated between bins
    fSys[idat][5].mult=fEnergyResolution;
    fSys[idat][5].type = MULT;
    fSys[idat][5].name = "UNCORR";

    // Fragmentation: asymmetric, fully correlated between all bins
    symmetriseErrors(fFragmentationP,fFragmentationM,&stmp,&dtmp);
    shift[idat]+=dtmp;
    fSys[idat][6].mult=stmp;
    fSys[idat][6].type = MULT;
    fSys[idat][6].name = "CORR";

    // Photon Efficiency: asymmetric, fully correlated between all bins
    symmetriseErrors(fPhotonEfficiencyP,fPhotonEfficiencyM,&stmp,&dtmp);
    shift[idat]+=dtmp;
    fSys[idat][7].mult=stmp;
    fSys[idat][7].type = MULT;
    fSys[idat][7].name = "CORR";

    // ATLAS2011 Luminosity: symmetric, fully correlated between all bins
    fSys[idat][8].mult=fATLAS2011Luminosity;
    fSys[idat][8].type = MULT;
    fSys[idat][8].name = "ATLASLUMI11";

    // Shift of central values due to asymmetric uncertainties
    //fData[idat]+=shift*fData[idat]*0.01;

    // Compute absolute systematics for additive part
    for (int l = 0; l < fNSys; l++)
      fSys[idat][l].add = fSys[idat][l].mult*fData[idat]*1e-2;    

    // Kinematic variables
    
    fKin1[idat] = 1.37/2.;                            // Avg. eta_gamma (|eta_g|<1.37)
    fKin2[idat] = pow((mbin[idat] + mbin[idat+1]) * 0.5,2);  // Avg. Et of each bin
    fKin3[idat] = 7000.;                              // LHC 7 TeV
    
  }
  
  f1.close();

}



/*
 * ATLASPHT11ETFWD dataset
 * 
 * d sigma/dE_{T}^{gamma} [pb / GeV], 1.52 < |eta^{gamma}| < 2.37
 *
 * Correlation information:
 * 
 * EnergyScale+/- is to be treated fully correlated between bins in one 
 * region of pseudorapidity (|eta|<1.37 or 1.52<|eta|<2.37), but uncorrelated 
 * between different regions of pseudorapidity.
 *
 * AlternativeMC and EnergyResolution are to be treated fully uncorrelated 
 * between bins.
 *
 * All other sources are to be treated fully correlated between all bins.
 *
 */

void ATLASPHT11ETGFWDFilter::ReadData()
{
  // Opening files
  fstream f1;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/ATLASPHT11/ATLASPHT11ETGFWD.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double mbin[fNData+1], shift[fNData], stmp, dtmp;

  // Reading data

  string line;
  
  // Filtering data
  for (int idat = 0; idat < fNData; idat++) {

    getline(f1, line);

    double fIsolationCutP, fIsolationCutM;
    double fEnergyScaleP, fEnergyScaleM;
    double fAlternativeMC;
    double fBackgroundMethod, fBackgroundIsolation;
    double fEnergyResolution;
    double fFragmentationP, fFragmentationM;
    double fPhotonEfficiencyP, fPhotonEfficiencyM;
    double fATLAS2011Luminosity;
    
    istringstream lstream(line);
    lstream >> mbin[idat] >> mbin[idat+1] >> fData[idat] >> fStat[idat]
	    >> fIsolationCutP >> fIsolationCutM
	    >> fEnergyScaleP >> fEnergyScaleM
	    >> fAlternativeMC
	    >> fBackgroundMethod >> fBackgroundIsolation
	    >> fEnergyResolution
	    >> fFragmentationP >> fFragmentationM
	    >> fPhotonEfficiencyP >> fPhotonEfficiencyM
	    >> fATLAS2011Luminosity; 
    
    shift[idat] = 0;

    // Systematics are in %

    // Isolation Cut: asymmetric, fully correlated between all bins
    symmetriseErrors(fIsolationCutP,fIsolationCutM,&stmp,&dtmp);
    shift[idat]+=dtmp;
    fSys[idat][0].mult=stmp;
    fSys[idat][0].type = MULT;
    fSys[idat][0].name = "CORR";

    // Energy Scale: asymmetric, fully correlated between bins in one region
    //               of pseudorapidity (|eta|<1.37 or 1.52<|eta|<2.37), but
    //               uncorrelated between different regions of pseudorapidity
    symmetriseErrors(fEnergyScaleP,fEnergyScaleM,&stmp,&dtmp);
    shift[idat]+=dtmp;
    fSys[idat][1].mult=stmp;
    fSys[idat][1].type = MULT;
    fSys[idat][1].name = "CORR";

    // Alternative MC: symmetric, fully uncorrelated between bins
    fSys[idat][2].mult=fAlternativeMC;
    fSys[idat][2].type = MULT;
    fSys[idat][2].name = "UNCORR";

    // Background Method: symmetric, fully correlated between all bins
    fSys[idat][3].mult=fBackgroundMethod;
    fSys[idat][3].type = MULT;
    fSys[idat][3].name = "CORR";

    // Background Isolation: symmetric, fully correlated between all bins
    fSys[idat][4].mult=fBackgroundIsolation;
    fSys[idat][4].type = MULT;
    fSys[idat][4].name = "CORR";

    // Energy Resolution: symmetric, fully uncorrelated between bins
    fSys[idat][5].mult=fEnergyResolution;
    fSys[idat][5].type = MULT;
    fSys[idat][5].name = "UNCORR";

    // Fragmentation: asymmetric, fully correlated between all bins
    symmetriseErrors(fFragmentationP,fFragmentationM,&stmp,&dtmp);
    shift[idat]+=dtmp;
    fSys[idat][6].mult=stmp;
    fSys[idat][6].type = MULT;
    fSys[idat][6].name = "CORR";

    // Photon Efficiency: asymmetric, fully correlated between all bins
    symmetriseErrors(fPhotonEfficiencyP,fPhotonEfficiencyM,&stmp,&dtmp);
    shift[idat]+=dtmp;
    fSys[idat][7].mult=stmp;
    fSys[idat][7].type = MULT;
    fSys[idat][7].name = "CORR";

    // ATLAS2011 Luminosity: symmetric, fully correlated between all bins
    fSys[idat][8].mult=fATLAS2011Luminosity;
    fSys[idat][8].type = MULT;
    fSys[idat][8].name = "ATLASLUMI11";

    // Shift of central values due to asymmetric uncertainties
    //fData[idat]+=shift*fData[idat]*0.01;

    // Compute absolute systematics for additive part
    for (int l = 0; l < fNSys; l++)
      fSys[idat][l].add = fSys[idat][l].mult*fData[idat]*1e-2;    

    // Kinematic variables
    
    fKin1[idat] = (2.37-1.52)/2.;                     // Avg. eta_gamma (1.52<|eta_g|<2.37)
    fKin2[idat] = pow((mbin[idat] + mbin[idat+1]) * 0.5,2);  // Avg. Et of each bin
    fKin3[idat] = 7000.;                              // LHC 7 TeV
    
  }
  
  f1.close();

}


/*
 * ATLASPHT11ETAG dataset
 * 
 * d sigma/d eta^{gamma} [pb / GeV], E_T^g > 100 GeV
 *
 * Correlation information:
 * 
 * EnergyScale+/- is to be treated fully correlated between bins in one 
 * region of pseudorapidity (|eta|<1.37 or 1.52<|eta|<2.37), but uncorrelated 
 * between different regions of pseudorapidity.
 *
 * AlternativeMC and EnergyResolution are to be treated fully uncorrelated 
 * between bins.
 *
 * All other sources are to be treated fully correlated between all bins.
 *
 */

void ATLASPHT11ETAGFilter::ReadData()
{
  // Opening files
  fstream f1;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/ATLASPHT11/ATLASPHT11ETAG.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double mbin[fNData+1], shift[fNData], stmp, dtmp;

  // Reading data

  string line;
  
  // Filtering data
  for (int idat = 0; idat < fNData; idat++) {

    getline(f1, line);

    double fIsolationCutP, fIsolationCutM;
    double fEnergyScaleP, fEnergyScaleM;
    double fAlternativeMC;
    double fBackgroundMethod, fBackgroundIsolation;
    double fEnergyResolution;
    double fFragmentationP, fFragmentationM;
    double fPhotonEfficiencyP, fPhotonEfficiencyM;
    double fATLAS2011Luminosity;
    
    istringstream lstream(line);
    lstream >> mbin[idat] >> mbin[idat+1] >> fData[idat] >> fStat[idat]
	    >> fIsolationCutP >> fIsolationCutM
	    >> fEnergyScaleP >> fEnergyScaleM
	    >> fAlternativeMC
	    >> fBackgroundMethod >> fBackgroundIsolation
	    >> fEnergyResolution
	    >> fFragmentationP >> fFragmentationM
	    >> fPhotonEfficiencyP >> fPhotonEfficiencyM
	    >> fATLAS2011Luminosity; 
    
    shift[idat] = 0;

    // Systematics are in %

    // Isolation Cut: asymmetric, fully correlated between all bins
    symmetriseErrors(fIsolationCutP,fIsolationCutM,&stmp,&dtmp);
    shift[idat]+=dtmp;
    fSys[idat][0].mult=stmp;
    fSys[idat][0].type = MULT;
    fSys[idat][0].name = "CORR";

    // Energy Scale: asymmetric, fully correlated between bins in one region
    //               of pseudorapidity (|eta|<1.37 or 1.52<|eta|<2.37), but
    //               uncorrelated between different regions of pseudorapidity
    symmetriseErrors(fEnergyScaleP,fEnergyScaleM,&stmp,&dtmp);
    shift[idat]+=dtmp;
    fSys[idat][1].mult=stmp;
    fSys[idat][1].type = MULT;
    fSys[idat][1].name = "CORR";

    // Alternative MC: symmetric, fully uncorrelated between bins
    fSys[idat][2].mult=fAlternativeMC;
    fSys[idat][2].type = MULT;
    fSys[idat][2].name = "UNCORR";

    // Background Method: symmetric, fully correlated between all bins
    fSys[idat][3].mult=fBackgroundMethod;
    fSys[idat][3].type = MULT;
    fSys[idat][3].name = "CORR";

    // Background Isolation: symmetric, fully correlated between all bins
    fSys[idat][4].mult=fBackgroundIsolation;
    fSys[idat][4].type = MULT;
    fSys[idat][4].name = "CORR";

    // Energy Resolution: symmetric, fully uncorrelated between bins
    fSys[idat][5].mult=fEnergyResolution;
    fSys[idat][5].type = MULT;
    fSys[idat][5].name = "UNCORR";

    // Fragmentation: asymmetric, fully correlated between all bins
    symmetriseErrors(fFragmentationP,fFragmentationM,&stmp,&dtmp);
    shift[idat]+=dtmp;
    fSys[idat][6].mult=stmp;
    fSys[idat][6].type = MULT;
    fSys[idat][6].name = "CORR";

    // Photon Efficiency: asymmetric, fully correlated between all bins
    symmetriseErrors(fPhotonEfficiencyP,fPhotonEfficiencyM,&stmp,&dtmp);
    shift[idat]+=dtmp;
    fSys[idat][7].mult=stmp;
    fSys[idat][7].type = MULT;
    fSys[idat][7].name = "CORR";

    // ATLAS2011 Luminosity: symmetric, fully correlated between all bins
    fSys[idat][8].mult=fATLAS2011Luminosity;
    fSys[idat][8].type = MULT;
    fSys[idat][8].name = "ATLASLUMI11";

    // Shift of central values due to asymmetric uncertainties
    //fData[idat]+=shift*fData[idat]*0.01;

    // Compute absolute systematics for additive part
    for (int l = 0; l < fNSys; l++)
      fSys[idat][l].add = fSys[idat][l].mult*fData[idat]*1e-2;    

    // Kinematic variables
    
    fKin1[idat] = (mbin[idat] + mbin[idat+1]) * 0.5;  // Avg. eta of each bin
    fKin2[idat] = 100.*100.;                               // Min Et_gamma
    fKin3[idat] = 7000.;                              // LHC 7 TeV
    
  }
  
  f1.close();

}


