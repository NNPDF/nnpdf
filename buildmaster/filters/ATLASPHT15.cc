/*
 *
 * ATLASPHT15 - ATLAS inclusive photons 2015
 *
 * ATLAS isolated photon production, LHC 13 TeV, 3.2 fb^-1
 * Reference:  [arXiv:1701.06882] 
 * Datafiles from http://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/STDM-2016-08/
 * Systematics are provided already combined - these are used as no covariance matrix is available
 *
 */

#include "ATLASPHT15.h"

/*
 * ATLASPHT15ETCTR dataset
 * 
 * d sigma/dE_{T}^{gamma} [pb / GeV], |eta^{gamma}| < 0.6
 *
 * Correlation information:
 * Systematics are already combined
 *
 */

void ATLASPHT15Filter::ReadData()
{
  // Opening files
  fstream cent, fwd1, fwd2, fwd3;


  stringstream datafileCNTR("");
  datafileCNTR << dataPath() << "rawdata/ATLASPHT15_SF/ATLASPHT15ETGCTR.data";
  cent.open(datafileCNTR.str().c_str(), ios::in);

  if (cent.fail()) {
    cerr << "Error opening data file " << datafileCNTR.str() << endl;
    exit(-1);
  }

  stringstream datafileFWD1("");
  datafileFWD1 << dataPath() << "rawdata/ATLASPHT15_SF/ATLASPHT15ETGFWD1.data";
  fwd1.open(datafileFWD1.str().c_str(), ios::in);

  if (fwd1.fail()) {
    cerr << "Error opening data file " << datafileFWD1.str() << endl;
    exit(-1);
  }


  stringstream datafileFWD2("");
  datafileFWD2 << dataPath() << "rawdata/ATLASPHT15_SF/ATLASPHT15ETGFWD2.data";
  fwd2.open(datafileFWD2.str().c_str(), ios::in);

  if (fwd2.fail()) {
    cerr << "Error opening data file " << datafileFWD2.str() << endl;
    exit(-1);
  }
 stringstream datafileFWD3("");
  datafileFWD3 << dataPath() << "rawdata/ATLASPHT15_SF/ATLASPHT15ETGFWD3.data";
  fwd3.open(datafileFWD3.str().c_str(), ios::in);

  if (fwd3.fail()) {
    cerr << "Error opening data file " << datafileFWD3.str() << endl;
    exit(-1);
  }


  // Starting filter

  // Reading data

  
  // Filtering data 1st bin
  for (int idat = 0; idat < 14; idat++) {
    double upper, lower, stmp, datshift;
    string line; 
    getline(cent,line);
    istringstream lstream(line);

    double fSystP, fSystM;
    double fATLAS2015Luminosity;
    
	lstream >> lower >> upper >> fData[idat] >> fStat[idat]  
     >> fSystP >> fSystM
     >> fATLAS2015Luminosity; 



    // Convert stat to absolute val
    fStat[idat] = fStat[idat]*fData[idat]*1e-2;

    // and systematics
    fSystP = fSystP*fData[idat]*1e-2;
    fSystM = fSystM*fData[idat]*1e-2;

    // Total systematics
    symmetriseErrors(fSystP,fSystM,&stmp,&datshift);

    fSys[idat][0].add=stmp;
    fSys[idat][0].type = MULT;
    fSys[idat][0].name = "UNCORR";
    // Shift of central values due to asymmetric uncertainties
    fData[idat] += datshift;

    // convert to mult uncertainties
    fSys[idat][0].mult = fSys[idat][0].add/fData[idat]*1e2;


    // ATLAS2015 Luminosity: symmetric, fully correlated between all bins
    fSys[idat][1].mult=fATLAS2015Luminosity;
    fSys[idat][1].type = MULT;
    fSys[idat][1].name = "ATLASLUMI15";
    fSys[idat][1].add = fSys[idat][1].mult*fData[idat]*1e-2;

    // Kinematic variables
    
    fKin1[idat] = 0.6/2.;                            // Avg. eta_gamma (|eta_g|<0.6)
    fKin2[idat] = pow((upper + lower) * 0.5,2);		// Avg. Et of each bin
    fKin3[idat] = 13000.;                              // LHC 13 TeV
    
  }
  

  // 2nd bin
  for (int idat = 14; idat < 28; idat++) {
    double upper, lower, stmp, datshift;

    string line; 
    getline(fwd1,line);
    istringstream lstream(line);

    double fSystP, fSystM;
    double fATLAS2015Luminosity;
    
	lstream >> lower >> upper >> fData[idat] >> fStat[idat]
	>> fSystP >> fSystM
    >> fATLAS2015Luminosity; 

    // Convert to percentages
    //  Convert stat to absolute val
     
    fStat[idat] = fStat[idat]*fData[idat]*1e-2;

    // and systematics
    fSystP = fSystP*fData[idat]*1e-2;
    fSystM = fSystM*fData[idat]*1e-2;

    // Total systematics
    symmetriseErrors(fSystP,fSystM,&stmp,&datshift);

    fSys[idat][0].add=stmp;
    fSys[idat][0].type = MULT;
    fSys[idat][0].name = "UNCORR";
    // Shift of central values due to asymmetric uncertainties
    fData[idat] += datshift;
    // convert to mult uncertainties
    fSys[idat][0].mult = fSys[idat][0].add/fData[idat]*1e2;

    // ATLAS2015 Luminosity: symmetric, fully correlated between all bins
    fSys[idat][1].mult=fATLAS2015Luminosity;
    fSys[idat][1].type = MULT;
    fSys[idat][1].name = "ATLASLUMI15";
    fSys[idat][1].add = fSys[idat][1].mult*fData[idat]*1e-2;

    // Kinematic variables
    
    fKin1[idat] = (1.37-0.6)/2.;                     // Avg. eta_gamma (0.6<|eta_g|<1.37)
    fKin2[idat] = pow((upper + lower) * 0.5,2);		// Avg. Et of each bin
    fKin3[idat] = 13000.;                              // LHC 13 TeV
    
  }


  // 3rd bin
  for (int idat = 28; idat < 41; idat++) {
    double upper, lower, stmp, datshift;
 
    string line; 
    getline(fwd2,line);
    istringstream lstream(line);

    double fSystP, fSystM;
    double fATLAS2015Luminosity;
    
	lstream >> lower >> upper >> fData[idat] >> fStat[idat]  
    >> fSystP >> fSystM
    >> fATLAS2015Luminosity; 



    // Convert stat to absolute val
    fStat[idat] = fStat[idat]*fData[idat]*1e-2;

    // and systematics
    fSystP = fSystP*fData[idat]*1e-2;
    fSystM = fSystM*fData[idat]*1e-2;

    // Total systematics
    symmetriseErrors(fSystP,fSystM,&stmp,&datshift);

    fSys[idat][0].add=stmp;
    fSys[idat][0].type = MULT;
    fSys[idat][0].name = "UNCORR";
    // Shift of central values due to asymmetric uncertainties
    fData[idat] += datshift;
    // convert to mult uncertainties
    fSys[idat][0].mult = fSys[idat][0].add/fData[idat]*1e2;

    // ATLAS2015 Luminosity: symmetric, fully correlated between all bins
    fSys[idat][1].mult=fATLAS2015Luminosity;
    fSys[idat][1].type = MULT;
    fSys[idat][1].name = "ATLASLUMI15";
    fSys[idat][1].add = fSys[idat][1].mult*fData[idat]*1e-2;

    // Kinematic variables
    
    fKin1[idat] = (1.81-1.56)/2.;                     // Avg. eta_gamma (1.56<|eta_g|<1.81)
    fKin2[idat] = pow((upper + lower) * 0.5,2);		// Avg. Et of each bin
    fKin3[idat] = 13000.;                            // LHC 13 TeV
    
  }

 for (int idat = 41; idat < fNData; idat++) {
    double upper, lower, stmp, datshift;
 
    string line; 
    getline(fwd3,line);
    istringstream lstream(line);

    double fSystP, fSystM;
    double fATLAS2015Luminosity;
    
	lstream >> lower >> upper >> fData[idat] >> fStat[idat]        
    >> fSystP >> fSystM
    >> fATLAS2015Luminosity; 


    // Convert to percentages

    // Convert stat to absolute val

    fStat[idat] = fStat[idat]*fData[idat]*1e-2;

    // and systematics
    fSystP = fSystP*fData[idat]*1e-2;
    fSystM = fSystM*fData[idat]*1e-2;

    // Total systematics
    symmetriseErrors(fSystP,fSystM,&stmp,&datshift);

    fSys[idat][0].add=stmp;
    fSys[idat][0].type = MULT;
    fSys[idat][0].name = "UNCORR";
    // Shift of central values due to asymmetric uncertainties
    fData[idat] += datshift;
    // convert to mult uncertainties
    fSys[idat][0].mult = fSys[idat][0].add/fData[idat]*1e2;

    // ATLAS2015 Luminosity: symmetric, fully correlated between all bins
    fSys[idat][1].mult=fATLAS2015Luminosity;
    fSys[idat][1].type = MULT;
    fSys[idat][1].name = "ATLASLUMI15";
    fSys[idat][1].add = fSys[idat][1].mult*fData[idat]*1e-2;

    // Kinematic variables
    
    fKin1[idat] = (2.37-1.81)/2.;                     // Avg. eta_gamma (1.56<|eta_g|<1.81)
    fKin2[idat] = pow((upper + lower) * 0.5,2);		// Avg. Et of each bin
    fKin3[idat] = 13000.;                            // LHC 13 TeV
    
  }



  cent.close();
  fwd1.close();
  fwd2.close();
  fwd3.close();
}





