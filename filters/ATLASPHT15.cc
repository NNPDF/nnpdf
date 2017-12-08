/*
 *
 * ATLASPHT15 - ATLAS inclusive photons 2015
 *
 * ATLAS isolated photon production, LHC 13 TeV, 3.2 fb^-1
 * Reference:  [arXiv:1701.06882] 
 * 
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
  fstream f1, f2, f3, f4;


  stringstream datafileCNTR("");
  datafileCNTR << dataPath() << "rawdata/ATLASPHT15/ATLASPHT15ETGCTR.data";
  f1.open(datafileCNTR.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafileCNTR.str() << endl;
    exit(-1);
  }

  stringstream datafileFWD1("");
  datafileFWD1 << dataPath() << "rawdata/ATLASPHT15/ATLASPHT15ETGFWD1.data";
  f2.open(datafileFWD1.str().c_str(), ios::in);

  if (f2.fail()) {
    cerr << "Error opening data file " << datafileFWD1.str() << endl;
    exit(-1);
  }


  stringstream datafileFWD2("");
  datafileFWD2 << dataPath() << "rawdata/ATLASPHT15/ATLASPHT15ETGFWD2.data";
  f3.open(datafileFWD2.str().c_str(), ios::in);

  if (f3.fail()) {
    cerr << "Error opening data file " << datafileFWD2.str() << endl;
    exit(-1);
  }
 stringstream datafileFWD3("");
  datafileFWD3 << dataPath() << "rawdata/ATLASPHT15/ATLASPHT15ETGFWD3.data";
  f4.open(datafileFWD3.str().c_str(), ios::in);

  if (f4.fail()) {
    cerr << "Error opening data file " << datafileFWD3.str() << endl;
    exit(-1);
  }


  // Starting filter

  // Reading data

  
  // Filtering data 1st bin
  for (int idat = 0; idat < 14; idat++) {
    double upper, lower, stmp, dtmp;
    string line; 
    getline(f1,line);
    istringstream lstream(line);

    double fSystP, fSystM;
    double fATLAS2015Luminosity;
    
	lstream >> lower >> upper >> fData[idat] >> fStat[idat]  
     >> fSystP >> fSystM
     >> fATLAS2015Luminosity; 

    
      double shift = 0;


    // Convert stat to absolute val
    fStat[idat] = fStat[idat]*fData[idat]*1e-2;

     // fSystP = fSystP*100/fData[idat];
     // fSystM = fSystM*100/fData[idat];

    // Total systematics
    symmetriseErrors(fSystP,fSystM,&stmp,&dtmp);

    fSys[idat][0].mult=stmp;
    fSys[idat][0].type = MULT;
    fSys[idat][0].name = "UNCORR";
    // convert to additive uncertainties
    fSys[idat][0].add = fSys[idat][0].mult*fData[idat]*1e-2;    

    shift += dtmp;

    // ATLAS2015 Luminosity: symmetric, fully correlated between all bins
    fSys[idat][1].mult=fATLAS2015Luminosity;
    fSys[idat][1].type = MULT;
    fSys[idat][1].name = "ATLASLUMI15";
    fSys[idat][1].add = fSys[idat][1].mult*fData[idat]*1e-2;

    // Shift of central values due to asymmetric uncertainties
    // fData[idat]+=shift*fData[idat];

    fData[idat]*=(1.0 + shift*0.01); //Shift from asymmetric errors

    // Kinematic variables
    
    fKin1[idat] = 0.6/2.;                            // Avg. eta_gamma (|eta_g|<0.6)
    fKin2[idat] = pow((upper + lower) * 0.5,2);		// Avg. Et of each bin
    fKin3[idat] = 13000.;                              // LHC 13 TeV
    
  }
  

  // 2nd bin
  for (int idat = 14; idat < 28; idat++) {
    double upper, lower, stmp, dtmp;

    string line; 
    getline(f2,line);
    istringstream lstream(line);

    double fSystP, fSystM;
    double fATLAS2015Luminosity;
    
	lstream >> lower >> upper >> fData[idat] >> fStat[idat]
	>> fSystP >> fSystM
    >> fATLAS2015Luminosity; 

    
      double shift = 0;

    // Convert to percentages
   //  Convert stat to absolute val
     
    fStat[idat] = fStat[idat]*fData[idat]*1e-2;


// fSystP = fSystP*100/fData[idat];
     // fSystM = fSystM*100/fData[idat];

    // Total systematics
       symmetriseErrors(fSystP,fSystM,&stmp,&dtmp);

    fSys[idat][0].mult=stmp;
    fSys[idat][0].type = MULT;
    fSys[idat][0].name = "UNCORR";
    // convert to additive uncertainties
    fSys[idat][0].add = fSys[idat][0].mult*fData[idat]*1e-2;    

    shift += dtmp;

    // ATLAS2015 Luminosity: symmetric, fully correlated between all bins
    fSys[idat][1].mult=fATLAS2015Luminosity;
    fSys[idat][1].type = MULT;
    fSys[idat][1].name = "ATLASLUMI15";
    fSys[idat][1].add = fSys[idat][1].mult*fData[idat]*1e-2;
    // Shift of central values due to asymmetric uncertainties
    // fData[idat]+=shift*fData[idat];

    fData[idat]*=(1.0 + shift*0.01);

    // Kinematic variables
    
    fKin1[idat] = (1.37-0.6)/2.;                     // Avg. eta_gamma (0.6<|eta_g|<1.37)
    fKin2[idat] = pow((upper + lower) * 0.5,2);		// Avg. Et of each bin
    fKin3[idat] = 13000.;                              // LHC 13 TeV
    
  }


  // 3rd bin
  for (int idat = 28; idat < 41; idat++) {
    double upper, lower, stmp, dtmp;
 
    string line; 
    getline(f3,line);
    istringstream lstream(line);

    double fSystP, fSystM;
    double fATLAS2015Luminosity;
    
	lstream >> lower >> upper >> fData[idat] >> fStat[idat]  
    >> fSystP >> fSystM
    >> fATLAS2015Luminosity; 

    
    double shift = 0;

    // Convert stat to absolute val
    fStat[idat] = fStat[idat]*fData[idat]*1e-2;


 	// Convert to percentages

     // fSystP = fSystP*100/fData[idat];
     // fSystM = fSystM*100/fData[idat];

    // Total systematics
    symmetriseErrors(fSystP,fSystM,&stmp,&dtmp);

    fSys[idat][0].mult=stmp;
    fSys[idat][0].type = MULT;
    fSys[idat][0].name = "UNCORR";
    // convert to additive uncertainties
    fSys[idat][0].add = fSys[idat][0].mult*fData[idat]*1e-2;    

    shift += dtmp;

    // ATLAS2015 Luminosity: symmetric, fully correlated between all bins
    fSys[idat][1].mult=fATLAS2015Luminosity;
    fSys[idat][1].type = MULT;
    fSys[idat][1].name = "ATLASLUMI15";
    fSys[idat][1].add = fSys[idat][1].mult*fData[idat]*1e-2;

    // Shift of central values due to asymmetric uncertainties
    // fData[idat]+=shift*fData[idat];

    fData[idat]*=(1.0 + shift*0.01);

    // Kinematic variables
    
    fKin1[idat] = (1.81-1.56)/2.;                     // Avg. eta_gamma (1.56<|eta_g|<1.81)
    fKin2[idat] = pow((upper + lower) * 0.5,2);		// Avg. Et of each bin
    fKin3[idat] = 13000.;                            // LHC 13 TeV
    
  }

 for (int idat = 41; idat < fNData; idat++) {
    double upper, lower, stmp, dtmp;
 
    string line; 
    getline(f4,line);
    istringstream lstream(line);

    double fSystP, fSystM;
    double fATLAS2015Luminosity;
    
	lstream >> lower >> upper >> fData[idat] >> fStat[idat]        
    >> fSystP >> fSystM
    >> fATLAS2015Luminosity; 

    
    double shift = 0;

    // Convert to percentages

    // Convert stat to absolute val

    fStat[idat] = fStat[idat]*fData[idat]*1e-2;


     // fSystP = fSystP*100/fData[idat];
     // fSystM = fSystM*100/fData[idat];

    // Total systematics
       symmetriseErrors(fSystP,fSystM,&stmp,&dtmp);

    fSys[idat][0].mult=stmp;
    fSys[idat][0].type = MULT;
    fSys[idat][0].name = "UNCORR";
    // convert to additive uncertainties
    fSys[idat][0].add = fSys[idat][0].mult*fData[idat]*1e-2;    

    shift += dtmp;

    // ATLAS2015 Luminosity: symmetric, fully correlated between all bins
    fSys[idat][1].mult=fATLAS2015Luminosity;
    fSys[idat][1].type = MULT;
    fSys[idat][1].name = "ATLASLUMI15";
    fSys[idat][1].add = fSys[idat][1].mult*fData[idat]*1e-2;

    // Shift of central values due to asymmetric uncertainties
    // fData[idat]+=shift*fData[idat];

    fData[idat]*=(1.0 + shift*0.01);

    // Kinematic variables
    
    fKin1[idat] = (2.37-1.81)/2.;                     // Avg. eta_gamma (1.56<|eta_g|<1.81)
    fKin2[idat] = pow((upper + lower) * 0.5,2);		// Avg. Et of each bin
    fKin3[idat] = 13000.;                            // LHC 13 TeV
    
  }



  f1.close();
  f2.close();
  f3.close();
  f4.close();
}





