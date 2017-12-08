/*
 *
 * ATLASPHT12 - ATLAS inclusive photons 2012
 *
 * ATLAS isolated photon production, LHC 8 TeV, 20.2 fb^-1
 * Reference:  [arXiv:1605.03495] 
 * 
 *
 */

#include "ATLASPHT12.h"

/*
 * ATLASPHT12ETCTR dataset
 * 
 * d sigma/dE_{T}^{gamma} [pb / GeV], |eta^{gamma}| < 0.6
 *
 * Correlation information:
 * Systematics are already combined
 *
 */

void ATLASPHT12Filter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;


  stringstream datafileCNTR("");
  datafileCNTR << dataPath() << "rawdata/ATLASPHT12/ATLASPHT12ETGCTR.data";
  f1.open(datafileCNTR.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafileCNTR.str() << endl;
    exit(-1);
  }

  stringstream datafileFWD1("");
  datafileFWD1 << dataPath() << "rawdata/ATLASPHT12/ATLASPHT12ETGFWD1.data";
  f2.open(datafileFWD1.str().c_str(), ios::in);

  if (f2.fail()) {
    cerr << "Error opening data file " << datafileFWD1.str() << endl;
    exit(-1);
  }


  stringstream datafileFWD2("");
  datafileFWD2 << dataPath() << "rawdata/ATLASPHT12/ATLASPHT12ETGFWD2.data";
  f3.open(datafileFWD2.str().c_str(), ios::in);

  if (f3.fail()) {
    cerr << "Error opening data file " << datafileFWD2.str() << endl;
    exit(-1);
  }



  // Starting filter

  for (int idat = 0; idat < 18; idat++) {
    double upper, lower, stmp, dtmp;
    string line; 
    getline(f1,line);
    istringstream lstream(line);

    double fSystP, fSystM;
    double fATLAS2012Luminosity;
    
    lstream >> lower >> upper >> fData[idat] >> fStat[idat]  
    >> fSystP >> fSystM
    >> fATLAS2012Luminosity; 

    
    double shift = 0;


    // Convert to percentages

     fSystP = fSystP*100/fData[idat];
     fSystM = fSystM*100/fData[idat];

    // Total systematics
    symmetriseErrors(fSystP,fSystM,&stmp,&dtmp);

    // convert mult to a percentage
    fSys[idat][0].mult=stmp;
    fSys[idat][0].type = MULT;
    fSys[idat][0].name = "UNCORR";
    // convert back again
    fSys[idat][0].add = fSys[idat][0].mult*fData[idat]*1e-2;    

    shift += dtmp;

    // ATLAS2012 Luminosity: symmetric, fully correlated between all bins
    fSys[idat][1].mult=fATLAS2012Luminosity/fData[idat]*1e2;
    fSys[idat][1].type = MULT;
    fSys[idat][1].name = "ATLASLUMI12";
    fSys[idat][1].add = fSys[idat][1].mult*fData[idat]*1e-2;

    // Shift of central values due to asymmetric uncertainties
    // fData[idat]+=shift*fData[idat];

    fData[idat]*=(1.0 + shift*0.01); //Shift from asymmetric errors

    // Kinematic variables
    
    fKin1[idat] = 0.6/2.;                            // Avg. eta_gamma (|eta_g|<0.6)
    fKin2[idat] = pow((upper + lower) * 0.5,2);		// Avg. Et of each bin
    fKin3[idat] = 8000.;                              // LHC 8 TeV
    
  }
  

  // 2nd bin
  for (int idat = 18; idat < 35; idat++) {
    double upper, lower, stmp, dtmp;

    string line; 
    getline(f2,line);
    istringstream lstream(line);

    double fSystP, fSystM;
    double fATLAS2012Luminosity;
    
      lstream >> lower >> upper >> fData[idat] >> fStat[idat]  
      >> fSystP >> fSystM
      >> fATLAS2012Luminosity; 

    //    cout << mbin[idat] << "   " << mbin[idat+1] << endl;
    
      double shift = 0;

    // Convert to percentages

     fSystP = fSystP*100/fData[idat];
     fSystM = fSystM*100/fData[idat];

    // Total systematics
    symmetriseErrors(fSystP,fSystM,&stmp,&dtmp);

    // convert mult to a percentage
    fSys[idat][0].mult=stmp;
    fSys[idat][0].type = MULT;
    fSys[idat][0].name = "UNCORR";
    // convert back again
    fSys[idat][0].add = fSys[idat][0].mult*fData[idat]*1e-2;    

    shift += dtmp;

    // ATLAS2012 Luminosity: symmetric, fully correlated between all bins
    fSys[idat][1].mult=fATLAS2012Luminosity/fData[idat]*1e2;
    fSys[idat][1].type = MULT;
    fSys[idat][1].name = "ATLASLUMI12";
    fSys[idat][1].add = fSys[idat][1].mult*fData[idat]*1e-2;

    // Shift of central values due to asymmetric uncertainties
    // fData[idat]+=shift*fData[idat];

    fData[idat]*=(1.0 + shift*0.01);

    // Kinematic variables
    
    fKin1[idat] = (1.37-0.6)/2.;                     // Avg. eta_gamma (0.6<|eta_g|<1.37)
    fKin2[idat] = pow((upper + lower) * 0.5,2);		// Avg. Et of each bin
    fKin3[idat] = 8000.;                              // LHC 8 TeV
    
  }


  // 3rd bin
  for (int idat = 35; idat < fNData; idat++) {
    double upper, lower, stmp, dtmp;
 
    string line; 
    getline(f3,line);
    istringstream lstream(line);

    double fSystP, fSystM;
    double fATLAS2012Luminosity;
    
    lstream >> lower >> upper >> fData[idat] >> fStat[idat]  
    >> fSystP >> fSystM
    >> fATLAS2012Luminosity; 

    
      double shift = 0;

    // Convert to percentages

     fSystP = fSystP*100/fData[idat];
     fSystM = fSystM*100/fData[idat];

    // Total systematics
    symmetriseErrors(fSystP,fSystM,&stmp,&dtmp);

    // convert mult to a percentage
    fSys[idat][0].mult=stmp;
    fSys[idat][0].type = MULT;
    fSys[idat][0].name = "UNCORR";
    // convert back again
    fSys[idat][0].add = fSys[idat][0].mult*fData[idat]*1e-2;    

    shift += dtmp;

    // ATLAS2012 Luminosity: symmetric, fully correlated between all bins
    fSys[idat][1].mult=fATLAS2012Luminosity/fData[idat]*1e2;
    fSys[idat][1].type = MULT;
    fSys[idat][1].name = "ATLASLUMI12";
    fSys[idat][1].add = fSys[idat][1].mult*fData[idat]*1e-2;

    // Shift of central values due to asymmetric uncertainties
    // fData[idat]+=shift*fData[idat];

    fData[idat]*=(1.0 + shift*0.01);

    // Kinematic variables
    
    fKin1[idat] = (1.81-1.56)/2.;                     // Avg. eta_gamma (1.56<|eta_g|<1.81)
    fKin2[idat] = pow((upper + lower) * 0.5,2);		// Avg. Et of each bin
    fKin3[idat] = 8000.;                            // LHC 8 TeV
    
  }


  f1.close();
  f2.close();
  f3.close();
}





