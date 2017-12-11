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
  fstream cent, fwd1, fwd2;


  stringstream datafileCNTR("");
  datafileCNTR << dataPath() << "rawdata/ATLASPHT12/ATLASPHT12ETGCTR.data";
  cent.open(datafileCNTR.str().c_str(), ios::in);

  if (cent.fail()) {
    cerr << "Error opening data file " << datafileCNTR.str() << endl;
    exit(-1);
  }

  stringstream datafileFWD1("");
  datafileFWD1 << dataPath() << "rawdata/ATLASPHT12/ATLASPHT12ETGFWD1.data";
  fwd1.open(datafileFWD1.str().c_str(), ios::in);

  if (fwd1.fail()) {
    cerr << "Error opening data file " << datafileFWD1.str() << endl;
    exit(-1);
  }


  stringstream datafileFWD2("");
  datafileFWD2 << dataPath() << "rawdata/ATLASPHT12/ATLASPHT12ETGFWD2.data";
  fwd2.open(datafileFWD2.str().c_str(), ios::in);

  if (fwd2.fail()) {
    cerr << "Error opening data file " << datafileFWD2.str() << endl;
    exit(-1);
  }



  // Starting filter

  for (int idat = 0; idat < 18; idat++) {
    double upper, lower, stmp, datshift;
    string line; 
    getline(cent,line);
    istringstream lstream(line);

    double fSystP, fSystM;
    double fATLAS2012Luminosity;
    
    lstream >> lower >> upper >> fData[idat] >> fStat[idat]  
    >> fSystP >> fSystM
    >> fATLAS2012Luminosity; 

    // Convert to percentages

     fSystP = fSystP*100/fData[idat];
     fSystM = fSystM*100/fData[idat];

    // Total systematics
    symmetriseErrors(fSystP,fSystM,&stmp,&datshift);

    // convert mult to a percentage
    fSys[idat][0].mult=stmp;
    fSys[idat][0].type = MULT;
    fSys[idat][0].name = "UNCORR";
    // convert back again
    fSys[idat][0].add = fSys[idat][0].mult*fData[idat]*1e-2;    


    // ATLAS2012 Luminosity: symmetric, fully correlated between all bins
    fSys[idat][1].mult=fATLAS2012Luminosity/fData[idat]*1e2;
    fSys[idat][1].type = MULT;
    fSys[idat][1].name = "ATLASLUMI12";
    fSys[idat][1].add = fSys[idat][1].mult*fData[idat]*1e-2;

    // Shift of central values due to asymmetric uncertainties
    // fData[idat]+=shift*fData[idat];

    fData[idat]*=(1.0 + datshift*0.01); //Shift from asymmetric errors

    // Kinematic variables
    
    fKin1[idat] = 0.6/2.;                            // Avg. eta_gamma (|eta_g|<0.6)
    fKin2[idat] = pow((upper + lower) * 0.5,2);		// Avg. Et of each bin
    fKin3[idat] = 8000.;                              // LHC 8 TeV
    
  }
  

  // 2nd bin
  for (int idat = 18; idat < 35; idat++) {
    double upper, lower, stmp, datshift;

    string line; 
    getline(fwd1,line);
    istringstream lstream(line);

    double fSystP, fSystM;
    double fATLAS2012Luminosity;
    
      lstream >> lower >> upper >> fData[idat] >> fStat[idat]  
      >> fSystP >> fSystM
      >> fATLAS2012Luminosity; 

    //    cout << mbin[idat] << "   " << mbin[idat+1] << endl;
    

    // Convert to percentages

     fSystP = fSystP*100/fData[idat];
     fSystM = fSystM*100/fData[idat];

    // Total systematics
    symmetriseErrors(fSystP,fSystM,&stmp,&datshift);

    // convert mult to a percentage
    fSys[idat][0].mult=stmp;
    fSys[idat][0].type = MULT;
    fSys[idat][0].name = "UNCORR";
    // convert back again
    fSys[idat][0].add = fSys[idat][0].mult*fData[idat]*1e-2;    


    // ATLAS2012 Luminosity: symmetric, fully correlated between all bins
    fSys[idat][1].mult=fATLAS2012Luminosity/fData[idat]*1e2;
    fSys[idat][1].type = MULT;
    fSys[idat][1].name = "ATLASLUMI12";
    fSys[idat][1].add = fSys[idat][1].mult*fData[idat]*1e-2;

    // Shift of central values due to asymmetric uncertainties
    // fData[idat]+=shift*fData[idat];

    fData[idat]*=(1.0 + datshift*0.01);

    // Kinematic variables
    
    fKin1[idat] = (1.37-0.6)/2.;                     // Avg. eta_gamma (0.6<|eta_g|<1.37)
    fKin2[idat] = pow((upper + lower) * 0.5,2);		// Avg. Et of each bin
    fKin3[idat] = 8000.;                              // LHC 8 TeV
    
  }


  // 3rd bin
  for (int idat = 35; idat < fNData; idat++) {
    double upper, lower, stmp, datshift;
 
    string line; 
    getline(fwd2,line);
    istringstream lstream(line);

    double fSystP, fSystM;
    double fATLAS2012Luminosity;
    
    lstream >> lower >> upper >> fData[idat] >> fStat[idat]  
    >> fSystP >> fSystM
    >> fATLAS2012Luminosity; 

   
    // Convert to percentages

     fSystP = fSystP*100/fData[idat];
     fSystM = fSystM*100/fData[idat];

    // Total systematics
    symmetriseErrors(fSystP,fSystM,&stmp,&datshift);

    // convert mult to a percentage
    fSys[idat][0].mult=stmp;
    fSys[idat][0].type = MULT;
    fSys[idat][0].name = "UNCORR";
    // convert back again
    fSys[idat][0].add = fSys[idat][0].mult*fData[idat]*1e-2;    

    // ATLAS2012 Luminosity: symmetric, fully correlated between all bins
    fSys[idat][1].mult=fATLAS2012Luminosity/fData[idat]*1e2;
    fSys[idat][1].type = MULT;
    fSys[idat][1].name = "ATLASLUMI12";
    fSys[idat][1].add = fSys[idat][1].mult*fData[idat]*1e-2;

    // Shift of central values due to asymmetric uncertainties

    fData[idat]*=(1.0 + datshift*0.01);

    // Kinematic variables
    
    fKin1[idat] = (1.81-1.56)/2.;                     // Avg. eta_gamma (1.56<|eta_g|<1.81)
    fKin2[idat] = pow((upper + lower) * 0.5,2);		// Avg. Et of each bin
    fKin3[idat] = 8000.;                            // LHC 8 TeV
    
  }


  cent.close();
  fwd1.close();
  fwd2.close();
}





