/*
 * NB: Not currently suitable for use in NNPDF4.0 because no systematics breakdown available.
 */

#include "ATLAS_SINGLETOP_TCH_DIFF_8TEV.h"

//==================================================================

// B - UNNORMALISED distributions

// 5) Distribution differential in modulus of top quark rapidity

void ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_RAPFilter::ReadData()
{
  // Create stream to read data file
  fstream f1;

  // Data files
  stringstream datafile("");
  string filename1;
  filename1 = "ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_RAP";
  datafile << dataPath()
           << "rawdata/" << filename1 << "/" << filename1 << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Start filter of data
  string line;

  // Skip over first two lines
  for (int i=0; i<2; i++)
  {
    getline(f1,line);
  }  

  for (int i=0; i<fNData; i++)
  {
    double rap_top; // Rapidity of top quark
    double rap_top_low, rap_top_high; // Limits of bin
    double stat; // Absolute statistical uncertainty
    double sys1, sys2; // Systematic uncertainties from file
    double sys1_mult, sys2_mult; // Multiplicative systematic uncertainties
    double up, down, sigma, datshift; // Arguments of symmetriseErrors

    getline(f1,line);
    istringstream lstream(line);

    lstream >> rap_top_low >> rap_top_high;
    lstream >> fData[i] >> stat;
    lstream >> sys1 >> sys2;

    rap_top = 0.5*(rap_top_low + rap_top_high);
    fStat[i] = (stat/fData[i])*100; // Store statistical uncertainty as percentage value

    fKin1[i] = rap_top;
    fKin2[i] = Mt*Mt; // Top mass squared
    fKin3[i] = 8000; // Centre of mass energy in GeV

    // Check whether numbers in file are percentage or absolute
    sys1_mult = (sys1/fData[i])*100;
    sys2_mult = (sys2/fData[i])*100;
    if (sys1_mult < 0) {up=sys2_mult; down=sys1_mult;}
    else {up=sys1_mult; down=sys2_mult;}
    symmetriseErrors(up, down, &sigma, &datshift);
    fSys[i][0].mult = sigma;
    fSys[i][0].add = fSys[i][0].mult*fData[i]/100;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CORR";
 
    fData[i] *= (1.0 + datshift*0.01); // Shift of central value due to asymmetric errors
  }
}
