/*
NB: Not currently suitable for use in NNPDF4.0 because no systematics breakdown available.

Differential cross section measurements of the single top and single antitop quark in the t-channel @LHC ATLAS 8 TeV

LHC-ATLAS 8 TeV
---------------

Selected events contain exactly one electron or muon, exactly two jets (exactly one of which must be b-tagged), and E_T^{miss} > 30 GeV.
Archived as: https://arxiv.org/pdf/1702.02859v3.pdf
Published in: Eur. Phys. J. C 77 (2017) 531

Distributions are converted, where necessary, so that they have the following dimensions:
Absolute transverse momentum: pb/GeV
Absolute rapidity: pb
Normalised transverse momentum: 1/GeV
Normalised rapidity: -
*/

#include "ATLAS_SINGLETOP_TCH_DIFF_8TEV.h"

// A - UNNORMALISED distributions

// 1) Normalised distribution differential in modulus of top quark rapidity

void ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_RAP_NORMFilter::ReadData()
{
  // Create stream to read data file
  fstream f1;

  // Data files
  stringstream datafile("");
  string filename1;
  filename1 = "ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_RAP";
  datafile << dataPath()
           << "rawdata/" << filename1 << "/" << filename1 << "_NORM.data";
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
    double sys1, sys2; // Systematic uncertainties from file
    double sys1_mult, sys2_mult; // Multiplicative systematic uncertainties
    double up, down, sigma, datshift; // Arguments of symmetriseErrors

    getline(f1,line);
    istringstream lstream(line);

    lstream >> rap_top_low >> rap_top_high;
    lstream >> fData[i] >> fStat[i];
    lstream >> sys1 >> sys2;

    fData[i] /= 1000;
    fStat[i] /= 1000;
    sys1 /= 1000;
    sys2 /= 1000;

    rap_top = 0.5*(rap_top_low + rap_top_high);
    
    fKin1[i] = rap_top;
    fKin2[i] = Mt*Mt; // Top mass squared
    fKin3[i] = 8000; // Centre of mass energy in GeV

    // Check whether numbers in file are percentage or absolute
    sys1_mult = (sys1/fData[i])*100;
    sys2_mult = (sys2/fData[i])*100;
    if (sys1_mult < 0)
    {
      up=sys2_mult;
      down=sys1_mult;
    }
    else
    {
      up=sys1_mult;
      down=sys2_mult;
    }
    symmetriseErrors(up, down, &sigma, &datshift);
    fSys[i][0].mult = sigma;
    fSys[i][0].add = fSys[i][0].mult*fData[i]/100;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CORR";
 
    fData[i] *= (1.0 + datshift*0.01); // Shift of central value due to asymmetric errors
  }
}

//==============================================================================

// 2) Normalised distribution differential in modulus of antitop quark rapidity

void ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_RAP_NORMFilter::ReadData()
{
  // Create stream to read data file
  fstream f1;

  // Data files
  stringstream datafile("");
  string filename1;
  filename1 = "ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_RAP";
  datafile << dataPath()
           << "rawdata/" << filename1 << "/" << filename1 << "_NORM.data";
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
    double sys1, sys2; // Systematic uncertainties from file
    double sys1_mult, sys2_mult; // Multiplicative systematic uncertainties
    double up, down, sigma, datshift; // Arguments of symmetriseErrors

    getline(f1,line);
    istringstream lstream(line);

    lstream >> rap_top_low >> rap_top_high;
    lstream >> fData[i] >> fStat[i];
    lstream >> sys1 >> sys2;

    fData[i] /= 1000;
    fStat[i] /= 1000;
    sys1 /= 1000;
    sys2 /= 1000;

    rap_top = 0.5*(rap_top_low + rap_top_high);

    fKin1[i] = rap_top;
    fKin2[i] = Mt*Mt; // Top mass squared
    fKin3[i] = 8000; // Centre of mass energy in GeV

    // Check whether numbers in file are percentage or absolute
    sys1_mult = (sys1/fData[i])*100;
    sys2_mult = (sys2/fData[i])*100;
    if (sys1_mult < 0)
    {
      up=sys2_mult;
      down=sys1_mult;
    }
    else
    {
      up=sys1_mult;
      down=sys2_mult;
    }
    symmetriseErrors(up, down, &sigma, &datshift);
    fSys[i][0].mult = sigma;
    fSys[i][0].add = fSys[i][0].mult*fData[i]/100;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CORR";
 
    fData[i] *= (1.0 + datshift*0.01); // Shift of central value due to asymmetric errors
  }
}

//==========================================================================

// 3) Normalised distribution differential in top quark transverse momentum

void ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_PT_NORMFilter::ReadData()
{
  // Create stream to read data file
  fstream f1;

  // Data files
  stringstream datafile("");
  string filename1;
  filename1 = "ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_PT";
  datafile << dataPath()
           << "rawdata/" << filename1 << "/" << filename1 << "_NORM.data";
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
    double pt_top; // Transverse momentum of top quark
    double pt_top_low, pt_top_high; // Limits of bin
    double sys1, sys2; // Systematic uncertainties from file
    double sys1_mult, sys2_mult; // Multiplicative systematic uncertainties
    double up, down, sigma, datshift; // Arguments of symmetriseErrors

    getline(f1,line);
    istringstream lstream(line);

    lstream >> pt_top_low >> pt_top_high;
    lstream >> fData[i] >> fStat[i];
    lstream >> sys1 >> sys2;

    // Convert to 1/GeV
    fData[i] /= 1000;
    fStat[i] /= 1000;
    sys1 /= 1000;
    sys2 /= 1000;

    pt_top = 0.5*(pt_top_low + pt_top_high);

    fKin1[i] = pt_top;
    fKin2[i] = Mt*Mt; // Top mass squared
    fKin3[i] = 8000; // Centre of mass energy in GeV

    // Check whether numbers in file are percentage or absolute
    sys1_mult = (sys1/fData[i])*100;
    sys2_mult = (sys2/fData[i])*100;
    if (sys1_mult < 0)
    {
      up=sys2_mult;
      down=sys1_mult;
    }
    else
    {
      up=sys1_mult;
      down=sys2_mult;
    }
    symmetriseErrors(up, down, &sigma, &datshift);
    fSys[i][0].mult = sigma;
    fSys[i][0].add = fSys[i][0].mult*fData[i]/100;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CORR";
 
    fData[i] *= (1.0 + datshift*0.01); // Shift of central value due to asymmetric errors
  }
}

//==============================================================================

// 4) Normalised distribution differential in antitop quark transverse momentum

void ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_PT_NORMFilter::ReadData()
{
  // Create stream to read data file
  fstream f1;

  // Data files
  stringstream datafile("");
  string filename1;
  filename1 = "ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_PT";
  datafile << dataPath()
           << "rawdata/" << filename1 << "/" << filename1 << "_NORM.data";
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
    double pt_top; // Transverse momentum of top quark
    double pt_top_low, pt_top_high; // Limits of bin
    double sys1, sys2; // Systematic uncertainties from file
    double sys1_mult, sys2_mult; // Multiplicative systematic uncertainties
    double up, down, sigma, datshift; // Arguments of symmetriseErrors

    getline(f1,line);
    istringstream lstream(line);

    lstream >> pt_top_low >> pt_top_high;
    lstream >> fData[i] >> fStat[i];
    lstream >> sys1 >> sys2;

    // Convert to 1/GeV
    fData[i] /= 1000;
    fStat[i] /= 1000;
    sys1 /= 1000;
    sys2 /= 1000;

    pt_top = 0.5*(pt_top_low + pt_top_high);

    fKin1[i] = pt_top;
    fKin2[i] = Mt*Mt; // Top mass squared
    fKin3[i] = 8000; // Centre of mass energy in GeV

    // Check whether numbers in file are percentage or absolute
    sys1_mult = (sys1/fData[i])*100;
    sys2_mult = (sys2/fData[i])*100;
    if (sys1_mult < 0)
    {
      up=sys2_mult;
      down=sys1_mult;
    }
    else
    {
      up=sys1_mult;
      down=sys2_mult;
    }
    symmetriseErrors(up, down, &sigma, &datshift);
    fSys[i][0].mult = sigma;
    fSys[i][0].add = fSys[i][0].mult*fData[i]/100;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CORR";
 
    fData[i] *= (1.0 + datshift*0.01); // Shift of central value due to asymmetric errors
  }
}

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
    double sys1, sys2; // Systematic uncertainties from file
    double sys1_mult, sys2_mult; // Multiplicative systematic uncertainties
    double up, down, sigma, datshift; // Arguments of symmetriseErrors

    getline(f1,line);
    istringstream lstream(line);

    lstream >> rap_top_low >> rap_top_high;
    lstream >> fData[i] >> fStat[i];
    lstream >> sys1 >> sys2;

    rap_top = 0.5*(rap_top_low + rap_top_high);

    fKin1[i] = rap_top;
    fKin2[i] = Mt*Mt; // Top mass squared
    fKin3[i] = 8000; // Centre of mass energy in GeV

    // Check whether numbers in file are percentage or absolute
    sys1_mult = (sys1/fData[i])*100;
    sys2_mult = (sys2/fData[i])*100;
    if (sys1_mult < 0)
    {
      up=sys2_mult;
      down=sys1_mult;
    }
    else
    {
      up=sys1_mult;
      down=sys2_mult;
    }
    symmetriseErrors(up, down, &sigma, &datshift);
    fSys[i][0].mult = sigma;
    fSys[i][0].add = fSys[i][0].mult*fData[i]/100;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CORR";
 
    fData[i] *= (1.0 + datshift*0.01); // Shift of central value due to asymmetric errors
  }
}

//==================================================================

// 6) Distribution differential in modulus of antitop quark rapidity

void ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_RAPFilter::ReadData()
{
  // Create stream to read data file
  fstream f1;

  // Data files
  stringstream datafile("");
  string filename1;
  filename1 = "ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_RAP";
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
    double sys1, sys2; // Systematic uncertainties from file
    double sys1_mult, sys2_mult; // Multiplicative systematic uncertainties
    double up, down, sigma, datshift; // Arguments of symmetriseErrors

    getline(f1,line);
    istringstream lstream(line);

    lstream >> rap_top_low >> rap_top_high;
    lstream >> fData[i] >> fStat[i];
    lstream >> sys1 >> sys2;

    rap_top = 0.5*(rap_top_low + rap_top_high);

    fKin1[i] = rap_top;
    fKin2[i] = Mt*Mt; // Top mass squared
    fKin3[i] = 8000; // Centre of mass energy in GeV

    // Check whether numbers in file are percentage or absolute
    sys1_mult = (sys1/fData[i])*100;
    sys2_mult = (sys2/fData[i])*100;
    if (sys1_mult < 0)
    {
      up=sys2_mult;
      down=sys1_mult;
    }
    else
    {
      up=sys1_mult;
      down=sys2_mult;
    }
    symmetriseErrors(up, down, &sigma, &datshift);
    fSys[i][0].mult = sigma;
    fSys[i][0].add = fSys[i][0].mult*fData[i]/100;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CORR";
 
    fData[i] *= (1.0 + datshift*0.01); // Shift of central value due to asymmetric errors
  }
}

//==================================================================

// 7) Distribution differential in top quark transverse momentum

void ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_PTFilter::ReadData()
{
  // Create stream to read data file
  fstream f1;

  // Data files
  stringstream datafile("");
  string filename1;
  filename1 = "ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_PT";
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
    double pt_top; // Transverse momentum of top quark
    double pt_top_low, pt_top_high; // Limits of bin
    double sys1, sys2; // Systematic uncertainties from file
    double sys1_mult, sys2_mult; // Multiplicative systematic uncertainties
    double up, down, sigma, datshift; // Arguments of symmetriseErrors

    getline(f1,line);
    istringstream lstream(line);

    lstream >> pt_top_low >> pt_top_high;
    lstream >> fData[i] >> fStat[i];
    lstream >> sys1 >> sys2;

    // Convert to pb/GeV
    fData[i] /= 1000;
    fStat[i] /= 1000;
    sys1 /= 1000;
    sys2 /= 1000;

    pt_top = 0.5*(pt_top_low + pt_top_high);

    fKin1[i] = pt_top;
    fKin2[i] = Mt*Mt; // Top mass squared
    fKin3[i] = 8000; // Centre of mass energy in GeV

    // Check whether numbers in file are percentage or absolute
    sys1_mult = (sys1/fData[i])*100;
    sys2_mult = (sys2/fData[i])*100;
    if (sys1_mult < 0)
    {
      up=sys2_mult;
      down=sys1_mult;
    }
    else
    {
      up=sys1_mult;
      down=sys2_mult;
    }
    symmetriseErrors(up, down, &sigma, &datshift);
    fSys[i][0].mult = sigma;
    fSys[i][0].add = fSys[i][0].mult*fData[i]/100;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CORR";
 
    fData[i] *= (1.0 + datshift*0.01); // Shift of central value due to asymmetric errors
  }
}

//==================================================================

// 8) Distribution differential in antitop quark transverse momentum

void ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_PTFilter::ReadData()
{
  // Create stream to read data file
  fstream f1;

  // Data files
  stringstream datafile("");
  string filename1;
  filename1 = "ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_PT";
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
    double pt_top; // Transverse momentum of top quark
    double pt_top_low, pt_top_high; // Limits of bin
    double sys1, sys2; // Systematic uncertainties from file
    double sys1_mult, sys2_mult; // Multiplicative systematic uncertainties
    double up, down, sigma, datshift; // Arguments of symmetriseErrors

    getline(f1,line);
    istringstream lstream(line);

    lstream >> pt_top_low >> pt_top_high;
    lstream >> fData[i] >> fStat[i];
    lstream >> sys1 >> sys2;

    // Convert to pb/GeV
    fData[i] /= 1000;
    fStat[i] /= 1000;
    sys1 /= 1000;
    sys2 /= 1000;

    pt_top = 0.5*(pt_top_low + pt_top_high);

    fKin1[i] = pt_top;
    fKin2[i] = Mt*Mt; // Top mass squared
    fKin3[i] = 8000; // Centre of mass energy in GeV

    // Check whether numbers in file are percentage or absolute
    sys1_mult = (sys1/fData[i])*100;
    sys2_mult = (sys2/fData[i])*100;
    if (sys1_mult < 0)
    {
      up=sys2_mult;
      down=sys1_mult;
    }
    else
    {
      up=sys1_mult;
      down=sys2_mult;
    }
    symmetriseErrors(up, down, &sigma, &datshift);
    fSys[i][0].mult = sigma;
    fSys[i][0].add = fSys[i][0].mult*fData[i]/100;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CORR";
 
    fData[i] *= (1.0 + datshift*0.01); // Shift of central value due to asymmetric errors
  }
}
