/*
Differential cross section measurements of the single top and single antitop quark in the t-channel @LHC ATLAS 7 TeV

LHC-ATLAS 7 TeV
---------------

Selected events contain one charged lepton, large missing transverse momentum, and two or three jets (L = 4.59 1/fb)
Archived as: https://arxiv.org/pdf/1406.7844v2.pdf
Published in: Physics Review D 90, 112006 (https://journals.aps.org/prd/abstract/10.1103/PhysRevD.90.112006)

Eight distributions are implemented here. These are normalised and unnormalised distributions differential in:
1) Top quark absolute rapidity
2) Antitop quark absolute rapidity
3) Top quark transverse momentum
4) Antitop quark transverse momentum

Description of raw data:
Cross sections and percentage statistical uncertainties are taken from Tables VI and VII of the paper.
The breakdowns of systematic uncertainties are taken from Tables IX-XVI of the paper.
Statistical correlation matrices for calculating bin-wise correlations of the statistical uncertainties
are taken from Figures 17 and 18 in the paper. HOWEVER, note that 17a and 17d, i.e. the matrices for the
pT(t) and y(tbar) distributions, are not positive semidefinite so cannot be used in the code. The experementalists
have been told about this, but while we wait for corrected matrices, bin-wise correlations are not computed for
these two distributions. 

Distributions are converted, where necessary, so that they have the following dimensions:
Absolute transverse momentum: pb/GeV
Absolute rapidity: pb
Normalised transverse momentum: 1/GeV
Normalised rapidity: -

Note that the data files can be found in the supplemental material here:
https://journals.aps.org/prd/abstract/10.1103/PhysRevD.90.112006

Notes:
1) The number of systematic uncertainties considered in the code is distribution-dependent.
2) All systematics are assumed to be multiplicative.
3) All systematics are treated as CORR (i.e. correlated), except from the luminosity uncertainty for the
unnormalised distributions which are treated as ATLASLUMI11 (i.e. ATLAS luminosity for the 2011 data set).
*/

#include "ATLAS_SINGLETOP_TCH_DIFF_7TEV.h"
#include "NNPDF/utils.h"

// A - NORMALISED distributions

// 1) Distribution differential in modulus of top quark rapidity
void ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_RAP_NORMFilter::ReadData()
{
  // Create streams to read data files
  fstream f1, f2, f3;

  // Data files
  stringstream datafile("");
  string filename1;
  filename1 = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_RAP";
  datafile << dataPath()
           << "rawdata/" << filename1 << "/" << filename1 << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream sysfile("");
  string filename2;
  filename2 = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_RAP_NORM_SYS_BREAKDOWN";
  sysfile << dataPath()
          << "rawdata/" << filename1 << "/" << filename2 << ".data";
  f2.open(sysfile.str().c_str(), ios::in);

  if (f2.fail())
  {
    cerr << "Error opening data file " << sysfile.str() << endl;
    exit(-1);
  }

  // Open correlation matrix file
  stringstream corrfile("");
  string filename3;
  filename3 = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_RAP_NORM";
  corrfile << dataPath()
           << "rawdata/" << filename1 << "/" << filename3 << ".corr";
  f3.open(corrfile.str().c_str(), ios::in);

  if (f3.fail())
  {
    cerr << "Error opening data file " << corrfile.str() << endl;
    exit(-1);
  }

  // Start filter of data
  string line;

  // Initialise array to store additive stat. uncerts.
  std::vector<double> fstat_additive(fNData);

  // Skip over first ten lines
  for (int i=0; i<10; i++)
  {
    getline(f1,line);
  }
  
  for (int i=0; i<fNData; i++)
  {
    double rap_top; // Rapidity of top quark
    double rap_top_low, rap_top_high; // Limits of bin
    double fstat_percentage; // Percentage statistical uncertainty
    string unneeded_info;

    getline(f1,line);
    istringstream lstream(line);

    lstream >> rap_top_low >> unneeded_info >> rap_top_high;
    rap_top = 0.5*(rap_top_low + rap_top_high);

    // Skip over next eight elements of line
    for (int j=0; j<8; j++)
    {
      lstream >> unneeded_info;
    }

    fKin1[i] = rap_top;
    fKin2[i] = Mt*Mt; // Top mass squared
    fKin3[i] = 7000; // Centre of mass energy in GeV

    lstream >> fData[i]; // Value of bin
    lstream >> unneeded_info >> fstat_percentage;
    fstat_additive[i] = fstat_percentage*fData[i]/100;

    fStat[i] = 1e-10; // Set stat. error to zero to avoid double counting when using artificial systematics
  }

  // Read statistical correlation matrix
  // Skip over first ten lines
  for (int i=0; i<10; i++)
  {
    getline(f3,line);
  }

  double** covmat = new double*[fNData];
  NNPDF::matrix<double> corrmat(fNData, fNData);
  for (int i=0; i<fNData; i++)
  {
    string unneeded_info;
    covmat[i] = new double[fNData];
    getline(f3,line);
    istringstream lstream(line);
    lstream >> unneeded_info >> unneeded_info >> unneeded_info;
    for (int j=0; j<fNData; j++)
    {
      lstream >> corrmat(i,j) >> unneeded_info;
      covmat[i][j] = corrmat(i,j) * fstat_additive[i] * fstat_additive[j];
    }
  }

  // Generate artificial systematics
  double** syscor = new double*[fNData];
  for (int i=0; i<fNData; i++)
    syscor[i] = new double[fNData];

  if (!genArtSys(fNData,covmat,syscor))
  {
    cerr << " in " << fSetName << endl;
    exit(-1);
  }

  // Assign artificial systematics
  for (int i=0; i<fNData; i++)
  {
    for (int j=0; j<fNData; j++)
      {
        fSys[i][j].add = syscor[i][j];
        fSys[i][j].mult = fSys[i][j].add*100/fData[i];
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }
  }

  // Read file with systematic uncertainty breakdown
  // Skip over first 20 lines (including stat. uncert.)
  for (int i=0; i<20; i++)
  {
    getline(f2,line);
  }
  
  double sys1, sys2, up, down, sigma, datshift;
  std::vector<double> shift(fNData);

  for (int j=fNData; j<fNSys; j++)
  {
    string unneeded_info;

    getline(f2,line);
    istringstream lstream(line);

    if (j==fNData+2) // Deal with aysymmetric errors
    {
      for (int i=0; i<3; i++)
      {
        lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
        if (sys1 < 0) 
        {
          up=sys2;
          down=sys1;
        }
        else
        {
          up=sys1;
          down=sys2;
        }
        symmetriseErrors(up, down, &sigma, &datshift);
        fSys[i][j].mult = sigma;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";    
        shift[i] += datshift; 
      }
  
      lstream >> fSys[3][j].mult;
      fSys[3][j].add = fSys[3][j].mult*fData[3]/100;
      fSys[3][j].type = MULT;
      fSys[3][j].name = "CORR";
    }
    else if (j==fNData+3) // Deal with asymmetric errors
    {
      for (int i=0; i<2; i++)
      {
        lstream >> fSys[i][j].mult >> unneeded_info;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }

      lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
      if (sys1 < 0)
      {
        up=sys2;
        down=sys1;
      }
      else
      {
        up=sys1;
        down=sys2;
      }
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[2][j].mult = sigma;
      fSys[2][j].add = fSys[2][j].mult*fData[2]/100;
      fSys[2][j].type = MULT;
      fSys[2][j].name = "CORR";
      shift[2] += datshift;

      lstream >> fSys[3][j].mult;
      fSys[3][j].add = fSys[3][j].mult*fData[3]/100;
      fSys[3][j].type = MULT;
      fSys[3][j].name = "CORR";
    }
    else if (j==fNData+5) // Deal with asymmetric errors
    {
      for (int i=0; i<2; i++)
      {
        lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
        if (sys1 < 0)
        {
          up=sys2;
          down=sys1;
        }
        else
        {
          up=sys1;
          down=sys2;
        }
        symmetriseErrors(up, down, &sigma, &datshift);
        fSys[i][j].mult = sigma;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
        shift[i] += datshift;
      }

      lstream >> fSys[2][j].mult >> unneeded_info;
      fSys[2][j].add = fSys[2][j].mult*fData[2]/100;
      fSys[2][j].type = MULT;
      fSys[2][j].name = "CORR";

      lstream >> sys1 >> unneeded_info >> sys2;
      if (sys1 < 0)
      {
        up=sys2;
        down=sys1;
      }
      else
      {
        up=sys1;
        down=sys2;
      }
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[3][j].mult = sigma;
      fSys[3][j].add = fSys[3][j].mult*fData[3]/100;
      fSys[3][j].type = MULT;
      fSys[3][j].name = "CORR";
      shift[3] += datshift;
    }
    else // Deal with lines that contain no asymmetric errors
    {
      for (int i=0; i<fNData; i++)
      {
        lstream >> fSys[i][j].mult >> unneeded_info;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }
    }
  }

  for (int i=0; i<fNData; i++)
  {
    fData[i] *= (1.0 + shift[i]*0.01); // Shift of central value due to asymmetric errors
  }

  // Clean-up
  for (int i=0; i<fNData; i++)
    delete[] syscor[i];

  delete[] syscor;

  for(int i=0; i<fNData; i++)
    delete[] covmat[i];

  delete[] covmat;
}

//==================================================================

// 2) Distribution differential in modulus of antitop quark rapidity
void ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_RAP_NORMFilter::ReadData()
{
  // Create streams to read data files
  fstream f1, f2;

  // Data files
  stringstream datafile("");
  string filename1;
  filename1 = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_RAP";
  datafile << dataPath()
           << "rawdata/" << filename1 << "/" << filename1 << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream sysfile("");
  string filename2;
  filename2 = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_RAP_NORM_SYS_BREAKDOWN";
  sysfile << dataPath()
          << "rawdata/" << filename1 << "/" << filename2 << ".data";
  f2.open(sysfile.str().c_str(), ios::in);

  if (f2.fail())
  {
    cerr << "Error opening data file " << sysfile.str() << endl;
    exit(-1);
  }

  // Start filter of data
  string line;

  // Skip over first ten lines
  for (int i=0; i<10; i++)
  {
    getline(f1,line);
  }
  
  for (int i=0; i<fNData; i++)
  {
    double rap_top; // Rapidity of top quark
    double rap_top_high, rap_top_low; // Limits of bin
    string unneeded_info;

    getline(f1,line);
    istringstream lstream(line);

    lstream >> rap_top_low >> unneeded_info >> rap_top_high;
    rap_top = 0.5*(rap_top_low + rap_top_high);

    // Skip over next eight elements of line
    for (int j=0; j<8; j++)
    {
      lstream >> unneeded_info;
    }

    fKin1[i] = rap_top;
    fKin2[i] = Mt*Mt; // Top mass squared
    fKin3[i] = 7000; // Centre of mass energy in GeV

    lstream >> fData[i]; // Value of bin
    // We do not read in stat. correlation matrix in this filter because it is not positive semi-definite
    lstream >> unneeded_info >> fStat[i]; // Statistical (percentage) uncertainty
    fStat[i] = fStat[i]*fData[i]/100.;
  }

  // Read file with systematic uncertainty breakdown
  // Skip over first 20 lines (including stat. uncert.)
  for (int i=0; i<20; i++)
  {
    getline(f2,line);
  }
  
  double sys1, sys2, up, down, sigma, datshift;
  std::vector<double> shift(fNData);

  for (int j=0; j<fNSys; j++)
  {
    string unneeded_info;

    getline(f2,line);
    istringstream lstream(line);

    if (j==2) // Deal with aysymmetric errors
    {
      lstream >> fSys[0][j].mult >> unneeded_info;
      fSys[0][j].add = fSys[0][j].mult*fData[0]/100;
      fSys[0][j].type = MULT;
      fSys[0][j].name = "CORR";

      lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
      if (sys1 < 0)
      {
        up=sys2;
        down=sys1;
      }
      else
      {
        up=sys1;
        down=sys2;
      }
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[1][j].mult = sigma;
      fSys[1][j].add = fSys[1][j].mult*fData[1]/100;
      fSys[1][j].type = MULT;
      fSys[1][j].name = "CORR";    
      shift[1] += datshift; 
  
      lstream >> fSys[2][j].mult >> unneeded_info;
      fSys[2][j].add = fSys[2][j].mult*fData[2]/100;
      fSys[2][j].type = MULT;
      fSys[2][j].name = "CORR";

      lstream >> sys1 >> unneeded_info >> sys2;
      if (sys1 < 0)
      {
        up=sys2;
        down=sys1;
      }
      else
      {
        up=sys1;
        down=sys2;
      }
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[3][j].mult = sigma;
      fSys[3][j].add = fSys[3][j].mult*fData[3]/100;
      fSys[3][j].type = MULT;
      fSys[3][j].name = "CORR";    
      shift[3] += datshift;
    }
    else if (j==4) // Deal with asymmetric errors
    {
      for (int i=0; i<fNData; i++)
      {
        lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
        if (sys1 < 0)
        {
          up=sys2;
          down=sys1;
        }
        else
        {
          up=sys1;
          down=sys2;
        }
        symmetriseErrors(up, down, &sigma, &datshift);
        fSys[i][j].mult = sigma;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";    
        shift[i] += datshift; 
      }
    }
    else // Deal with lines that contain no asymmetric errors
    {
      for (int i=0; i<fNData; i++)
      {
        lstream >> fSys[i][j].mult >> unneeded_info;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }
    }
  }

  for (int i=0; i<fNData; i++)
  {
    fData[i] *= (1.0 + shift[i]*0.01); // Shift of central value due to asymmetric errors
  }
}

//==================================================================

// 3) Distribution differential in top quark transverse momentum
void ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_PT_NORMFilter::ReadData()
{
  // Create streams to read data files
  fstream f1, f2;

  // Data files
  stringstream datafile("");
  string filename1;
  filename1 = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_PT";
  datafile << dataPath()
           << "rawdata/" << filename1 << "/" << filename1 << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream sysfile("");
  string filename2;
  filename2 = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_PT_NORM_SYS_BREAKDOWN";
  sysfile << dataPath()
          << "rawdata/" << filename1 << "/" << filename2 << ".data";
  f2.open(sysfile.str().c_str(), ios::in);

  if (f2.fail())
  {
    cerr << "Error opening data file " << sysfile.str() << endl;
    exit(-1);
  }

  // Start filter of data
  string line;

  // Skip over first ten lines
  for (int i=0; i<10; i++)
  {
    getline(f1,line);
  }
  
  for (int i=0; i<fNData; i++)
  {
    double pt_top; // Transverse momentum of top quark
    double pt_top_high, pt_top_low; // Limits of bin
    string unneeded_info;

    getline(f1,line);
    istringstream lstream(line);

    lstream >> pt_top_low >> unneeded_info >> pt_top_high;
    pt_top = 0.5*(pt_top_low + pt_top_high);

    // Skip over next eight elements of line
    for (int j=0; j<8; j++)
    {
      lstream >> unneeded_info;
    }

    fKin1[i] = pt_top;
    fKin2[i] = Mt*Mt; // Top mass squared
    fKin3[i] = 7000; // Centre of mass energy in GeV

    lstream >> fData[i]; // Value of bin
    fData[i] /= 1000; // Convert to 1/GeV
    // We do not read in stat. correlation matrix in this filter because it is not positive semi-definite
    lstream >> unneeded_info >> fStat[i]; // Statistical (percentage) uncertainty
    fStat[i] = fStat[i]*fData[i]/100.;
  }

  // Read file with systematic uncertainty breakdown
  // Skip over first 20 lines (including stat. uncert.)
  for (int i=0; i<20; i++)
  {
    getline(f2,line);
  }
  
  double sys1, sys2, up, down, sigma, datshift;
  std::vector<double> shift(fNData);

  for (int j=0; j<fNSys; j++)
  {
    string unneeded_info;

    getline(f2,line);
    istringstream lstream(line);

    if (j==2) // Deal with aysymmetric errors
    {
      for (int i=0; i<3; i++)
      {
        lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
        if (sys1 < 0)
        {
          up=sys2;
          down=sys1;
        }
        else
        {
          up=sys1;
          down=sys2;
        }
        symmetriseErrors(up, down, &sigma, &datshift);
        fSys[i][j].mult = sigma;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";    
        shift[i] += datshift;
      }

      lstream >> fSys[3][j].mult >> unneeded_info;
      fSys[3][j].add = fSys[3][j].mult*fData[3]/100;
      fSys[3][j].type = MULT;
      fSys[3][j].name = "CORR";

      lstream >> sys1 >> unneeded_info >> sys2;
      if (sys1 < 0)
      {
        up=sys2;
        down=sys1;
      }
      else
      {
        up=sys1;
        down=sys2;
      }
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[4][j].mult = sigma;
      fSys[4][j].add = fSys[4][j].mult*fData[4]/100;
      fSys[4][j].type = MULT;
      fSys[4][j].name = "CORR";
      shift[4] += datshift;
    }
    else if (j==3) // Deal with asymmetric errors
    {
      for (int i=0; i<3; i++)
      {
        lstream >> fSys[i][j].mult >> unneeded_info;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }
 
      lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
      if (sys1 < 0)
      {
        up=sys2;
        down=sys1;
      }
      else
      {
        up=sys1;
        down=sys2;
      }
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[3][j].mult = sigma;
      fSys[3][j].add = fSys[3][j].mult*fData[3]/100;
      fSys[3][j].type = MULT;
      fSys[3][j].name = "CORR";
      shift[3] += datshift;

      lstream >> fSys[4][j].mult;
      fSys[4][j].add = fSys[4][j].mult*fData[4]/100;
      fSys[4][j].type = MULT;
      fSys[4][j].name = "CORR";
    }
    else if (j==6) // Deal with asymmetric errors
    {
      for (int i=0; i<2; i++)
      {
        lstream >> fSys[i][j].mult >> unneeded_info;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }

      for (int i=2; i<4; i++)
      {
        lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
        if (sys1 < 0)
        {
          up=sys2;
          down=sys1;
        }
        else
        {
          up=sys1;
          down=sys2;
        }
        symmetriseErrors(up, down, &sigma, &datshift);
        fSys[i][j].mult = sigma;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
        shift[i] += datshift;
      }

      lstream >> fSys[4][j].mult;
      fSys[4][j].add = fSys[4][j].mult*fData[4]/100;
      fSys[4][j].type = MULT;
      fSys[4][j].name = "CORR";
    }
    else if (j==7) // Deal with asymmetric errors
    {
      lstream >> fSys[0][j].mult >> unneeded_info;
      fSys[0][j].add = fSys[0][j].mult*fData[0]/100;
      fSys[0][j].type = MULT;
      fSys[0][j].name = "CORR";

      lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
      if (sys1 < 0)
      {
        up=sys2;
        down=sys1;
      }
      else
      {
        up=sys1;
        down=sys2;
      }
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[1][j].mult = sigma;
      fSys[1][j].add = fSys[1][j].mult*fData[1]/100;
      fSys[1][j].type = MULT;
      fSys[1][j].name = "CORR";
      shift[1] += datshift;

      lstream >> fSys[2][j].mult >> unneeded_info;
      fSys[2][j].add = fSys[2][j].mult*fData[2]/100;
      fSys[2][j].type = MULT;
      fSys[2][j].name = "CORR";

      for (int i=3; i<fNData; i++)
      {
        lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
        if (sys1 < 0)
        {
          up=sys2;
          down=sys1;
        }
        else
        {
          up=sys1;
          down=sys2;
        }
        symmetriseErrors(up, down, &sigma, &datshift);
        fSys[i][j].mult = sigma;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
        shift[i] += datshift;
      }
    }
    else // Deal with lines that contain no asymmetric errors
    {
      for (int i=0; i<fNData; i++)
      {
        lstream >> fSys[i][j].mult >> unneeded_info;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }
    }
  }

  for (int i=0; i<fNData; i++)
  {
    fData[i] *= (1.0 + shift[i]*0.01); // Shift of central value due to asymmetric errors
  }
}

//==================================================================

// 4) Distribution differential in antitop quark transverse momentum
void ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_PT_NORMFilter::ReadData()
{
  // Create streams to read data files
  fstream f1, f2, f3;

  // Data files
  stringstream datafile("");
  string filename1;
  filename1 = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_PT";
  datafile << dataPath()
           << "rawdata/" << filename1 << "/" << filename1 << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream sysfile("");
  string filename2;
  filename2 = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_PT_NORM_SYS_BREAKDOWN";
  sysfile << dataPath()
          << "rawdata/" << filename1 << "/" << filename2 << ".data";
  f2.open(sysfile.str().c_str(), ios::in);

  if (f2.fail())
  {
    cerr << "Error opening data file " << sysfile.str() << endl;
    exit(-1);
  }

  // Open correlation matrix file
  stringstream corrfile("");
  string filename3;
  filename3 = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_PT_NORM";
  corrfile << dataPath()
           << "rawdata/" << filename1 << "/" << filename3 << ".corr";
  f3.open(corrfile.str().c_str(), ios::in);

  if (f3.fail())
  {
    cerr << "Error opening data file " << corrfile.str() << endl;
    exit(-1);
  }

  // Start filter of data
  string line;

  // Initialise array to store additive stat. uncerts.
  std::vector<double> fstat_additive(fNData);

  // Skip over first ten lines
  for (int i=0; i<10; i++)
  {
    getline(f1,line);
  }
  
  for (int i=0; i<fNData; i++)
  {
    double pt_top; // Transverse momentum of top quark
    double pt_top_high, pt_top_low; // Limits of bin
    double fstat_percentage; // Percentage statistical uncertainty
    string unneeded_info;

    getline(f1,line);
    istringstream lstream(line);

    lstream >> pt_top_low >> unneeded_info >> pt_top_high;
    pt_top = 0.5*(pt_top_low + pt_top_high);

    // Skip over next eight elements of line
    for (int j=0; j<8; j++)
    {
      lstream >> unneeded_info;
    }

    fKin1[i] = pt_top;
    fKin2[i] = Mt*Mt; // Top mass squared
    fKin3[i] = 7000; // Centre of mass energy in GeV

    lstream >> fData[i]; // Value of bin
    fData[i] /= 1000; // Convert to 1/GeV
    lstream >> unneeded_info >> fstat_percentage;
    fstat_additive[i] = fstat_percentage*fData[i]/100;

    fStat[i] = 1e-10; // Set stat. error to zero to avoid double counting when using artificial systematics
  }

  // Read statistical correlation matrix
  // Skip over first ten lines
  for (int i=0; i<10; i++)
  {
    getline(f3,line);
  }

  double** covmat = new double*[fNData];
  NNPDF::matrix<double> corrmat(fNData, fNData);
  for (int i=0; i<fNData; i++)
  {
    string unneeded_info;
    covmat[i] = new double[fNData];
    getline(f3,line);
    istringstream lstream(line);
    lstream >> unneeded_info >> unneeded_info >> unneeded_info;
    for (int j=0; j<fNData; j++)
    {
      lstream >> corrmat(i,j) >> unneeded_info;
      covmat[i][j] = corrmat(i,j) * fstat_additive[i] * fstat_additive[j];
    }
  }

  // Generate artificial systematics
  double** syscor = new double*[fNData];
  for (int i=0; i<fNData; i++)
    syscor[i] = new double[fNData];

  if (!genArtSys(fNData,covmat,syscor))
  {
    cerr << " in " << fSetName << endl;
    exit(-1);
  }

  // Assign artificial systematics
  for (int i=0; i<fNData; i++)
  {
    for (int j=0; j<fNData; j++)
      {
        fSys[i][j].add = syscor[i][j];
        fSys[i][j].mult = fSys[i][j].add*100/fData[i];
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }
  }

  // Read file with systematic uncertainty breakdown
  // Skip over first 20 lines (including stat. uncert.)
  for (int i=0; i<20; i++)
  {
    getline(f2,line);
  }
  
  double sys1, sys2, up, down, sigma, datshift;
  std::vector<double> shift(fNData);

  for (int j=fNData; j<fNSys; j++)
  {
    string unneeded_info;

    getline(f2,line);
    istringstream lstream(line);

    if (j==fNData+2) // Deal with aysymmetric errors
    {
      for (int i=0; i<(fNData-1); i++)
      {
        lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
        if (sys1 < 0)
        {
          up=sys2;
          down=sys1;
        }
        else
        {
          up=sys1;
          down=sys2;
        }
        symmetriseErrors(up, down, &sigma, &datshift);
        fSys[i][j].mult = sigma;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";    
        shift[i] += datshift; 
      }

      lstream >> fSys[4][j].mult >> unneeded_info;
      fSys[4][j].add = fSys[4][j].mult*fData[4]/100;
      fSys[4][j].type = MULT;
      fSys[4][j].name = "CORR";
    }
    else if (j==fNData+3) // Deal with asymmetric errors
    {
      for (int i=0; i<2; i++)
      {
        lstream >> fSys[i][j].mult >> unneeded_info;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }

      lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
      if (sys1 < 0)
      {
        up=sys2;
        down=sys1;
      }
      else
      {
        up=sys1;
        down=sys2;
      }
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[2][j].mult = sigma;
      fSys[2][j].add = fSys[2][j].mult*fData[2]/100;
      fSys[2][j].type = MULT;
      fSys[2][j].name = "CORR";    
      shift[2] += datshift;

      lstream >> fSys[3][j].mult >> unneeded_info;
      fSys[3][j].add = fSys[3][j].mult*fData[3]/100;
      fSys[3][j].type = MULT;
      fSys[3][j].name = "CORR";

      lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
      if (sys1 < 0)
      {
        up=sys2;
        down=sys1;
      }
      else
      {
        up=sys1;
        down=sys2;
      }
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[4][j].mult = sigma;
      fSys[4][j].add = fSys[4][j].mult*fData[4]/100;
      fSys[4][j].type = MULT;
      fSys[4][j].name = "CORR";    
      shift[4] += datshift;
    }
    else if (j==fNData+4) // Deal with asymmetric errors
    {
      for (int i=0; i<(fNData-1); i++)
      {
        lstream >> fSys[i][j].mult >> unneeded_info;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }

      lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
      if (sys1 < 0)
      {
        up=sys2;
        down=sys1;
      }
      else
      {
        up=sys1;
        down=sys2;
      }
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[4][j].mult = sigma;
      fSys[4][j].add = fSys[4][j].mult*fData[4]/100;
      fSys[4][j].type = MULT;
      fSys[4][j].name = "CORR";    
      shift[4] += datshift;
    }
    else if (j==fNData+5) // Deal with asymmetric errors
    {
      lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
      if (sys1 < 0)
      {
        up=sys2;
        down=sys1;
      }
      else
      {
        up=sys1;
        down=sys2;
      }
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[0][j].mult = sigma;
      fSys[0][j].add = fSys[0][j].mult*fData[0]/100;
      fSys[0][j].type = MULT;
      fSys[0][j].name = "CORR";    
      shift[0] += datshift;

      for (int i=1; i<3; i++)
      {
        lstream >> fSys[i][j].mult >> unneeded_info;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }

      for (int i=3; i<fNData; i++)
      {
        lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
        if (sys1 < 0)
        {
          up=sys2;
          down=sys1;
        }
        else
        {
          up=sys1;
          down=sys2;
        }
        symmetriseErrors(up, down, &sigma, &datshift);
        fSys[i][j].mult = sigma;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";    
        shift[i] += datshift;
      }
    }
    else if (j==fNData+7) // Deal with asymmetric errors
    {
      for (int i=0; i<fNData; i++)
      {
        lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
        if (sys1 < 0)
        {
          up=sys2;
          down=sys1;
        }
        else
        {
          up=sys1;
          down=sys2;
        }
        symmetriseErrors(up, down, &sigma, &datshift);
        fSys[i][j].mult = sigma;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";    
        shift[i] += datshift; 
      }
    }
    else if (j==fNData+8) // Deal with asymmetric errors
    {
      lstream >> fSys[0][j].mult >> unneeded_info;
      fSys[0][j].add = fSys[0][j].mult*fData[0]/100;
      fSys[0][j].type = MULT;
      fSys[0][j].name = "CORR";

      lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
      if (sys1 < 0)
      {
        up=sys2;
        down=sys1;
      }
      else
      {
        up=sys1;
        down=sys2;
      }
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[1][j].mult = sigma;
      fSys[1][j].add = fSys[1][j].mult*fData[1]/100;
      fSys[1][j].type = MULT;
      fSys[1][j].name = "CORR";    
      shift[1] += datshift;

      lstream >> fSys[2][j].mult >> unneeded_info;
      fSys[2][j].add = fSys[2][j].mult*fData[2]/100;
      fSys[2][j].type = MULT;
      fSys[2][j].name = "CORR";

      lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
      if (sys1 < 0)
      {
        up=sys2;
        down=sys1;
      }
      else
      {
        up=sys1;
        down=sys2;
      }
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[3][j].mult = sigma;
      fSys[3][j].add = fSys[3][j].mult*fData[3]/100;
      fSys[3][j].type = MULT;
      fSys[3][j].name = "CORR";    
      shift[3] += datshift;

      lstream >> fSys[4][j].mult >> unneeded_info;
      fSys[4][j].add = fSys[4][j].mult*fData[4]/100;
      fSys[4][j].type = MULT;
      fSys[4][j].name = "CORR";
    }
    else // Deal with lines that contain no asymmetric errors
    {
      for (int i=0; i<fNData; i++)
      {
        lstream >> fSys[i][j].mult >> unneeded_info;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }
    }
  }

  for (int i=0; i<fNData; i++)
  {
    fData[i] *= (1.0 + shift[i]*0.01); // Shift of central value due to asymmetric errors
  }

  // Clean-up
  for (int i=0; i<fNData; i++)
    delete[] syscor[i];

  delete[] syscor;

  for(int i=0; i<fNData; i++)
    delete[] covmat[i];

  delete[] covmat;
}

//==================================================================

// B - UNNORMALISED distributions

// 5) Distribution differential in modulus of top quark rapidity
void ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_RAPFilter::ReadData()
{
  // Create streams to read data files
  fstream f1, f2, f3;

  // Data files
  stringstream datafile("");
  string filename1;
  filename1 = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_RAP";
  datafile << dataPath()
           << "rawdata/" << filename1 << "/" << filename1 << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream sysfile("");
  string filename2;
  filename2 = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_RAP_SYS_BREAKDOWN";
  sysfile << dataPath()
          << "rawdata/" << filename1 << "/" << filename2 << ".data";
  f2.open(sysfile.str().c_str(), ios::in);

  if (f2.fail())
  {
    cerr << "Error opening data file " << sysfile.str() << endl;
    exit(-1);
  }

  // Open correlation matrix file
  stringstream corrfile("");
  string filename3;
  filename3 = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_RAP";
  corrfile << dataPath()
           << "rawdata/" << filename1 << "/" << filename3 << ".corr";
  f3.open(corrfile.str().c_str(), ios::in);

  if (f3.fail())
  {
    cerr << "Error opening data file " << corrfile.str() << endl;
    exit(-1);
  }

  // Start filter of data
  string line;

  // Initialise array to store additive stat. uncerts.
  std::vector<double> fstat_additive(fNData);

  // Skip over first ten lines
  for (int i=0; i<10; i++)
  {
    getline(f1,line);
  }
  
  for (int i=0; i<fNData; i++)
  {
    double rap_top; // Rapidity of top quark
    double rap_top_low, rap_top_high; // Limits of bin
    double fstat_percentage; // Percentage statistical uncertainty
    string unneeded_info;

    getline(f1,line);
    istringstream lstream(line);

    lstream >> rap_top_low >> unneeded_info >> rap_top_high >> unneeded_info;
    rap_top = 0.5*(rap_top_low + rap_top_high);

    fKin1[i] = rap_top;
    fKin2[i] = Mt*Mt; // Top mass squared
    fKin3[i] = 7000; // Centre of mass energy in GeV

    lstream >> fData[i]; // Value of bin
    lstream >> unneeded_info >> fstat_percentage;
    fstat_additive[i] = fstat_percentage*fData[i]/100;

    fStat[i] = 1e-10; // Set stat. error to zero to avoid double counting when using artificial systematics
  }

  // Read statistical correlation matrix
  // Skip over first ten lines
  for (int i=0; i<10; i++)
  {
    getline(f3,line);
  }

  double** covmat = new double*[fNData];
  NNPDF::matrix<double> corrmat(fNData, fNData);
  for (int i=0; i<fNData; i++)
  {
    string unneeded_info;
    covmat[i] = new double[fNData];
    getline(f3,line);
    istringstream lstream(line);
    lstream >> unneeded_info >> unneeded_info >> unneeded_info;
    for (int j=0; j<fNData; j++)
    {
      lstream >> corrmat(i,j) >> unneeded_info;
      covmat[i][j] = corrmat(i,j) * fstat_additive[i] * fstat_additive[j];
    }
  }

  // Generate artificial systematics
  double** syscor = new double*[fNData];
  for (int i=0; i<fNData; i++)
    syscor[i] = new double[fNData];

  if (!genArtSys(fNData,covmat,syscor))
  {
    cerr << " in " << fSetName << endl;
    exit(-1);
  }

  // Assign artificial systematics
  for (int i=0; i<fNData; i++)
  {
    for (int j=0; j<fNData; j++)
      {
        fSys[i][j].add = syscor[i][j];
        fSys[i][j].mult = fSys[i][j].add*100/fData[i];
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }
  }

  // Read file with systematic uncertainty breakdown
  // Skip over first 19 lines (including stat. uncert.)
  for (int i=0; i<19; i++)
  {
    getline(f2,line);
  }
  
  double sys1, sys2, up, down, sigma, datshift;
  std::vector<double> shift(fNData);

  for (int j=fNData; j<fNSys; j++)
  {
    string unneeded_info;

    getline(f2,line);
    istringstream lstream(line);

    if (j==fNData+2 || j==fNData+3) // Deal with aysymmetric errors
    {
      lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
      if (sys1 < 0)
      {
        up=sys2;
        down=sys1;
      }
      else
      {
        up=sys1;
        down=sys2;
      }
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[0][j].mult = sigma;
      fSys[0][j].add = fSys[0][j].mult*fData[0]/100;
      fSys[0][j].type = MULT;
      fSys[0][j].name = "CORR";    
      shift[0] += datshift;

      lstream >> fSys[1][j].mult >> unneeded_info;
      fSys[1][j].add = fSys[1][j].mult*fData[1]/100;
      fSys[1][j].type = MULT;
      fSys[1][j].name = "CORR";

      lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
      if (sys1 < 0)
      {
        up=sys2;
        down=sys1;
      }
      else
      {
        up=sys1;
        down=sys2;
      }
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[2][j].mult = sigma;
      fSys[2][j].add = fSys[2][j].mult*fData[2]/100;
      fSys[2][j].type = MULT;
      fSys[2][j].name = "CORR";    
      shift[2] += datshift;

      lstream >> fSys[3][j].mult >> unneeded_info;
      fSys[3][j].add = fSys[3][j].mult*fData[3]/100;
      fSys[3][j].type = MULT;
      fSys[3][j].name = "CORR";
    }
    else if (j==fNData+10) // Deal with asymmetric errors
    {
      lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
      if (sys1 < 0)
      {
        up=sys2;
        down=sys1;
      }
      else
      {
        up=sys1;
        down=sys2;
      }
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[0][j].mult = sigma;
      fSys[0][j].add = fSys[0][j].mult*fData[0]/100;
      fSys[0][j].type = MULT;
      fSys[0][j].name = "CORR";    
      shift[0] += datshift;

      for (int i=1; i<fNData; i++)
      {
        lstream >> fSys[i][j].mult >> unneeded_info;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }
    }
    else if (j==fNData+12) // Deal with luminosity uncertainties
    {
      for (int i=0; i<fNData; i++)
      {
        lstream >> fSys[i][j].mult >> unneeded_info;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "ATLASLUMI11";
      }
    }
    else // Deal with lines that contain no asymmetric errors
    {
      for (int i=0; i<fNData; i++)
      {
        lstream >> fSys[i][j].mult >> unneeded_info;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }
    }
  }

  for (int i=0; i<fNData; i++)
  {
    fData[i] *= (1.0 + shift[i]*0.01); // Shift of central value due to asymmetric errors
  }

  // Clean-up
  for (int i=0; i<fNData; i++)
    delete[] syscor[i];

  delete[] syscor;

  for(int i=0; i<fNData; i++)
    delete[] covmat[i];

  delete[] covmat;
}

//==================================================================

// 6) Distribution differential in modulus of antitop quark rapidity
void ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_RAPFilter::ReadData()
{
  // Create streams to read data files
  fstream f1, f2, f3;

  // Data files
  stringstream datafile("");
  string filename1;
  filename1 = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_RAP";
  datafile << dataPath()
           << "rawdata/" << filename1 << "/" << filename1 << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream sysfile("");
  string filename2;
  filename2 = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_RAP_SYS_BREAKDOWN";
  sysfile << dataPath()
          << "rawdata/" << filename1 << "/" << filename2 << ".data";
  f2.open(sysfile.str().c_str(), ios::in);

  if (f2.fail())
  {
    cerr << "Error opening data file " << sysfile.str() << endl;
    exit(-1);
  }

  // Open correlation matrix file
  stringstream corrfile("");
  string filename3;
  filename3 = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_RAP";
  corrfile << dataPath()
           << "rawdata/" << filename1 << "/" << filename3 << ".corr";
  f3.open(corrfile.str().c_str(), ios::in);

  if (f3.fail())
  {
    cerr << "Error opening data file " << corrfile.str() << endl;
    exit(-1);
  }

  // Start filter of data
  string line;

  // Initialise array to store additive stat. uncerts.
  std::vector<double> fstat_additive(fNData);

  // Skip over first ten lines
  for (int i=0; i<10; i++)
  {
    getline(f1,line);
  }
  
  for (int i=0; i<fNData; i++)
  {
    double rap_top; // Rapidity of top quark
    double rap_top_low, rap_top_high; // Limits of bin
    double fstat_percentage; // Percentage statistical uncertainty
    string unneeded_info;

    getline(f1,line);
    istringstream lstream(line);

    lstream >> rap_top_low >> unneeded_info >> rap_top_high >> unneeded_info;
    rap_top = 0.5*(rap_top_low + rap_top_high);

    fKin1[i] = rap_top;
    fKin2[i] = Mt*Mt; // Top mass squared
    fKin3[i] = 7000; // Centre of mass energy in GeV

    lstream >> fData[i]; // Value of bin
    lstream >> unneeded_info >> fstat_percentage;
    fstat_additive[i] = fstat_percentage*fData[i]/100;

    fStat[i] = 1e-10; // Set stat. error to zero to avoid double counting when using artificial systematics
  }

  // Read statistical correlation matrix
  // Skip over first ten lines
  for (int i=0; i<10; i++)
  {
    getline(f3,line);
  }

  double** covmat = new double*[fNData];
  NNPDF::matrix<double> corrmat(fNData, fNData);
  for (int i=0; i<fNData; i++)
  {
    string unneeded_info;
    covmat[i] = new double[fNData];
    getline(f3,line);
    istringstream lstream(line);
    lstream >> unneeded_info >> unneeded_info >> unneeded_info;
    for (int j=0; j<fNData; j++)
    {
      lstream >> corrmat(i,j) >> unneeded_info;
      covmat[i][j] = corrmat(i,j) * fstat_additive[i] * fstat_additive[j];
    }
  }

  // Generate artificial systematics
  double** syscor = new double*[fNData];
  for (int i=0; i<fNData; i++)
    syscor[i] = new double[fNData];

  if (!genArtSys(fNData,covmat,syscor))
  {
    cerr << " in " << fSetName << endl;
    exit(-1);
  }

  // Assign artificial systematics
  for (int i=0; i<fNData; i++)
  {
    for (int j=0; j<fNData; j++)
      {
        fSys[i][j].add = syscor[i][j];
        fSys[i][j].mult = fSys[i][j].add*100/fData[i];
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }
  }

  // Read file with systematic uncertainty breakdown
  // Skip over first 19 lines (including stat. uncert.)
  for (int i=0; i<19; i++)
  {
    getline(f2,line);
  }
  
  double sys1, sys2, up, down, sigma, datshift;
  std::vector<double> shift(fNData);

  for (int j=fNData; j<fNSys; j++)
  {
    string unneeded_info;

    getline(f2,line);
    istringstream lstream(line);

    if (j==fNData+2) // Deal with aysymmetric errors
    {
      lstream >> fSys[0][j].mult >> unneeded_info;
      fSys[0][j].add = fSys[0][j].mult*fData[0]/100;
      fSys[0][j].type = MULT;
      fSys[0][j].name = "CORR";

      lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
      if (sys1 < 0)
      {
        up=sys2;
        down=sys1;
      }
      else
      {
        up=sys1;
        down=sys2;
      }
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[1][j].mult = sigma;
      fSys[1][j].add = fSys[1][j].mult*fData[1]/100;
      fSys[1][j].type = MULT;
      fSys[1][j].name = "CORR";    
      shift[1] += datshift;

      lstream >> fSys[2][j].mult >> unneeded_info;
      fSys[2][j].add = fSys[2][j].mult*fData[2]/100;
      fSys[2][j].type = MULT;
      fSys[2][j].name = "CORR";

      lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
      if (sys1 < 0)
      {
        up=sys2;
        down=sys1;
      }
      else
      {
        up=sys1;
        down=sys2;
      }
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[3][j].mult = sigma;
      fSys[3][j].add = fSys[3][j].mult*fData[3]/100;
      fSys[3][j].type = MULT;
      fSys[3][j].name = "CORR";    
      shift[3] += datshift;
    }
    else if (j==fNData+8) // Deal with asymmetric errors
    {
      lstream >> fSys[0][j].mult >> unneeded_info;
      fSys[0][j].add = fSys[0][j].mult*fData[0]/100;
      fSys[0][j].type = MULT;
      fSys[0][j].name = "CORR";

      lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
      if (sys1 < 0)
      {
        up=sys2;
        down=sys1;
      }
      else
      {
        up=sys1;
        down=sys2;
      }
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[1][j].mult = sigma;
      fSys[1][j].add = fSys[1][j].mult*fData[1]/100;
      fSys[1][j].type = MULT;
      fSys[1][j].name = "CORR";    
      shift[1] += datshift;

      for (int i=2; i<fNData; i++)
      {
        lstream >> fSys[i][j].mult >> unneeded_info;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }
    }
    else if (j==fNData+10) // Deal with luminosity uncertainties
    {
      for (int i=0; i<fNData; i++)
      {
        lstream >> fSys[i][j].mult >> unneeded_info;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "ATLASLUMI11";
      }
    }
    else // Deal with lines that contain no asymmetric errors
    {
      for (int i=0; i<fNData; i++)
      {
        lstream >> fSys[i][j].mult >> unneeded_info;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }
    }
  }

  for (int i=0; i<fNData; i++)
  {
    fData[i] *= (1.0 + shift[i]*0.01); // Shift of central value due to asymmetric errors
  }

  // Clean-up
  for (int i=0; i<fNData; i++)
    delete[] syscor[i];

  delete[] syscor;

  for(int i=0; i<fNData; i++)
    delete[] covmat[i];

  delete[] covmat;
}

//==================================================================

// 7) Distribution differential in top quark transverse momentum
void ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_PTFilter::ReadData()
{
  // Create streams to read data files
  fstream f1, f2, f3;

  // Data files
  stringstream datafile("");
  string filename1;
  filename1 = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_PT";
  datafile << dataPath()
           << "rawdata/" << filename1 << "/" << filename1 << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream sysfile("");
  string filename2;
  filename2 = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_PT_SYS_BREAKDOWN";
  sysfile << dataPath()
          << "rawdata/" << filename1 << "/" << filename2 << ".data";
  f2.open(sysfile.str().c_str(), ios::in);

  if (f2.fail())
  {
    cerr << "Error opening data file " << sysfile.str() << endl;
    exit(-1);
  }

  // Open correlation matrix file
  stringstream corrfile("");
  string filename3;
  filename3 = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_PT";
  corrfile << dataPath()
           << "rawdata/" << filename1 << "/" << filename3 << ".corr";
  f3.open(corrfile.str().c_str(), ios::in);

  if (f3.fail())
  {
    cerr << "Error opening data file " << corrfile.str() << endl;
    exit(-1);
  }

  // Start filter of data
  string line;

  // Initialise array to store additive stat. uncerts.
  std::vector<double> fstat_additive(fNData);

  // Skip over first ten lines
  for (int i=0; i<10; i++)
  {
    getline(f1,line);
  }
  
  for (int i=0; i<fNData; i++)
  {
    double pt_top; // Transverse momentum of top quark
    double pt_top_high, pt_top_low; // Limits of bin
    double fstat_percentage; // Percentage statistical uncertainty
    string unneeded_info;

    getline(f1,line);
    istringstream lstream(line);

    lstream >> pt_top_low >> unneeded_info >> pt_top_high >> unneeded_info;
    pt_top = 0.5*(pt_top_low + pt_top_high);

    fKin1[i] = pt_top;
    fKin2[i] = Mt*Mt; // Top mass squared
    fKin3[i] = 7000; // Centre of mass energy in GeV

    lstream >> fData[i]; // Value of bin
    fData[i] /= 1000; // Convert to fb/GeV
    lstream >> unneeded_info >> fstat_percentage;
    fstat_additive[i] = fstat_percentage*fData[i]/100;

    fStat[i] = 1e-10; // Set stat. error to zero to avoid double counting when using artificial systematics
  }

  // Read statistical correlation matrix
  // Skip over first ten lines
  for (int i=0; i<10; i++)
  {
    getline(f3,line);
  }

  double** covmat = new double*[fNData];
  NNPDF::matrix<double> corrmat(fNData, fNData);
  for (int i=0; i<fNData; i++)
  {
    string unneeded_info;
    covmat[i] = new double[fNData];
    getline(f3,line);
    istringstream lstream(line);
    lstream >> unneeded_info >> unneeded_info >> unneeded_info;
    for (int j=0; j<fNData; j++)
    {
      lstream >> corrmat(i,j) >> unneeded_info;
      covmat[i][j] = corrmat(i,j) * fstat_additive[i] * fstat_additive[j];
    }
  }

  // Generate artificial systematics
  double** syscor = new double*[fNData];
  for (int i=0; i<fNData; i++)
    syscor[i] = new double[fNData];

  if (!genArtSys(fNData,covmat,syscor))
  {
    cerr << " in " << fSetName << endl;
    exit(-1);
  }

  // Assign artificial systematics
  for (int i=0; i<fNData; i++)
  {
    for (int j=0; j<fNData; j++)
      {
        fSys[i][j].add = syscor[i][j];
        fSys[i][j].mult = fSys[i][j].add*100/fData[i];
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }
  }

  // Read file with systematic uncertainty breakdown
  // Skip over first 19 lines (including stat. uncert.)
  for (int i=0; i<19; i++)
  {
    getline(f2,line);
  }
  
  double sys1, sys2, up, down, sigma, datshift;
  std::vector<double> shift(fNData);

  for (int j=fNData; j<fNSys; j++)
  {
    string unneeded_info;

    getline(f2,line);
    istringstream lstream(line);

    if (j==fNData+2) // Deal with aysymmetric errors
    {
      lstream >> fSys[0][j].mult >> unneeded_info;
      fSys[0][j].add = fSys[0][j].mult*fData[0]/100;
      fSys[0][j].type = MULT;
      fSys[0][j].name = "CORR";

      for (int i=1; i<3; i++)
      {
        lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
        if (sys1 < 0)
        {
          up=sys2;
          down=sys1;
        }
        else
        {
          up=sys1;
          down=sys2;
        }
        symmetriseErrors(up, down, &sigma, &datshift);
        fSys[i][j].mult = sigma;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";    
        shift[i] += datshift; 
      }

      lstream >> fSys[3][j].mult >> unneeded_info;
      fSys[3][j].add = fSys[3][j].mult*fData[3]/100;
      fSys[3][j].type = MULT;
      fSys[3][j].name = "CORR";

      lstream >> sys1 >> unneeded_info >> sys2;
      if (sys1 < 0)
      {
        up=sys2;
        down=sys1;
      }
      else
      {
        up=sys1;
        down=sys2;
      }
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[4][j].mult = sigma;
      fSys[4][j].add = fSys[4][j].mult*fData[4]/100;
      fSys[4][j].type = MULT;
      fSys[4][j].name = "CORR";    
      shift[4] += datshift; 
    }
    else if (j==fNData+3) // Deal with aysymmetric errors
    {
      lstream >> fSys[0][j].mult >> unneeded_info;
      fSys[0][j].add = fSys[0][j].mult*fData[0]/100;
      fSys[0][j].type = MULT;
      fSys[0][j].name = "CORR";

      lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
      if (sys1 < 0)
      {
        up=sys2;
        down=sys1;
      }
      else
      {
        up=sys1;
        down=sys2;
      }
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[1][j].mult = sigma;
      fSys[1][j].add = fSys[1][j].mult*fData[1]/100;
      fSys[1][j].type = MULT;
      fSys[1][j].name = "CORR";    
      shift[1] += datshift; 

      for (int i=2; i<fNData; i++)
      {
        lstream >> fSys[i][j].mult >> unneeded_info;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }
    }
    else if (j==fNData+10) // Deal with asymmetric errors
    {
      for (int i=0; i<2; i++)
      {
        lstream >> fSys[i][j].mult >> unneeded_info;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }

      for (int i=2; i<fNData; i++)
      {
        lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
        if (sys1 < 0)
        {
          up=sys2;
          down=sys1;
        }
        else
        {
          up=sys1;
          down=sys2;
        }
        symmetriseErrors(up, down, &sigma, &datshift);
        fSys[i][j].mult = sigma;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";    
        shift[i] += datshift;
      }
    }
    else if (j==fNData+11) // Deal with asymmetric errors
    {
      for (int i=0; i<(fNData-1); i++)
      {
        lstream >> fSys[i][j].mult >> unneeded_info;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }

      lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
      if (sys1 < 0)
      {
        up=sys2;
        down=sys1;
      }
      else
      {
        up=sys1;
        down=sys2;
      }
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[4][j].mult = sigma;
      fSys[4][j].add = fSys[4][j].mult*fData[4]/100;
      fSys[4][j].type = MULT;
      fSys[4][j].name = "CORR";    
      shift[4] += datshift;
    }
    else if (j==fNData+13) // Deal with luminosity uncertainties
    {
      for (int i=0; i<fNData; i++)
      {
        lstream >> fSys[i][j].mult >> unneeded_info;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "ATLASLUMI11";
      }
    }
    else // Deal with lines that contain no asymmetric errors
    {
      for (int i=0; i<fNData; i++)
      {
        lstream >> fSys[i][j].mult >> unneeded_info;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }
    }
  }

  for (int i=0; i<fNData; i++)
  {
    fData[i] *= (1.0 + shift[i]*0.01); // Shift of central value due to asymmetric errors
  }

  // Clean-up
  for (int i=0; i<fNData; i++)
    delete[] syscor[i];

  delete[] syscor;

  for(int i=0; i<fNData; i++)
    delete[] covmat[i];

  delete[] covmat;
}

//==================================================================

// 8) Distribution differential in antitop quark transverse momentum
void ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_PTFilter::ReadData()
{
  // Create streams to read data files
  fstream f1, f2, f3;

  // Data files
  stringstream datafile("");
  string filename1;
  filename1 = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_PT";
  datafile << dataPath()
           << "rawdata/" << filename1 << "/" << filename1 << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream sysfile("");
  string filename2;
  filename2 = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_PT_SYS_BREAKDOWN";
  sysfile << dataPath()
          << "rawdata/" << filename1 << "/" << filename2 << ".data";
  f2.open(sysfile.str().c_str(), ios::in);

  if (f2.fail())
  {
    cerr << "Error opening data file " << sysfile.str() << endl;
    exit(-1);
  }

  // Open correlation matrix file
  stringstream corrfile("");
  string filename3;
  filename3 = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_PT";
  corrfile << dataPath()
           << "rawdata/" << filename1 << "/" << filename3 << ".corr";
  f3.open(corrfile.str().c_str(), ios::in);

  if (f3.fail())
  {
    cerr << "Error opening data file " << corrfile.str() << endl;
    exit(-1);
  }

  // Start filter of data
  string line;

  // Initialise array to store additive stat. uncerts.
  std::vector<double> fstat_additive(fNData);

  // Skip over first ten lines
  for (int i=0; i<10; i++)
  {
    getline(f1,line);
  }
  
  for (int i=0; i<fNData; i++)
  {
    double pt_top; // Transverse momentum of top quark
    double pt_top_high, pt_top_low; // Limits of bin
    double fstat_percentage; // Percentage statistical uncertainty
    string unneeded_info;

    getline(f1,line);
    istringstream lstream(line);

    lstream >> pt_top_low >> unneeded_info >> pt_top_high >> unneeded_info;
    pt_top = 0.5*(pt_top_low + pt_top_high);

    fKin1[i] = pt_top;
    fKin2[i] = Mt*Mt; // Top mass squared
    fKin3[i] = 7000; // Centre of mass energy in GeV

    lstream >> fData[i]; // Value of bin
    fData[i] /= 1000; // Convert to fb/GeV
    lstream >> unneeded_info >> fstat_percentage;
    fstat_additive[i] = fstat_percentage*fData[i]/100;

    fStat[i] = 1e-10; // Set stat. error to zero to avoid double counting when using artificial systematics
  }

  // Read statistical correlation matrix
  // Skip over first ten lines
  for (int i=0; i<10; i++)
  {
    getline(f3,line);
  }

  double** covmat = new double*[fNData];
  NNPDF::matrix<double> corrmat(fNData, fNData);
  for (int i=0; i<fNData; i++)
  {
    string unneeded_info;
    covmat[i] = new double[fNData];
    getline(f3,line);
    istringstream lstream(line);
    lstream >> unneeded_info >> unneeded_info >> unneeded_info;
    for (int j=0; j<fNData; j++)
    {
      lstream >> corrmat(i,j) >> unneeded_info;
      covmat[i][j] = corrmat(i,j) * fstat_additive[i] * fstat_additive[j];
    }
  }

  // Generate artificial systematics
  double** syscor = new double*[fNData];
  for (int i=0; i<fNData; i++)
    syscor[i] = new double[fNData];

  if (!genArtSys(fNData,covmat,syscor))
  {
    cerr << " in " << fSetName << endl;
    exit(-1);
  }

  // Assign artificial systematics
  for (int i=0; i<fNData; i++)
  {
    for (int j=0; j<fNData; j++)
      {
        fSys[i][j].add = syscor[i][j];
        fSys[i][j].mult = fSys[i][j].add*100/fData[i];
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }
  }

  // Read file with systematic uncertainty breakdown
  // Skip over first 19 lines (including stat. uncert.)
  for (int i=0; i<19; i++)
  {
    getline(f2,line);
  }
  
  double sys1, sys2, up, down, sigma, datshift;
  std::vector<double> shift(fNData);

  for (int j=fNData; j<fNSys; j++)
  {
    string unneeded_info;

    getline(f2,line);
    istringstream lstream(line);

    if (j==fNData+2) // Deal with aysymmetric errors
    {
      for (int i=0; i<(fNData-1); i++)
      {
        lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
        if (sys1 < 0)
        {
          up=sys2;
          down=sys1;
        }
        else
        {
          up=sys1;
          down=sys2;
        }
        symmetriseErrors(up, down, &sigma, &datshift);
        fSys[i][j].mult = sigma;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";    
        shift[i] += datshift; 
      }

      lstream >> fSys[4][j].mult >> unneeded_info;
      fSys[4][j].add = fSys[4][j].mult*fData[4]/100;
      fSys[4][j].type = MULT;
      fSys[4][j].name = "CORR";
    }
    else if (j==fNData+8) // Deal with asymmetric errors
    {
      for (int i=0; i<2; i++)
      {
        lstream >> fSys[i][j].mult >> unneeded_info;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }

      lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
      if (sys1 < 0) {up=sys2; down=sys1;}
      else {up=sys1; down=sys2;}
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[2][j].mult = sigma;
      fSys[2][j].add = fSys[2][j].mult*fData[2]/100;
      fSys[2][j].type = MULT;
      fSys[2][j].name = "CORR";    
      shift[2] += datshift;

      for (int i=3; i<fNData; i++)
      {
        lstream >> fSys[i][j].mult >> unneeded_info;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }
    }
    else if (j==fNData+11) // Deal with luminosity uncertainties
    {
      for (int i=0; i<fNData; i++)
      {
        lstream >> fSys[i][j].mult >> unneeded_info;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "ATLASLUMI11";
      }
    }
    else // Deal with lines that contain no asymmetric errors
    {
      for (int i=0; i<fNData; i++)
      {
        lstream >> fSys[i][j].mult >> unneeded_info;
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }
    }
  }

  for (int i=0; i<fNData; i++)
  {
    fData[i] *= (1.0 + shift[i]*0.01); // Shift of central value due to asymmetric errors
  }

  // Clean-up
  for (int i=0; i<fNData; i++)
    delete[] syscor[i];

  delete[] syscor;

  for(int i=0; i<fNData; i++)
    delete[] covmat[i];

  delete[] covmat;
}
