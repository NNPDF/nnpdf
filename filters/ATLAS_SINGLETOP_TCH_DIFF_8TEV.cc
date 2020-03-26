/*
Differential cross section measurements of the single top and single antitop quark in the t-channel @LHC ATLAS 8 TeV

LHC-ATLAS 8 TeV
---------------

Selected events contain exactly one electron or muon, exactly two jets (exactly one of which must be b-tagged), and E_T^{miss} > 30 GeV.
Archived as: https://arxiv.org/pdf/1702.02859v3.pdf
Published in: Eur. Phys. J. C 77 (2017) 531

Eight distributions are implemented here. These are normalised and unnormalised
distributions differential in:
1) Top quark absolute rapidity
2) Antitop quark absolute rapidity
3) Top quark transverse momentum
4) Antitop quark transverse momentum

Description of raw data:
Tables containing central values are taken from Tables 18-21 of the paper
The rest of the information is taken from the auxiliary material: https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/TOPQ-2015-05/
The statistical correlation matrices are taken from Figures 6a-d and 7a-d
The uncertainty breakdowns are taken from Tables 19-26

Distributions are converted, where necessary, so that they have the following dimensions:
Absolute transverse momentum: pb/GeV
Absolute rapidity: pb
Normalised transverse momentum: 1/GeV
Normalised rapidity: -

Notes:
1) The number of systematic uncertainties is distribution-dependent.
2) All real systematic uncertainties are treated as being multiplicative.
3) All real systematic uncertainties, except the lumi. uncertainty, are
   treated following the prescription defined by Eq. 6 in https://arxiv.org/pdf/1703.01630.pdf.
   That is, each positive and negative variation is treated separately.
   This is done at the expense of symmetrising each asymmetric uncertainty.
4) All artificial systematic uncertainties (i.e. those generated using the
   statistical uncertainties and the statistical correlation matrix) are treated
   as being additive. Note that both the 'Data statistics' and 'Monte Carlo'
   statistics are included in the total statistical uncertainty.
5) All systematics are treated as CORR (i.e. correlated), except for the
   luminosity uncertainty for the unnormalised distributions which are treated
   as ATLASLUMI12 (i.e. ATLAS luminosity for the 2012 data set).
6) The last bin is removed from all the normalised distributions, because it
   is a linear combination of the other. This also removes the spurious
   feature of covariance matrices not being positive-semidefinite.
*/

#include "ATLAS_SINGLETOP_TCH_DIFF_8TEV.h"
#include "NNPDF/utils.h"

// A - UNNORMALISED distributions

// 1) Normalised distribution differential in modulus of top quark rapidity

void ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_RAP_NORMFilter::ReadData()
{
  // Create stream to read data file
  fstream f1, f2, f3;

  // Data files
  stringstream datafile("");
  string dirname;
  dirname = "ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_RAP";
  datafile << dataPath()
           << "rawdata/" << dirname << "/" << dirname << "_NORM.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Open uncertainty breakdown file
  stringstream sysfile("");
  sysfile << dataPath()
           << "rawdata/" << dirname << "/" << dirname << "_NORM_SYS_BREAKDOWN.data";
  f2.open(sysfile.str().c_str(), ios::in);

  if (f2.fail())
  {
    cerr << "Error opening data file " << sysfile.str() << endl;
    exit(-1);
  }

  // Open correlation matrix file
  stringstream corrfile("");
  string filename3;
  corrfile << dataPath()
           << "rawdata/" << dirname << "/" << dirname << "_NORM.corr";
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

  // Skip over first two lines
  for (int i=0; i<2; i++)
  {
    getline(f1,line);
  }  

  for (int i=0; i<fNData; i++)
  {
    double rap_top; // Rapidity of top quark
    double rap_top_low, rap_top_high; // Limits of bin

    getline(f1,line);
    istringstream lstream(line);

    lstream >> rap_top_low >> rap_top_high;
    lstream >> fData[i];

    // Convert to correct units
    fData[i] /= 1000;

    rap_top = 0.5*(rap_top_low + rap_top_high);
    
    fKin1[i] = rap_top;
    fKin2[i] = Mt*Mt; // Top mass squared
    fKin3[i] = 8000; // Centre of mass energy in GeV

    fStat[i] = 1e-10; // Set stat. error to zero to avoid double counting
  }

  // Read file with systematic uncertainty breakdown
  // Skip over first line
  getline(f2,line);

  double sys1, sys2;
  string unneeded_info;
  const int uncerts=16; // Number of uncertainty sources (inc. stat. errors)

  std::vector<double> fdata_stat(fNData); // Stat. data error
  std::vector<double> fmc_stat(fNData); // Stat. Monte Carlo error

  for (int j=0; j<uncerts; j++)
  {
    getline(f2,line);
    istringstream lstream(line);

    for (int i=0; i<fNData; i++)
    {
      if (j==0) // Read stat. data error
        lstream >> fdata_stat[i] >> unneeded_info;
      else if (j==1) // Read stat. MC error
        lstream >> fmc_stat[i] >> unneeded_info;
      else // Read systematic errors
      {
        lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;

        if (sys1<0. && sys2<0.)
          sys2=0;
        else if (sys1>0. && sys2>0.)
          sys1=0.;

        sys1=sys1/sqrt(2.);
        sys2=sys2/sqrt(2.);

        fSys[i][fNData+2*j-4].mult = sys1;
        fSys[i][fNData+2*j-4].add  = fSys[i][fNData+2*j-4].mult*fData[i]/100;
        fSys[i][fNData+2*j-4].type = MULT;
        fSys[i][fNData+2*j-4].name = "CORR";

        fSys[i][fNData+2*j+1-4].mult = sys2;
        fSys[i][fNData+2*j+1-4].add  = fSys[i][fNData+2*j+1-4].mult*fData[i]/100;
        fSys[i][fNData+2*j+1-4].type = MULT;
        fSys[i][fNData+2*j+1-4].name = "CORR";
      }
    }
  }

  // Compute total additive statistical error
  for (int i=0; i<fNData; i++)
  {
    fstat_additive[i] = (fdata_stat[i] + fmc_stat[i])*fData[i]/100;
  }

  // Read statistical correlation matrix
  // Skip over first line
  getline(f3,line);

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
      fSys[i][j].add  = syscor[i][j];
      fSys[i][j].mult = fSys[i][j].add*100/fData[i];
      fSys[i][j].type = ADD;
      fSys[i][j].name = "CORR";
    }
  }

  // Clean-up
  for (int i=0; i<fNData; i++)
    delete[] syscor[i];

  delete[] syscor;

  for (int i=0; i<fNData; i++)
    delete[] covmat[i];

  delete[] covmat;
}

//==============================================================================

// 2) Normalised distribution differential in modulus of antitop quark rapidity

void ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_RAP_NORMFilter::ReadData()
{
  // Create stream to read data file
  fstream f1, f2, f3;

  // Data files
  stringstream datafile("");
  string dirname;
  dirname = "ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_RAP";
  datafile << dataPath()
           << "rawdata/" << dirname << "/" << dirname << "_NORM.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Open uncertainty breakdown file
  stringstream sysfile("");
  sysfile << dataPath()
           << "rawdata/" << dirname << "/" << dirname << "_NORM_SYS_BREAKDOWN.data";
  f2.open(sysfile.str().c_str(), ios::in);

  if (f2.fail())
  {
    cerr << "Error opening data file " << sysfile.str() << endl;
    exit(-1);
  }

  // Open correlation matrix file
  stringstream corrfile("");
  string filename3;
  corrfile << dataPath()
           << "rawdata/" << dirname << "/" << dirname << "_NORM.corr";
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

  // Skip over first two lines
  for (int i=0; i<2; i++)
  {
    getline(f1,line);
  }  

  for (int i=0; i<fNData; i++)
  {
    double rap_top; // Rapidity of top quark
    double rap_top_low, rap_top_high; // Limits of bin

    getline(f1,line);
    istringstream lstream(line);

    lstream >> rap_top_low >> rap_top_high;
    lstream >> fData[i];

    // Convert to correct units
    fData[i] /= 1000;

    rap_top = 0.5*(rap_top_low + rap_top_high);

    fKin1[i] = rap_top;
    fKin2[i] = Mt*Mt; // Top mass squared
    fKin3[i] = 8000; // Centre of mass energy in GeV

    fStat[i] = 1e-10; // Set stat. error to zero to avoid double counting
  }

  // Read file with systematic uncertainty breakdown
  // Skip over first line
  getline(f2,line);

  double sys1, sys2;
  string unneeded_info;
  const int uncerts=16; // Number of uncertainty sources (inc. stat. errors)

  std::vector<double> fdata_stat(fNData); // Stat. data error
  std::vector<double> fmc_stat(fNData); // Stat. Monte Carlo error

  for (int j=0; j<uncerts; j++)
  {
    getline(f2,line);
    istringstream lstream(line);

    for (int i=0; i<fNData; i++)
    {
      if (j==0) // Read stat. data error
        lstream >> fdata_stat[i] >> unneeded_info;
      else if (j==1) // Read stat. MC error
        lstream >> fmc_stat[i] >> unneeded_info;
      else // Read systematic errors
      {
        lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;

        if (sys1<0. && sys2<0.)
          sys2=0;
        else if (sys1>0. && sys2>0.)
          sys1=0.;

        sys1=sys1/sqrt(2.);
        sys2=sys2/sqrt(2.);

        fSys[i][fNData+2*j-4].mult = sys1;
        fSys[i][fNData+2*j-4].add  = fSys[i][fNData+2*j-4].mult*fData[i]/100;
        fSys[i][fNData+2*j-4].type = MULT;
        fSys[i][fNData+2*j-4].name = "CORR";

        fSys[i][fNData+2*j+1-4].mult = sys2;
        fSys[i][fNData+2*j+1-4].add  = fSys[i][fNData+2*j+1-4].mult*fData[i]/100;
        fSys[i][fNData+2*j+1-4].type = MULT;
        fSys[i][fNData+2*j+1-4].name = "CORR";
      }
    }
  }

  // Compute total additive statistical error
  for (int i=0; i<fNData; i++)
  {
    fstat_additive[i] = (fdata_stat[i] + fmc_stat[i])*fData[i]/100;
  }

  // Read statistical correlation matrix
  // Skip over first line
  getline(f3,line);

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
      fSys[i][j].add  = syscor[i][j];
      fSys[i][j].mult = fSys[i][j].add*100/fData[i];
      fSys[i][j].type = ADD;
      fSys[i][j].name = "CORR";
    }
  }

  // Clean-up
  for (int i=0; i<fNData; i++)
    delete[] syscor[i];

  delete[] syscor;

  for (int i=0; i<fNData; i++)
    delete[] covmat[i];

  delete[] covmat;
}

//==========================================================================

// 3) Normalised distribution differential in top quark transverse momentum

void ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_PT_NORMFilter::ReadData()
{
  // Create stream to read data file
  fstream f1, f2, f3;

  // Data files
  stringstream datafile("");
  string dirname;
  dirname = "ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_PT";
  datafile << dataPath()
           << "rawdata/" << dirname << "/" << dirname << "_NORM.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Open uncertainty breakdown file
  stringstream sysfile("");
  sysfile << dataPath()
           << "rawdata/" << dirname << "/" << dirname << "_NORM_SYS_BREAKDOWN.data";
  f2.open(sysfile.str().c_str(), ios::in);

  if (f2.fail())
  {
    cerr << "Error opening data file " << sysfile.str() << endl;
    exit(-1);
  }

  // Open correlation matrix file
  stringstream corrfile("");
  string filename3;
  corrfile << dataPath()
           << "rawdata/" << dirname << "/" << dirname << "_NORM.corr";
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

  // Skip over first two lines
  for (int i=0; i<2; i++)
  {
    getline(f1,line);
  }  

  for (int i=0; i<fNData; i++)
  {
    double pt_top; // Transverse momentum of top quark
    double pt_top_low, pt_top_high; // Limits of bin

    getline(f1,line);
    istringstream lstream(line);

    lstream >> pt_top_low >> pt_top_high;
    lstream >> fData[i];

    // Convert to 1/GeV
    fData[i] /= 1000;

    pt_top = 0.5*(pt_top_low + pt_top_high);

    fKin1[i] = pt_top;
    fKin2[i] = Mt*Mt; // Top mass squared
    fKin3[i] = 8000; // Centre of mass energy in GeV

    fStat[i] = 1e-10; // Set stat. error to zero to avoid double counting
  }

  // Read file with systematic uncertainty breakdown
  // Skip over first line
  getline(f2,line);

  double sys1, sys2;
  string unneeded_info;
  const int uncerts=17; // Number of uncertainty sources (inc. stat. errors)

  std::vector<double> fdata_stat(fNData); // Stat. data error
  std::vector<double> fmc_stat(fNData); // Stat. Monte Carlo error

  for (int j=0; j<uncerts; j++)
  {
    getline(f2,line);
    istringstream lstream(line);

    for (int i=0; i<fNData; i++)
    {
      if (j==0) // Read stat. data error
        lstream >> fdata_stat[i] >> unneeded_info;
      else if (j==1) // Read stat. MC error
        lstream >> fmc_stat[i] >> unneeded_info;
      else // Read systematic errors
      {
        lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;

        if (sys1<0. && sys2<0.)
          sys2=0;
        else if (sys1>0. && sys2>0.)
          sys1=0.;

        sys1=sys1/sqrt(2.);
        sys2=sys2/sqrt(2.);

        fSys[i][fNData+2*j-4].mult = sys1;
        fSys[i][fNData+2*j-4].add  = fSys[i][fNData+2*j-4].mult*fData[i]/100;
        fSys[i][fNData+2*j-4].type = MULT;
        fSys[i][fNData+2*j-4].name = "CORR";

        fSys[i][fNData+2*j+1-4].mult = sys2;
        fSys[i][fNData+2*j+1-4].add  = fSys[i][fNData+2*j+1-4].mult*fData[i]/100;
        fSys[i][fNData+2*j+1-4].type = MULT;
        fSys[i][fNData+2*j+1-4].name = "CORR";
      }
    }
  }

  // Compute total additive statistical error
  for (int i=0; i<fNData; i++)
  {
    fstat_additive[i] = (fdata_stat[i] + fmc_stat[i])*fData[i]/100;
  }

  // Read statistical correlation matrix
  // Skip over first line
  getline(f3,line);

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
      fSys[i][j].add  = syscor[i][j];
      fSys[i][j].mult = fSys[i][j].add*100/fData[i];
      fSys[i][j].type = ADD;
      fSys[i][j].name = "CORR";
    }
  }

  // Clean-up
  for (int i=0; i<fNData; i++)
    delete[] syscor[i];

  delete[] syscor;

  for (int i=0; i<fNData; i++)
    delete[] covmat[i];

  delete[] covmat;
}

//==============================================================================

// 4) Normalised distribution differential in antitop quark transverse momentum

void ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_PT_NORMFilter::ReadData()
{
  // Create stream to read data file
  fstream f1, f2, f3;

  // Data files
  stringstream datafile("");
  string dirname;
  dirname = "ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_PT";
  datafile << dataPath()
           << "rawdata/" << dirname << "/" << dirname << "_NORM.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Open uncertainty breakdown file
  stringstream sysfile("");
  sysfile << dataPath()
           << "rawdata/" << dirname << "/" << dirname << "_NORM_SYS_BREAKDOWN.data";
  f2.open(sysfile.str().c_str(), ios::in);

  if (f2.fail())
  {
    cerr << "Error opening data file " << sysfile.str() << endl;
    exit(-1);
  }

  // Open correlation matrix file
  stringstream corrfile("");
  string filename3;
  corrfile << dataPath()
           << "rawdata/" << dirname << "/" << dirname << "_NORM.corr";
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

  // Skip over first two lines
  for (int i=0; i<2; i++)
  {
    getline(f1,line);
  }  

  for (int i=0; i<fNData; i++)
  {
    double pt_top; // Transverse momentum of top quark
    double pt_top_low, pt_top_high; // Limits of bin

    getline(f1,line);
    istringstream lstream(line);

    lstream >> pt_top_low >> pt_top_high;
    lstream >> fData[i];

    // Convert to 1/GeV
    fData[i] /= 1000;

    pt_top = 0.5*(pt_top_low + pt_top_high);

    fKin1[i] = pt_top;
    fKin2[i] = Mt*Mt; // Top mass squared
    fKin3[i] = 8000; // Centre of mass energy in GeV

    fStat[i] = 1e-10; // Set stat. error to zero to avoid double counting
  }

  // Read file with systematic uncertainty breakdown
  // Skip over first line
  getline(f2,line);

  double sys1, sys2;
  string unneeded_info;
  const int uncerts=17; // Number of uncertainty sources (inc. stat. errors)

  std::vector<double> fdata_stat(fNData); // Stat. data error
  std::vector<double> fmc_stat(fNData); // Stat. Monte Carlo error

  for (int j=0; j<uncerts; j++)
  {
    getline(f2,line);
    istringstream lstream(line);

    for (int i=0; i<fNData; i++)
    {
      if (j==0) // Read stat. data error
        lstream >> fdata_stat[i] >> unneeded_info;
      else if (j==1) // Read stat. MC error
        lstream >> fmc_stat[i] >> unneeded_info;
      else // Read systematic errors
      {
        lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;

        if (sys1<0. && sys2<0.)
          sys2=0;
        else if (sys1>0. && sys2>0.)
          sys1=0.;

        sys1=sys1/sqrt(2.);
        sys2=sys2/sqrt(2.);

        fSys[i][fNData+2*j-4].mult = sys1;
        fSys[i][fNData+2*j-4].add  = fSys[i][fNData+2*j-4].mult*fData[i]/100;
        fSys[i][fNData+2*j-4].type = MULT;
        fSys[i][fNData+2*j-4].name = "CORR";

        fSys[i][fNData+2*j+1-4].mult = sys2;
        fSys[i][fNData+2*j+1-4].add  = fSys[i][fNData+2*j+1-4].mult*fData[i]/100;
        fSys[i][fNData+2*j+1-4].type = MULT;
        fSys[i][fNData+2*j+1-4].name = "CORR";
      }
    }
  }

  // Compute total additive statistical error
  for (int i=0; i<fNData; i++)
  {
    fstat_additive[i] = (fdata_stat[i] + fmc_stat[i])*fData[i]/100;
  }

  // Read statistical correlation matrix
  // Skip over first line
  getline(f3,line);

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
      fSys[i][j].add  = syscor[i][j];
      fSys[i][j].mult = fSys[i][j].add*100/fData[i];
      fSys[i][j].type = ADD;
      fSys[i][j].name = "CORR";
    }
  }

  // Clean-up
  for (int i=0; i<fNData; i++)
    delete[] syscor[i];

  delete[] syscor;

  for (int i=0; i<fNData; i++)
    delete[] covmat[i];

  delete[] covmat;
}

//==================================================================

// B - UNNORMALISED distributions

// 5) Distribution differential in modulus of top quark rapidity

void ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_RAPFilter::ReadData()
{
  // Create stream to read data file
  fstream f1, f2, f3;

  // Data files
  stringstream datafile("");
  string dirname;
  dirname = "ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_RAP";
  datafile << dataPath()
           << "rawdata/" << dirname << "/" << dirname << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Open uncertainty breakdown file
  stringstream sysfile("");
  sysfile << dataPath()
           << "rawdata/" << dirname << "/" << dirname << "_SYS_BREAKDOWN.data";
  f2.open(sysfile.str().c_str(), ios::in);

  if (f2.fail())
  {
    cerr << "Error opening data file " << sysfile.str() << endl;
    exit(-1);
  }

  // Open correlation matrix file
  stringstream corrfile("");
  string filename3;
  corrfile << dataPath()
           << "rawdata/" << dirname << "/" << dirname << ".corr";
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

  // Skip over first two lines
  for (int i=0; i<2; i++)
  {
    getline(f1,line);
  }  

  for (int i=0; i<fNData; i++)
  {
    double rap_top; // Rapidity of top quark
    double rap_top_low, rap_top_high; // Limits of bin

    getline(f1,line);
    istringstream lstream(line);

    lstream >> rap_top_low >> rap_top_high;
    lstream >> fData[i];

    rap_top = 0.5*(rap_top_low + rap_top_high);

    fKin1[i] = rap_top;
    fKin2[i] = Mt*Mt; // Top mass squared
    fKin3[i] = 8000; // Centre of mass energy in GeV

    fStat[i] = 1e-10; // Set stat. error to zero to avoid double counting
  }

  // Read file with systematic uncertainty breakdown
  // Skip over first line
  getline(f2,line);

  double sys1, sys2;
  double sys;
  string unneeded_info;
  const int uncerts=18; // Number of uncertainty sources (inc. stat. errors)

  std::vector<double> fdata_stat(fNData); // Stat. data error
  std::vector<double> fmc_stat(fNData); // Stat. Monte Carlo error

  for (int j=0; j<uncerts; j++)
  {
    getline(f2,line);
    istringstream lstream(line);

    for (int i=0; i<fNData; i++)
    {
      if (j==0) // Read stat. data error
        lstream >> fdata_stat[i] >> unneeded_info;
      else if (j==1) // Read stat. MC error
        lstream >> fmc_stat[i] >> unneeded_info;
      else if (j==uncerts-1) // Read luminosity uncertainty
      {
        lstream >> sys >> unneeded_info;

        fSys[i][fNData+2*j-4].mult = sys;
        fSys[i][fNData+2*j-4].add  = fSys[i][fNData+2*j-4].mult*fData[i]/100;
        fSys[i][fNData+2*j-4].type = MULT;
        fSys[i][fNData+2*j-4].name = "ATLASLUMI12";
      }
      else // Read systematic errors
      {
        lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;

        if (sys1<0. && sys2<0.)
          sys2=0;
        else if (sys1>0. && sys2>0.)
          sys1=0.;

        sys1=sys1/sqrt(2.);
        sys2=sys2/sqrt(2.);

        fSys[i][fNData+2*j-4].mult = sys1;
        fSys[i][fNData+2*j-4].add  = fSys[i][fNData+2*j-4].mult*fData[i]/100;
        fSys[i][fNData+2*j-4].type = MULT;
        fSys[i][fNData+2*j-4].name = "CORR";

        fSys[i][fNData+2*j+1-4].mult = sys2;
        fSys[i][fNData+2*j+1-4].add  = fSys[i][fNData+2*j+1-4].mult*fData[i]/100;
        fSys[i][fNData+2*j+1-4].type = MULT;
        fSys[i][fNData+2*j+1-4].name = "CORR";
      }
    }
  }

  // Compute total additive statistical error
  for (int i=0; i<fNData; i++)
  {
    fstat_additive[i] = (fdata_stat[i] + fmc_stat[i])*fData[i]/100;
  }

  // Read statistical correlation matrix
  // Skip over first line
  getline(f3,line);

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
      fSys[i][j].add  = syscor[i][j];
      fSys[i][j].mult = fSys[i][j].add*100/fData[i];
      fSys[i][j].type = ADD;
      fSys[i][j].name = "CORR";
    }
  }

  // Clean-up
  for (int i=0; i<fNData; i++)
    delete[] syscor[i];

  delete[] syscor;

  for (int i=0; i<fNData; i++)
    delete[] covmat[i];

  delete[] covmat;
}

//==================================================================

// 6) Distribution differential in modulus of antitop quark rapidity

void ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_RAPFilter::ReadData()
{
  // Create stream to read data file
  fstream f1, f2, f3;

  // Data files
  stringstream datafile("");
  string dirname;
  dirname = "ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_RAP";
  datafile << dataPath()
           << "rawdata/" << dirname << "/" << dirname << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Open uncertainty breakdown file
  stringstream sysfile("");
  sysfile << dataPath()
           << "rawdata/" << dirname << "/" << dirname << "_SYS_BREAKDOWN.data";
  f2.open(sysfile.str().c_str(), ios::in);

  if (f2.fail())
  {
    cerr << "Error opening data file " << sysfile.str() << endl;
    exit(-1);
  }

  // Open correlation matrix file
  stringstream corrfile("");
  string filename3;
  corrfile << dataPath()
           << "rawdata/" << dirname << "/" << dirname << ".corr";
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

  // Skip over first two lines
  for (int i=0; i<2; i++)
  {
    getline(f1,line);
  }  

  for (int i=0; i<fNData; i++)
  {
    double rap_top; // Rapidity of top quark
    double rap_top_low, rap_top_high; // Limits of bin

    getline(f1,line);
    istringstream lstream(line);

    lstream >> rap_top_low >> rap_top_high;
    lstream >> fData[i];

    rap_top = 0.5*(rap_top_low + rap_top_high);

    fKin1[i] = rap_top;
    fKin2[i] = Mt*Mt; // Top mass squared
    fKin3[i] = 8000; // Centre of mass energy in GeV

    fStat[i] = 1e-10; // Set stat. error to zero to avoid double counting
  }

  // Read file with systematic uncertainty breakdown
  // Skip over first line
  getline(f2,line);

  double sys1, sys2;
  double sys;
  string unneeded_info;
  const int uncerts=19; // Number of uncertainty sources (inc. stat. errors)

  std::vector<double> fdata_stat(fNData); // Stat. data error
  std::vector<double> fmc_stat(fNData); // Stat. Monte Carlo error

  for (int j=0; j<uncerts; j++)
  {
    getline(f2,line);
    istringstream lstream(line);

    for (int i=0; i<fNData; i++)
    {
      if (j==0) // Read stat. data error
        lstream >> fdata_stat[i] >> unneeded_info;
      else if (j==1) // Read stat. MC error
        lstream >> fmc_stat[i] >> unneeded_info;
      else if (j==uncerts-1) // Read luminosity uncertainty
      {
        lstream >> sys >> unneeded_info;

        fSys[i][fNData+2*j-4].mult = sys;
        fSys[i][fNData+2*j-4].add  = fSys[i][fNData+2*j-4].mult*fData[i]/100;
        fSys[i][fNData+2*j-4].type = MULT;
        fSys[i][fNData+2*j-4].name = "ATLASLUMI12";
      }
      else // Read systematic errors
      {
        lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;

        if (sys1<0. && sys2<0.)
          sys2=0;
        else if (sys1>0. && sys2>0.)
          sys1=0.;

        sys1=sys1/sqrt(2.);
        sys2=sys2/sqrt(2.);

        fSys[i][fNData+2*j-4].mult = sys1;
        fSys[i][fNData+2*j-4].add  = fSys[i][fNData+2*j-4].mult*fData[i]/100;
        fSys[i][fNData+2*j-4].type = MULT;
        fSys[i][fNData+2*j-4].name = "CORR";

        fSys[i][fNData+2*j+1-4].mult = sys2;
        fSys[i][fNData+2*j+1-4].add  = fSys[i][fNData+2*j+1-4].mult*fData[i]/100;
        fSys[i][fNData+2*j+1-4].type = MULT;
        fSys[i][fNData+2*j+1-4].name = "CORR";
      }
    }
  }

  // Compute total additive statistical error
  for (int i=0; i<fNData; i++)
  {
    fstat_additive[i] = (fdata_stat[i] + fmc_stat[i])*fData[i]/100;
  }

  // Read statistical correlation matrix
  // Skip over first line
  getline(f3,line);

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
      fSys[i][j].add  = syscor[i][j];
      fSys[i][j].mult = fSys[i][j].add*100/fData[i];
      fSys[i][j].type = ADD;
      fSys[i][j].name = "CORR";
    }
  }

  // Clean-up
  for (int i=0; i<fNData; i++)
    delete[] syscor[i];

  delete[] syscor;

  for (int i=0; i<fNData; i++)
    delete[] covmat[i];

  delete[] covmat;
}

//==================================================================

// 7) Distribution differential in top quark transverse momentum

void ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_PTFilter::ReadData()
{
  // Create stream to read data file
  fstream f1, f2, f3;

  // Data files
  stringstream datafile("");
  string dirname;
  dirname = "ATLAS_SINGLETOP_TCH_DIFF_8TEV_T_PT";
  datafile << dataPath()
           << "rawdata/" << dirname << "/" << dirname << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Open uncertainty breakdown file
  stringstream sysfile("");
  sysfile << dataPath()
           << "rawdata/" << dirname << "/" << dirname << "_SYS_BREAKDOWN.data";
  f2.open(sysfile.str().c_str(), ios::in);

  if (f2.fail())
  {
    cerr << "Error opening data file " << sysfile.str() << endl;
    exit(-1);
  }

  // Open correlation matrix file
  stringstream corrfile("");
  string filename3;
  corrfile << dataPath()
           << "rawdata/" << dirname << "/" << dirname << ".corr";
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

  // Skip over first two lines
  for (int i=0; i<2; i++)
  {
    getline(f1,line);
  }  

  for (int i=0; i<fNData; i++)
  {
    double pt_top; // Transverse momentum of top quark
    double pt_top_low, pt_top_high; // Limits of bin

    getline(f1,line);
    istringstream lstream(line);

    lstream >> pt_top_low >> pt_top_high;
    lstream >> fData[i];

    // Convert to pb/GeV
    fData[i] /= 1000;

    pt_top = 0.5*(pt_top_low + pt_top_high);

    fKin1[i] = pt_top;
    fKin2[i] = Mt*Mt; // Top mass squared
    fKin3[i] = 8000; // Centre of mass energy in GeV

    fStat[i] = 1e-10; // Set stat. error to zero to avoid double counting
  }

  // Read file with systematic uncertainty breakdown
  // Skip over first line
  getline(f2,line);

  double sys1, sys2;
  double sys;
  string unneeded_info;
  const int uncerts=19; // Number of uncertainty sources (inc. stat. errors)

  std::vector<double> fdata_stat(fNData); // Stat. data error
  std::vector<double> fmc_stat(fNData); // Stat. Monte Carlo error

  for (int j=0; j<uncerts; j++)
  {
    getline(f2,line);
    istringstream lstream(line);

    for (int i=0; i<fNData; i++)
    {
      if (j==0) // Read stat. data error
        lstream >> fdata_stat[i] >> unneeded_info;
      else if (j==1) // Read stat. MC error
        lstream >> fmc_stat[i] >> unneeded_info;
      else if (j==uncerts-1) // Read luminosity uncertainty
      {
        lstream >> sys >> unneeded_info;

        fSys[i][fNData+2*j-4].mult = sys;
        fSys[i][fNData+2*j-4].add  = fSys[i][fNData+2*j-4].mult*fData[i]/100;
        fSys[i][fNData+2*j-4].type = MULT;
        fSys[i][fNData+2*j-4].name = "ATLASLUMI12";
      }
      else // Read systematic errors
      {
        lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;

        if (sys1<0. && sys2<0.)
          sys2=0;
        else if (sys1>0. && sys2>0.)
          sys1=0.;

        sys1=sys1/sqrt(2.);
        sys2=sys2/sqrt(2.);

        fSys[i][fNData+2*j-4].mult = sys1;
        fSys[i][fNData+2*j-4].add  = fSys[i][fNData+2*j-4].mult*fData[i]/100;
        fSys[i][fNData+2*j-4].type = MULT;
        fSys[i][fNData+2*j-4].name = "CORR";

        fSys[i][fNData+2*j+1-4].mult = sys2;
        fSys[i][fNData+2*j+1-4].add  = fSys[i][fNData+2*j+1-4].mult*fData[i]/100;
        fSys[i][fNData+2*j+1-4].type = MULT;
        fSys[i][fNData+2*j+1-4].name = "CORR";
      }
    }
  }

  // Compute total additive statistical error
  for (int i=0; i<fNData; i++)
  {
    fstat_additive[i] = (fdata_stat[i] + fmc_stat[i])*fData[i]/100;
  }

  // Read statistical correlation matrix
  // Skip over first line
  getline(f3,line);

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
      fSys[i][j].add  = syscor[i][j];
      fSys[i][j].mult = fSys[i][j].add*100/fData[i];
      fSys[i][j].type = ADD;
      fSys[i][j].name = "CORR";
    }
  }

  // Clean-up
  for (int i=0; i<fNData; i++)
    delete[] syscor[i];

  delete[] syscor;

  for (int i=0; i<fNData; i++)
    delete[] covmat[i];

  delete[] covmat;
}

//==================================================================

// 8) Distribution differential in antitop quark transverse momentum

void ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_PTFilter::ReadData()
{
  // Create stream to read data file
  fstream f1, f2, f3;

  // Data files
  stringstream datafile("");
  string dirname;
  dirname = "ATLAS_SINGLETOP_TCH_DIFF_8TEV_TBAR_PT";
  datafile << dataPath()
           << "rawdata/" << dirname << "/" << dirname << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Open uncertainty breakdown file
  stringstream sysfile("");
  sysfile << dataPath()
           << "rawdata/" << dirname << "/" << dirname << "_SYS_BREAKDOWN.data";
  f2.open(sysfile.str().c_str(), ios::in);

  if (f2.fail())
  {
    cerr << "Error opening data file " << sysfile.str() << endl;
    exit(-1);
  }

  // Open correlation matrix file
  stringstream corrfile("");
  string filename3;
  corrfile << dataPath()
           << "rawdata/" << dirname << "/" << dirname << ".corr";
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

  // Skip over first two lines
  for (int i=0; i<2; i++)
  {
    getline(f1,line);
  }  

  for (int i=0; i<fNData; i++)
  {
    double pt_top; // Transverse momentum of top quark
    double pt_top_low, pt_top_high; // Limits of bin

    getline(f1,line);
    istringstream lstream(line);

    lstream >> pt_top_low >> pt_top_high;
    lstream >> fData[i];

    // Convert to pb/GeV
    fData[i] /= 1000;

    pt_top = 0.5*(pt_top_low + pt_top_high);

    fKin1[i] = pt_top;
    fKin2[i] = Mt*Mt; // Top mass squared
    fKin3[i] = 8000; // Centre of mass energy in GeV

    fStat[i] = 1e-10; // Set stat. error to zero to avoid double counting
  }

  // Read file with systematic uncertainty breakdown
  // Skip over first line
  getline(f2,line);

  double sys1, sys2;
  double sys;
  string unneeded_info;
  const int uncerts=19; // Number of uncertainty sources (inc. stat. errors)

  std::vector<double> fdata_stat(fNData); // Stat. data error
  std::vector<double> fmc_stat(fNData); // Stat. Monte Carlo error

  for (int j=0; j<uncerts; j++)
  {
    getline(f2,line);
    istringstream lstream(line);

    for (int i=0; i<fNData; i++)
    {
      if (j==0) // Read stat. data error
        lstream >> fdata_stat[i] >> unneeded_info;
      else if (j==1) // Read stat. MC error
        lstream >> fmc_stat[i] >> unneeded_info;
      else if (j==uncerts-1) // Read luminosity uncertainty
      {
        lstream >> sys >> unneeded_info;

        fSys[i][fNData+2*j-4].mult = sys;
        fSys[i][fNData+2*j-4].add  = fSys[i][fNData+2*j-4].mult*fData[i]/100;
        fSys[i][fNData+2*j-4].type = MULT;
        fSys[i][fNData+2*j-4].name = "ATLASLUMI12";
      }
      else // Read systematic errors
      {
        lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;

        if (sys1<0. && sys2<0.)
          sys2=0;
        else if (sys1>0. && sys2>0.)
          sys1=0.;

        sys1=sys1/sqrt(2.);
        sys2=sys2/sqrt(2.);

        fSys[i][fNData+2*j-4].mult = sys1;
        fSys[i][fNData+2*j-4].add  = fSys[i][fNData+2*j-4].mult*fData[i]/100;
        fSys[i][fNData+2*j-4].type = MULT;
        fSys[i][fNData+2*j-4].name = "CORR";

        fSys[i][fNData+2*j+1-4].mult = sys2;
        fSys[i][fNData+2*j+1-4].add  = fSys[i][fNData+2*j+1-4].mult*fData[i]/100;
        fSys[i][fNData+2*j+1-4].type = MULT;
        fSys[i][fNData+2*j+1-4].name = "CORR";
      }
    }
  }

  // Compute total additive statistical error
  for (int i=0; i<fNData; i++)
  {
    fstat_additive[i] = (fdata_stat[i] + fmc_stat[i])*fData[i]/100;
  }

  // Read statistical correlation matrix
  // Skip over first line
  getline(f3,line);

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
      fSys[i][j].add  = syscor[i][j];
      fSys[i][j].mult = fSys[i][j].add*100/fData[i];
      fSys[i][j].type = ADD;
      fSys[i][j].name = "CORR";
    }
  }

  // Clean-up
  for (int i=0; i<fNData; i++)
    delete[] syscor[i];

  delete[] syscor;

  for (int i=0; i<fNData; i++)
    delete[] covmat[i];

  delete[] covmat;
}
