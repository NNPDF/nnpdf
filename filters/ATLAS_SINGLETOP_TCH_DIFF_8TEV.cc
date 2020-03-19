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
    double sigma, datshift; // Outputs of symmetriseErrors

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

    symmetriseErrors(sys1, sys2, &sigma, &datshift);

    fData[i] += datshift; // Shift of central value due to asymmetric errors

    fSys[i][0].add = sigma;
    fSys[i][0].mult = (fSys[i][0].add/fData[i])*100;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CORR";
  }
}
