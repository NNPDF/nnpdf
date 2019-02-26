#include "ATLAS_SINGLETOP_TCH_DIFF_7TEV.h"

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
  double fstat_additive[4];

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

    fStat[i] = 0; // Set stat. error to zero to avoid double counting when using artificial systematics
  }

  // Read statistical correlation matrix
  // Skip over first ten lines
  for (int i=0; i<10; i++)
  {
    getline(f3,line);
  }

  double** covmat = new double*[fNData];
  double corrmat[fNData][fNData];
  for (int i=0; i<fNData; i++)
  {
    string unneeded_info;
    covmat[i] = new double[fNData];
    getline(f3,line);
    istringstream lstream(line);
    lstream >> unneeded_info >> unneeded_info >> unneeded_info;
    for (int j=0; j<fNData; j++)
    {
      lstream >> corrmat[i][j] >> unneeded_info;
      covmat[i][j] = corrmat[i][j] * fstat_additive[i] * fstat_additive[j];
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
  double shift[4] = {0, 0, 0, 0};

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
        if (sys1 < 0) {up=sys2; down=sys1;}
        else {up=sys1; down=sys2;}
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
      if (sys1 < 0) {up=sys2; down=sys1;}
      else {up=sys1; down=sys2;}
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
        if (sys1 < 0) {up=sys2; down=sys1;}
        else {up=sys1; down=sys2;}
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
      if (sys1 < 0) {up=sys2; down=sys1;}
      else {up=sys1; down=sys2;}
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[3][j].mult = sigma;
      fSys[3][j].add = fSys[3][j].mult*fData[3]/100;
      fSys[3][j].type = MULT;
      fSys[3][j].name = "CORR";
      shift[3] += datshift;
    }
    else // Deal with lines that contain no asymmetric errors
    {
      lstream >> fSys[0][j].mult >> unneeded_info >> fSys[1][j].mult >> unneeded_info >> fSys[2][j].mult >> unneeded_info >> fSys[3][j].mult;
      for (int i=0; i<4; i++)
      {
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }
    }
  }

  for (int i=0; i<4; i++)
  {
    fData[i] *= (1.0 + shift[i]*0.01); // Shift of central value due to asymmetric errors
  }
}

//==================================================================

// 2) Distribution differential in modulus of antitop quark rapidity
void ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_RAP_NORMFilter::ReadData()
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
  filename2 = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_RAP_NORM_SYS_BREAKDOWN";
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
  filename3 = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_TBAR_RAP_NORM";
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
  double fstat_additive[4];

  // Skip over first ten lines
  for (int i=0; i<10; i++)
  {
    getline(f1,line);
  }
  
  for (int i=0; i<fNData; i++)
  {
    double rap_top; // Rapidity of top quark
    double rap_top_high, rap_top_low; // Limits of bin
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

    fStat[i] = 0; // Set stat. error to zero to avoid double counting when using artificial systematics
  }

  // Read statistical correlation matrix
  // Skip over first ten lines
  for (int i=0; i<10; i++)
  {
    getline(f3,line);
  }

  double** covmat = new double*[fNData];
  double corrmat[fNData][fNData];
  for (int i=0; i<fNData; i++)
  {
    string unneeded_info;
    covmat[i] = new double[fNData];
    getline(f3,line);
    istringstream lstream(line);
    lstream >> unneeded_info >> unneeded_info >> unneeded_info;
    for (int j=0; j<fNData; j++)
    {
      lstream >> corrmat[i][j] >> unneeded_info;
      covmat[i][j] = corrmat[i][j] * fstat_additive[i] * fstat_additive[j];
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
  double shift[4] = {0, 0, 0, 0};

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
      if (sys1 < 0) {up=sys2; down=sys1;}
      else {up=sys1; down=sys2;}
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
      if (sys1 < 0) {up=sys2; down=sys1;}
      else {up=sys1; down=sys2;}
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[3][j].mult = sigma;
      fSys[3][j].add = fSys[3][j].mult*fData[3]/100;
      fSys[3][j].type = MULT;
      fSys[3][j].name = "CORR";    
      shift[3] += datshift;
    }
    else if (j==fNData+4) // Deal with asymmetric errors
    {
      for (int i=0; i<4; i++)
      {
        lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
        if (sys1 < 0) {up=sys2; down=sys1;}
        else {up=sys1; down=sys2;}
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
      lstream >> fSys[0][j].mult >> unneeded_info >> fSys[1][j].mult >> unneeded_info >> fSys[2][j].mult >> unneeded_info >> fSys[3][j].mult;
      for (int i=0; i<4; i++)
      {
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }
    }
  }

  for (int i=0; i<4; i++)
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
  double fstat_additive[fNData];

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
    lstream >> unneeded_info >> fstat_percentage;
    fstat_additive[i] = fstat_percentage*fData[i]/100;

    fStat[i] = 0; // Set stat. error to zero to avoid double counting when using artificial systematics
  }

  // Read statistical correlation matrix
  // Skip over first ten lines
  for (int i=0; i<10; i++)
  {
    getline(f3,line);
  }

  double** covmat = new double*[fNData];
  double corrmat[fNData][fNData];
  for (int i=0; i<fNData; i++)
  {
    string unneeded_info;
    covmat[i] = new double[fNData];
    getline(f3,line);
    istringstream lstream(line);
    lstream >> unneeded_info >> unneeded_info >> unneeded_info;
    for (int j=0; j<fNData; j++)
    {
      lstream >> corrmat[i][j] >> unneeded_info;
      covmat[i][j] = corrmat[i][j] * fstat_additive[i] * fstat_additive[j];
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
  double shift[fNData];

  for (int j=fNData; j<fNSys; j++)
  {
    string unneeded_info;

    getline(f2,line);
    istringstream lstream(line);

    if (j==fNData+2) // Deal with aysymmetric errors
    {
    }
    else if (j==fNData+3) // Deal with asymmetric errors
    {
    }
    else if (j==fNData+4) // Deal with asymmetric errors
    {
    }
    else if (j==fNData+5) // Deal with asymmetric errors
    {
    }
    else if (j==fNData+7) // Deal with asymmetric errors
    {
    }
    else if (j==fNData+8) // Deal with asymmetric errors
    {
    }
    else // Deal with lines that contain no asymmetric errors
    {
      lstream >> fSys[0][j].mult >> unneeded_info >> fSys[1][j].mult >> unneeded_info >> fSys[2][j].mult >> unneeded_info >> fSys[3][j].mult;
      for (int i=0; i<4; i++)
      {
        fSys[i][j].add = fSys[i][j].mult*fData[i]/100;
        fSys[i][j].type = MULT;
        fSys[i][j].name = "CORR";
      }
    }
  }

  for (int i=0; i<4; i++)
  {
    fData[i] *= (1.0 + shift[i]*0.01); // Shift of central value due to asymmetric errors
  }
}
