#include "ATLAS_SINGLETOP_TCH_DIFF_7TEV.h"

void ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_RAP_NORMFilter::ReadData()
{
  // Create streams to read data files
  fstream f1, f2, f3;

  // Data files
  stringstream datafile1("");
  string filename1;
  filename1 = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_RAP";
  datafile1 << dataPath()
           << "rawdata/" << filename1 << "/" << filename1 << ".data";
  f1.open(datafile1.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile1.str() << endl;
    exit(-1);
  }

  stringstream datafile2("");
  string filename2;
  filename2 = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_RAP_NORM_SYS_BREAKDOWN";
  datafile2 << dataPath()
           << "rawdata/" << filename1 << "/" << filename2 << ".data";
  f2.open(datafile2.str().c_str(), ios::in);

  if (f2.fail())
  {
    cerr << "Error opening data file " << datafile2.str() << endl;
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
// double fstat_additive[4] = {0, 0, 0, 0};

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

    lstream >> rap_top_high >> unneeded_info >> rap_top_low;
    rap_top = 0.5*(rap_top_high + rap_top_low);

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
      cout << corrmat[i][j] << "\t";
      cout << covmat[i][j] << "\t";
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

  // Artificial systematics
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

  // Skip over first 20 lines (including stat. uncert.)
  for (int i=0; i<20; i++)
  {
    getline(f2,line);
  }
  
  double sys1, sys2, up, down, sigma, datshift;
  double shift[4] = {0, 0, 0, 0};

  for (int i=fNData; i<fNSys; i++)
  {
    string unneeded_info;

    getline(f2,line);
    istringstream lstream(line);

    if (i==fNData+2) // Deal with aysymmetric errors
    {
      for (int j=0; j<3; j++)
      {
        lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
        if (sys1 < 0) {up=sys2; down=sys1;}
        else {up=sys1; down=sys2;}
        symmetriseErrors(up, down, &sigma, &datshift);
        fSys[j][i].mult = sigma;
        fSys[j][i].add = fSys[j][i].mult*fData[j]/100;
        fSys[j][i].type = MULT;
        fSys[j][i].name = "CORR";    
        shift[j] += datshift; 
      }
  
      lstream >> fSys[3][i].mult;
      fSys[3][i].add = fSys[3][i].mult*fData[3]/100;
      fSys[3][i].type = MULT;
      fSys[3][i].name = "CORR";
    }
    else if (i==fNData+3) // Deal with asymmetric errors
    {
      for (int j=0; j<2; j++)
      {
        lstream >> fSys[j][i].mult >> unneeded_info;
        fSys[j][i].add = fSys[j][i].mult*fData[j]/100;
        fSys[j][i].type = MULT;
        fSys[j][i].name = "CORR";
      }

      lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
      if (sys1 < 0) {up=sys2; down=sys1;}
      else {up=sys1; down=sys2;}
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[2][i].mult = sigma;
      fSys[2][i].add = fSys[2][i].mult*fData[2]/100;
      fSys[2][i].type = MULT;
      fSys[2][i].name = "CORR";
      shift[2] += datshift;

      lstream >> fSys[3][i].mult;
      fSys[3][i].add = fSys[3][i].mult*fData[3]/100;
      fSys[3][i].type = MULT;
      fSys[3][i].name = "CORR";
    }
    else if (i==fNData+5) // Deal with asymmetric errors
    {
      for (int j=0; j<2; j++)
      {
        lstream >> sys1 >> unneeded_info >> sys2 >> unneeded_info;
        if (sys1 < 0) {up=sys2; down=sys1;}
        else {up=sys1; down=sys2;}
        symmetriseErrors(up, down, &sigma, &datshift);
        fSys[j][i].mult = sigma;
        fSys[j][i].add = fSys[j][i].mult*fData[j]/100;
        fSys[j][i].type = MULT;
        fSys[j][i].name = "CORR";
        shift[j] += datshift;
      }

      lstream >> fSys[2][i].mult >> unneeded_info;
      fSys[2][i].add = fSys[2][i].mult*fData[2]/100;
      fSys[2][i].type = MULT;
      fSys[2][i].name = "CORR";

      lstream >> sys1 >> unneeded_info >> sys2;
      if (sys1 < 0) {up=sys2; down=sys1;}
      else {up=sys1; down=sys2;}
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[3][i].mult = sigma;
      fSys[3][i].add = fSys[3][i].mult*fData[3]/100;
      fSys[3][i].type = MULT;
      fSys[3][i].name = "CORR";
      shift[3] += datshift;
    }
    else
    {
      lstream >> fSys[0][i].mult >> unneeded_info >> fSys[1][i].mult >> unneeded_info >> fSys[2][i].mult >> unneeded_info >> fSys[3][i].mult;
      for (int j=0; j<4; j++)
      {
        fSys[j][i].add = fSys[j][i].mult*fData[j]/100;
        fSys[j][i].type = MULT;
        fSys[j][i].name = "CORR";
      }
    }
  }

  for (int i=0; i<4; i++)
  {
    fData[i] *= (1.0 + shift[i]*0.01); // Shift of central value due to asymmetric errors
  }
}
