/*
 * Differential distribution measurements in p_T and abs(y) of the top/antitop quark for single top production data from ATLAS at sqrt(s) = 7 TeV
 * Normalised and unnormalised distributions
 * Measurements are taken from a data set with an integrated luminosity of 4.59 1/fb
 *
 * Archived as: 1406.7844v2
 * DOI: https://doi.org/10.1103/PhysRevD.90.112006
 *
 * Data and correlation information taken from supplemental material of paper: https://journals.aps.org/prd/supplemental/10.1103/PhysRevD.90.112006
 */

#include "ATLAS_SINGLETOP_TCH_DIFF_7TEV.h"

void ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_RAP_NORMFilter::ReadData()
{
  // Create streams to read data files
  fstream f1, f2, f3;

  // Open data files
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

  if (f1.fail())
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
  f2.open(corrfile.str().c_str(), ios::in);

  if (f3.fail())
  {
    cerr << "Error opening data file " << corrfile.str() << endl;
    exit(-1);
  }

  // Read statistical correlation matrix
  // NEED TO CONVERT TO COVARIANCE MATRIX
//  string line;
  
  // Skip over first ten lines
//  for (int i=0; i<10; i++)
//  {
//    getline(f3,line);
//  }

//  double** covmat = new double*[fNData];
//  double corrmat_element;
//  for (int i=0; i<fNData; i++)
//  {
//    string unneeded_info;
//    covmat[i] = new double[fNData];
//    getline(f3,line);
//    istringstream lstream(line);
//    lstream >> unneeded_info >> unneeded_info >> unneeded_info;
//    for (int j=0; j<fNData; j++)
//    {
//      lstream >> corrmat_element >> unneeded_info;
//      covmat[i][j] = corrmat_element * fstat_additive[i] * fstat_additive[j];
//    }
//  }

  // Generate artificial systematics
//  double** syscor = new double*[fNData];
//  for (int i = 0; i < fNData; i++)
//    syscor[i] = new double[fNData];
//
//  if (!genArtSys(fNData,covmat,syscor))
//  {
//    cerr << " in " << fSetName << endl;
//    exit(-1);
//  }

  // Start filter of data
  string line;
  double fstat_additive[4] = {0, 0, 0, 0};

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
    fstat_additive[i] = fstat_percentage*fData[j]/100;

//    // Artificial systematics
//    for (int j=0; j<fNData; j++)
//    {
//      fSys[i][j].add = syscor[i][j];
//      fSys[i][j].mult = fSys[i][j].add*100/fData[i];
//      fSys[i][j].type = "MULT";
//      fSys[i][j].name = "CORR";
//    }
  }

  // Read statistical correlation matrix
  // NEED TO CONVERT TO COVARIANCE MATRIX
  string line; 

  // Skip over first ten lines
  for (int i=0; i<10; i++)
  {
    getline(f3,line);
  }

  double** covmat = new double*[fNData];
  double corrmat_element;
  for (int i=0; i<fNData; i++)
  {
    string unneeded_info;
    covmat[i] = new double[fNData];
    getline(f3,line);
    istringstream lstream(line);
    lstream >> unneeded_info >> unneeded_info >> unneeded_info;
    for (int j=0; j<fNData; j++)
    {
      lstream >> corrmat_element >> unneeded_info;
      covmat[i][j] = corrmat_element * fstat_additive[i] * fstat_additive[j];
    }

  // Generate artificial systematics
  double** syscor = new double*[fNData];
  for (int i = 0; i < fNData; i++)
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
        fSys[i][j].type = "MULT";
        fSys[i][j].name = "CORR";
      }
  }

  // Set statistical uncertainty to zero because we use artificial systematics
  // Streamline by setting in above loop later
  for (int i=0; i<4; i++)
  {
    fStat[i] = 0;
  }

  // Skip over first 20 lines (including stat. uncert.)
  for (int i=0; i<20; i++)
  {
    getline(f2,line);
  }
  
//  int fnuncerts;
//  fnuncerts = 7;
  double up, down, sigma, datshift;
  double shift[4] = {0, 0, 0, 0};

//  shift = 0;
 
  for (int i=fNData; i<fNSys; i++)
  {
//    double fstat_percentage1, fstat_percentage2, fstat_percentage3, fstat_percentage4; // Percentage statistical uncertainties for each bin

//    lstream >> fstat_percentage1 >> unneeded_info >> fstat_percentage2 >> unneeded_info >> fstat_percentage3 >> unneeded_info >> fstat_percentage4;

    string unneeded_info;

    getline(f2,line);
    istringstream lstream(line);

    if (i==(fNData+2))
    {
      for (int j=0; j<3; j++)
      {
        lstream >> up >> unneeded_info >> down >> unneeded_info;
        symmetriseErrors(up, down, &sigma, &datshift);
        fSys[j][i].mult = sigma;
        fSys[j][i].add = fSys[j][i].mult*fData[j]/100;
        fSys[j][i].type = MULT;
        fSys[j][i].name = "CORR";    
        shift[j] += datshift; 
      }
      
      lstream >> fSys[3][i].mult >> unneeded_info;
      fSys[3][i].add = fSys[3][i].mult*fData[3]/100;
      fSys[3][i].type = MULT;
      fSys[3][i].name = "CORR";
    }
    else if (i==(fNData+5))
    {
      for (int j=0; j<2; j++)
      {
        lstream >> up >> unneeded_info >> down >> unneeded_info;
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

      lstream >> up >> unneeded_info >> down;
      symmetriseErrors(up, down, &sigma, &datshift);
      fSys[3][i].mult = sigma;
      fSys[3][i].add = fSys[3][i].mult*fData[3]/100;
      fSys[3][i].type = MULT;
      fSys[3][i].name = "CORR";
      shift[3] += datshift;
    }
    else
    {
      lstream >> fSys[0][i].mult >> unneeded_info >> fSys[1][i].mult >> unneeded_info >> fSys[2][i].mult >> unneeded_info >> fSys[3][i].mult >> unneeded_info;
      for (int j=0; j<4; j++)
      {
        fSys[j][i].add = fSys[j][i].mult*fData[j]/100;
        fSys[j][i].type = MULT;
        fSys[j][i].name = "CORR";
      }
    }
//    for (int j=0; j<4; j++)
//       {
//         lstream >> fSys[j][i].mult >> unneeded_info;
//         fSys[j][i].add = fSys[j][i].mult*fData[j]/100;
//         fSys[j][i].type = MULT;
//         fSys[j][i].name = "CORR";
//       }
  }

  for (int i=0; i<4; i++)
  {
    fData[i] *= (1.0 + shift[i]*0.01); // Shift of central value due to asymmetric errors
  }
}
