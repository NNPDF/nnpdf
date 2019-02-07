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
  // Create stream to read data file
  fstream f1;

  // Data file
  stringstream datafile("");
  string filename;
  filename = "ATLAS_SINGLETOP_TCH_DIFF_7TEV_T_RAP";
  datafile << dataPath()
           << "rawdata/" << filename << "/" << filename << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
  {
    cerr << "Error opening data file " << datafile.str() << endl;
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
    lstream >> unneeded_info;
    lstream >> fstat_percentage;
    fStat[i] = fstat_percentage*fData[i]/100; // Absolute statistical uncertainty
    lstream >> unneeded_info >> unneeded_info;

    lstream >> fSys[i][0].mult; // Percentage total systematic uncertainty
    fSys[i][0].add = fSys[i][0].mult*fData[i]/100; // Absolute total systematic uncertainty
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CORR";
  }
}
