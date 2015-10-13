/*

  Include some general comments here
  - paper references etc

  */

#include "ATLAS2011JETS.h"

void ATLASR06JETS2011Filter::Loop(fstream & file, double rap, int nDataMin, int nDataMax)
{

  // Reading files
  string line;
  double dummylines = 31;
  double ptmin, ptmax; 
  double dummy;
  double extra; 
  double NPcorr;
  double sigma;
  double stat1, stat2;

  //double lcorr = 1.0187; // correction factor due to luminosity upgrade
  for (int i = 0; i < dummylines; i++)
  {
    string dummyLine;
    getline(file, dummyLine);
  }

  for (int i = nDataMin; i < nDataMax; i++)  // Loop over datapoints
  {
   // Jets usually y, pT^2, s - midpoints of bins
    fKin1[i] = rap;
    file >> dummy;
    file >> ptmin;
    file >> ptmax;
    file >> extra;
    file >> NPcorr;
    file >> sigma;
    file >> stat1;
    file >> stat2;

    fKin2[i] = pow((ptmax+ptmin)/2.,2);
    fKin3[i] = 7000; 

    fData[i] = sigma;// fill datapoints
    fData[i] *= NPcorr;// apply corrections (i.e. non-perturbative)

    fStat[i] = sqrt(pow(stat1,2)+pow(stat2,2));//statistical uncertainty (symmetric)
    fStat[i] = fStat[i]*fData[i]*1e-2;

    for (int l = 0; l < fNSys; ++l)
    {
      file >> fSys[i][l].mult;
      fSys[i][l].type = MULT;
      fSys[i][l].name = "UNCORR";
      fSys[i][l].add = fSys[i][l].mult*fData[i]*1e-2;
    }

  }

}

void ATLASR06JETS2011Filter::ReadData()
{
 
  // Open and read rawdata here
  fstream f1, f2, f3, f4, f5, f6;

  stringstream rapbin1("");
  rapbin1 << dataPath() << "rawdata/"
  << fSetName << "/incjets_7000_R06_rapidity_y_00_05.dat";
  f1.open(rapbin1.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << rapbin1.str() << endl;
    exit(-1);
  }

  stringstream rapbin2("");
  rapbin2 << dataPath() << "rawdata/"
  << fSetName << "/incjets_7000_R06_rapidity_y_05_10.dat";
  f2.open(rapbin2.str().c_str(), ios::in);

  if (f2.fail()) {
    cerr << "Error opening data file " << rapbin2.str() << endl;
    exit(-1);
  }

  stringstream rapbin3("");
  rapbin3 << dataPath() << "rawdata/"
  << fSetName << "/incjets_7000_R06_rapidity_y_10_15.dat";
  f3.open(rapbin3.str().c_str(), ios::in);

  if (f3.fail()) {
    cerr << "Error opening data file " << rapbin3.str() << endl;
    exit(-1);
  }
  
  stringstream rapbin4("");
  rapbin4 << dataPath() << "rawdata/"
  << fSetName << "/incjets_7000_R06_rapidity_y_15_20.dat";
  f4.open(rapbin4.str().c_str(), ios::in);

  if (f4.fail()) {
    cerr << "Error opening data file " << rapbin4.str() << endl;
    exit(-1);
  }

  stringstream rapbin5("");
  rapbin5 << dataPath() << "rawdata/"
  << fSetName << "/incjets_7000_R06_rapidity_y_20_25.dat";
  f5.open(rapbin5.str().c_str(), ios::in);

  if (f5.fail()) {
    cerr << "Error opening data file " << rapbin5.str() << endl;
    exit(-1);
  }

  stringstream rapbin6("");
  rapbin6 << dataPath() << "rawdata/"
  << fSetName << "/incjets_7000_R06_rapidity_y_25_30.dat";
  f6.open(rapbin6.str().c_str(), ios::in);

  if (f6.fail()) {
    cerr << "Error opening data file " << rapbin6.str() << endl;
    exit(-1);
  }

  //now read the data
  Loop(f1, 0.25, 0, 31);
  Loop(f2, 0.75, 31, 29+31);
  Loop(f3, 1.25, 31+29, 29+31+26);
  Loop(f4, 1.75, 31+29+26, 29+31+26+23);
  Loop(f5, 2.25, 31+29+26+23, 29+31+26+23+19);
  Loop(f6, 2.75, 31+29+26+23+19, 29+31+26+23+19+12);
  
}
