/*

  Include some general comments here
  - paper references etc

  */

#include "ATLAS2011JETS.h"

void ATLAS1JET11Filter::Loop(fstream & file, fstream & fileew, double rap, int nDataMin, int nDataMax, int rap_bin)
{

  // Reading files
  string line;
  int dummylines = 31;
  double ptmin, ptmax; 
  double dummy;
  double extra; 
  double NPcorr;
  double sigma;
  double stat1, stat2;

  //EW corrections
  int dummyEWlines = 16;
  double * EWcorr = new double [nDataMax-nDataMin];

  double * NPerr = new double [nDataMax-nDataMin];
  
  for (int i = 0; i < dummyEWlines; i++)
  {
    string dummyLine;
    getline(fileew, dummyLine);
  }

  for (int i = nDataMin; i < nDataMax; i++)  // Loop over datapoints
  {
    double delta;
    double NPerr1;
    double NPerr2;
    fileew >> dummy;
    fileew >> dummy;
    fileew >> dummy;
    fileew >> NPerr1;
    fileew >> NPerr2;
    fileew >> EWcorr[i-nDataMin];
    fileew >> dummy;
    fileew >> dummy;
    symmetriseErrors(NPerr1,NPerr2,&NPerr[i-nDataMin],&delta);
    
  }

  //Raw data
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
    fData[i] *= EWcorr[i-nDataMin];// apply corrections (i.e. EW)

    fStat[i] = sqrt(pow(stat1,2)+pow(stat2,2));//statistical uncertainty (symmetric)
    fStat[i] = fStat[i]*fData[i]*1e-2;

    int Nsysprime = 67; //number of systematics in each rapidity bin (excluding lumi)
    //fSys is equal to 0 for the other rapidities
    
    for (int l = 0 ; l < fNSys-1; ++l)
    {
      fSys[i][l].mult = 0;
      fSys[i][l].add = 0;
      fSys[i][l].type = MULT;
      fSys[i][l].name = "CORR";
    }
    
    for (int l = rap_bin*Nsysprime; l < (rap_bin+1)*Nsysprime; ++l)
    {
      file >> fSys[i][l].mult;
      fSys[i][l].add = fSys[i][l].mult*fData[i]*1e-2;
    }

    //now fill the lumi 
    file >> fSys[i][fNSys-2].mult;
    fSys[i][fNSys-2].type = MULT;
    fSys[i][fNSys-2].add = fSys[i][fNSys-2].mult*fData[i]*1e-2;
    fSys[i][fNSys-2].name = "ATLASLUMI11";

    fSys[i][fNSys-1].mult = NPerr[i-nDataMin];
    fSys[i][fNSys-1].type = MULT;
    fSys[i][fNSys-1].name = "ATLAS1JET11_NP_err";
    fSys[i][fNSys-1].add = fSys[i][fNSys-1].mult*fData[i]*1e-2;

  }

  delete[] EWcorr;
  delete[] NPerr;
  
}

void ATLAS1JET11Filter::ReadData()
{
 
  // Open and read rawdata here
  fstream f1, f2, f3, f4, f5, f6;
  fstream few1, few2, few3, few4, few5, few6;

  stringstream rapbin1("");
  rapbin1 << dataPath() << "rawdata/"
  << fSetName << "/incjets_7000_R06_rapidity_y_00_05.dat";
  f1.open(rapbin1.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << rapbin1.str() << endl;
    exit(-1);
  }

  stringstream ewrapbin1("");
  ewrapbin1 << dataPath() << "rawdata/"
  << fSetName << "/NPandEW_R06_Eta1.txt";
  few1.open(ewrapbin1.str().c_str(), ios::in);

  if (few1.fail()) {
    cerr << "Error opening data file " << ewrapbin1.str() << endl;
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

  stringstream ewrapbin2("");
  ewrapbin2 << dataPath() << "rawdata/"
  << fSetName << "/NPandEW_R06_Eta2.txt";
  few2.open(ewrapbin2.str().c_str(), ios::in);

  if (few2.fail()) {
    cerr << "Error opening data file " << ewrapbin2.str() << endl;
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

  stringstream ewrapbin3("");
  ewrapbin3 << dataPath() << "rawdata/"
  << fSetName << "/NPandEW_R06_Eta3.txt";
  few3.open(ewrapbin3.str().c_str(), ios::in);

  if (few3.fail()) {
    cerr << "Error opening data file " << ewrapbin3.str() << endl;
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

  stringstream ewrapbin4("");
  ewrapbin4 << dataPath() << "rawdata/"
  << fSetName << "/NPandEW_R06_Eta4.txt";
  few4.open(ewrapbin4.str().c_str(), ios::in);

  if (few4.fail()) {
    cerr << "Error opening data file " << ewrapbin4.str() << endl;
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

  stringstream ewrapbin5("");
  ewrapbin5 << dataPath() << "rawdata/"
  << fSetName << "/NPandEW_R06_Eta5.txt";
  few5.open(ewrapbin5.str().c_str(), ios::in);

  if (few5.fail()) {
    cerr << "Error opening data file " << ewrapbin5.str() << endl;
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

  stringstream ewrapbin6("");
  ewrapbin6 << dataPath() << "rawdata/"
  << fSetName << "/NPandEW_R06_Eta6.txt";
  few6.open(ewrapbin6.str().c_str(), ios::in);

  if (few6.fail()) {
    cerr << "Error opening data file " << ewrapbin6.str() << endl;
    exit(-1);
  }

  //now read the data
  Loop(f1, few1, 0.25, 0, 31, 0);
  Loop(f2, few2, 0.75, 31, 29+31, 1);
  Loop(f3, few3, 1.25, 31+29, 29+31+26, 2);
  Loop(f4, few4, 1.75, 31+29+26, 29+31+26+23, 3);
  Loop(f5, few5, 2.25, 31+29+26+23, 29+31+26+23+19, 4);
  Loop(f6, few6, 2.75, 31+29+26+23+19, 29+31+26+23+19+12, 5);
  
}
