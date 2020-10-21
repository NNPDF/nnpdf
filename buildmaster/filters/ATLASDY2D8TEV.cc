#include "ATLASDY2D8TEV.h"

/**
 * Data from Table 3 of arXiv:1606.01736 from hepdata /HepDat/9203/d2-x1-y1
 * Using xfitter format file
 */
void ATLASDY2D8TEVFilter::ReadData()
{
  // Opening files
  ifstream f1, f2;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"<< fSetName << "/HMDYMeasurement.rapidity";
  f1.open(datafile.str().c_str());

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream datafile2("");
  datafile2 << dataPath() << "rawdata/"<< fSetName << "/hepdata.txt";
  f2.open(datafile2.str().c_str());

  if (f2.fail()) {
    cerr << "Error opening data file " << datafile2.str() << endl;
    exit(-1);
  }

  string line, kinematics, tmp, name;
  double mllbin_low, mllbin_high, ybin_low, ybin_high;
  const double s = 8000;

  // skip header
  for (int i = 0; i < 15; i++)
    getline(f1, line);

  for (int i = 0; i < 13; i++)
    getline(f2, kinematics);

  // Filter data file
  for (int idat = 0; idat < fNData; idat++)
  {
    getline(f1, line);
    getline(f2, kinematics);
    istringstream lstream(line);
    istringstream lkstream(kinematics);

    lstream >> fKin2[idat] >> fKin1[idat] >> fData[idat] >> fStat[idat];
    lkstream >> tmp >> mllbin_low >> mllbin_high >> tmp >> ybin_low >> ybin_high;

    // ugly normalization to match applgrid predictions
    fData[idat] *= 2 * (mllbin_high-mllbin_low);
    
    fStat[idat] = fStat[idat] * fData[idat] * 1e-2;

    for (int l = 0; l < fNSys; l++)
      {
        lstream >> fSys[idat][l].mult;
        fSys[idat][l].add = fSys[idat][l].mult * fData[idat] * 1e-2;
        fSys[idat][l].type = MULT;

        if (l == 0)            name = "UNCORR";
        else if (l == fNSys-1) name = "ATLASLUMI12";          
        else                   name = "CORR";
	
	fSys[idat][l].name = name;
      }

    fKin2[idat] *= fKin2[idat];
    fKin3[idat] = s;
  }

  f1.close();
  f2.close();
}

