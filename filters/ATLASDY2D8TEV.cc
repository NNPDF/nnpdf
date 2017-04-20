#include "ATLASDY2D8TEV.h"

/**
 * Data from Table 3 of arXiv:1606.01736 from hepdata /HepDat/9203/d2-x1-y1
 */
void ATLASDY2D8TEVFilter::ReadData()
{
  // Opening files
  ifstream f1;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"<< fSetName << "/hepdata.txt";
  f1.open(datafile.str().c_str());

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  string line, tmp;
  const double s = 8000;

  // skip header
  for (int i = 0; i < 13; i++)
    getline(f1, line);

  // Filter data file
  for (int idat = 0; idat < fNData; idat++)
  {
    getline(f1, line);
    istringstream lstream(line);

    lstream >> fKin2[idat] >> tmp >> tmp
            >> fKin1[idat] >> tmp >> tmp
            >> fData[idat] >> fStat[idat] >> tmp;

    for (int l = 0; l < fNSys; l++)
      {
        lstream >> fSys[idat][l].add >> tmp;
        fSys[idat][l].mult = fSys[idat][l].add*100/fData[idat];
        fSys[idat][l].type = MULT;

        if (l == 0)
          fSys[idat][l].name = "UNCORR";
        else if (l == fNSys-1)
          {
            fSys[idat][l].add = fabs(fSys[idat][l].add);
            fSys[idat][l].mult = fabs(fSys[idat][l].mult);
            fSys[idat][l].name = "ATLASLUMI12";
          }
        else
          fSys[idat][l].name = "CORR";
      }


    fKin2[idat] *= fKin2[idat];
    fKin3[idat] = s;
  }

  f1.close();
}

