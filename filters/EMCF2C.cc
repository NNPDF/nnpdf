/**
 * EMCF2C.cc
 *
 *  Production of Charmed Particles in 250-GeV mu+-iron interactions
 *  The European Muon Collaboration (J.J. Auber et al.)
 *  Nucl. Phys. B213 (1983) 31
 *
 */

#include "EMCF2C.h"

void EMCF2CFilter::ReadData()
{
  // Opening files
  fstream f1;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/EMCF2C.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double syscor[fNData][fNSys];

  for (int i = 0; i < fNData; i++)
  {
    for (int l = 0; l < fNSys; l++)
      syscor[i][l] = 0.0;
  }


/*
  // Reading data
  index = 0;
  string tmp;
  getline(f1,tmp);
  {
    for (int i = 0; i < fNData; i++)
    for (int j = 0; j < 52; j++)
      f1 >> datain[i][j];
    if (datain[i][0] != i+1) {
      cerr << "Mismatch" << endl;
      exit(-1);
    }
    index++;
  }

  if (index != fNData)
  {
    cerr << "Mismatch in the number of data points" << endl;
    exit(-1);
  }
*/
  // Filtering data
  double x[fNData], y[fNData], q2[fNData];
  double xmin, xmax, q2min, q2max;


  string tmp;
  getline(f1,tmp);
  getline(f1,tmp);
  for (int i = 0; i < fNData; i++)
  {
    f1 >> xmin >> xmax >> q2min >> q2max >> fData[i] >> fStat[i];
    x[i]   = (xmin + xmax)/2.;
    q2[i]  = (q2min + q2max)/2.;
    y[i]   = 0.0;

    // Uncorrelated systematics
//    fSys[i][1].mult = datain[i][6];
//    fSys[i][1].type = ADD;
//    fSys[i][1].name = "UNCORR";

    // Systematic errors
//    double sys_tot = 0.0;
//    for (int l = 0; l < fNSys-2; l++)
//    {
//      syscor[i][l] = datain[i][l+8];
//      sys_tot += syscor[i][l]*syscor[i][l];
//    }
//    sys_tot = sqrt(sys_tot);

    // Normalization uncertainty
//    fSys[i][0].mult = absnorexp;
//    fSys[i][0].type = MULT;
//    fSys[i][0].name = "CORR";

    // Check total error
//    double tot = datain[i][7];
//    double tot_tmp = sqrt(sys_tot*sys_tot + stat[i]*stat[i] + fSys[i][1].mult*fSys[i][1].mult);
//    if (fabs(tot-tot_tmp) > 1.4)
//    {
//      cerr << "Mismatch " << i << "\t" << tot << "\t" << tot_tmp << endl;
//      exit(-1);
//    }

    fKin1[i] = x[i];
    fKin2[i] = q2[i];
    fKin3[i] = y[i];

    //rescaling for BR - check
    fData[i] = fData[i]/0.8;
    sist = fData[i]*0.15;
    fSys[i][0].add = sist;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CORR_EMC";
    fSys[i][0].mult = fSys[i][0].add/(fData[i]*1e-2);
  }

  f1.close();
}
