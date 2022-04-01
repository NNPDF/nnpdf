/****************
 * LHCb 2fb^{-1}
 *
 * Z production > e+e- data from the LHCb experiment
 *                8 TeV data
 * Final data from the LHCb preprint: 1503.00963
 * Luminosity uncertainty is a 1.2% in all data points and is quoted separately
 * from other sources of systematic uncertainty
 */

#include "LHCb.h"

void LHCBZEE2FBFilter::ReadData()
{
  // Opening files
  fstream f1, f2;

  stringstream datafile1("");
  datafile1 << dataPath() << "rawdata/"
  << fSetName << "/lhcb_zrap.data";
  f1.open(datafile1.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile1.str() << endl;
    exit(-1);
  }

  stringstream datafile2("");
  datafile2 << dataPath() << "rawdata/"
  << fSetName << "/lhcb_zrap.corr";
  f2.open(datafile2.str().c_str(), ios::in);

  if (f2.fail()) {
    cerr << "Error opening data file " << datafile2.str() << endl;
    exit(-1);
  }


  // Initialize
  const double convfac = 1000.; // Must multiply from pb to fb
  double MZ2 = pow(MZ,2.0);
  double s = 8000.;

  double etamin,etamax;
  double sys_corr[fNData],lumi[fNData];//,tot_unc[fNData];
  double inmat[fNData][fNData];
  string line;

  for (int i = 0; i < fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);

      lstream >> etamin >> etamax;
      fKin1[i] = (etamin+etamax) * 0.5;       // <eta>
      fKin2[i] = MZ2;
      fKin3[i] = s;

      lstream >> fData[i];
      fData[i] *= convfac;
      // stat uncertainty
      lstream >> fStat[i];
      fStat[i] *= convfac;

      // Uncorrelated systematics
      lstream >> fSys[i][0].add;
      fSys[i][0].add *= convfac;
      fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
      fSys[i][0].type = MULT;
      fSys[i][0].name = "UNCORR";

      // Correlated and luminosity uncertainties
      lstream >> sys_corr[i] >> lumi[i];
      sys_corr[i] *= convfac;
      lumi[i] *= convfac;

      // Total uncertainty
      //tot_unc[i] = pow(sys_corr[i]*sys_corr[i]+fSys[i][0].add*fSys[i][0].add+fStat[i]*fStat[i],0.5);

      // Normalization uncertainty is defined both MULT (in %) and ADD (in abs value)
      fSys[i][fNSys-1].add   = lumi[i];
      fSys[i][fNSys-1].mult  = fSys[i][fNSys-1].add/fData[i]*1e2;
      fSys[i][fNSys-1].type = MULT;
      fSys[i][fNSys-1].name = "LHCBLUMI8TEV";
    }


  // Reading Covariance Matrix
  for (int i = 0; i < fNData; i++) {
    for (int j = 0; j < i+1; j++) {             // read only lower triangle
      f2 >> inmat[j][i];
    }
  }
  for (int i = 0; i < fNData; i++) {
    for (int j = i+1; j < fNData; j++) {
      inmat[j][i] = inmat[i][j];               // symmetrize
    }
  }

  //  Multiply by total systematic uncertainty
  double** covmat = new double*[fNData];
  for(int i = 0; i < fNData; i++)
  {
    covmat[i] = new double[fNData];
    for(int j = 0; j < fNData; j++) {
      covmat[i][j]=inmat[i][j]*sys_corr[i]*sys_corr[j];
      }
  }

  // Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
   {
     cerr << " in " << fSetName << endl;
     exit(-1);
   }

  for (int i = 0; i < fNData; i++)
    for (int l = 1; l < fNSys-1; l++)
    {
      fSys[i][l].add  = syscor[i][l-1];
      fSys[i][l].mult = fSys[i][l-1].add/fData[i]*1e2;
      fSys[i][l].type = MULT;
      fSys[i][l].name = "CORR";
    }

  for(int i = 0; i < fNData; i++) 
    {
      delete[] syscor[i];
      delete[] covmat[i];
    }
  delete[] syscor;
  delete[] covmat;
  
  f1.close();
  f2.close();
}

// LHCBZEE2FB with fix described in https://github.com/NNPDF/nnpdf/issues/904
void LHCBZEE2FB_40Filter::ReadData()
{
  // Opening files
  fstream f1, f2;

  stringstream datafile1("");
  // use same rawdata as original filter.
  datafile1 << dataPath() << "rawdata/LHCBZEE2FB/lhcb_zrap.data";
  f1.open(datafile1.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile1.str() << endl;
    exit(-1);
  }

  stringstream datafile2("");
  datafile2 << dataPath() << "rawdata/LHCBZEE2FB/lhcb_zrap.corr";
  f2.open(datafile2.str().c_str(), ios::in);

  if (f2.fail()) {
    cerr << "Error opening data file " << datafile2.str() << endl;
    exit(-1);
  }


  // Initialize
  const double convfac = 1000.; // Must multiply from pb to fb
  double MZ2 = pow(MZ,2.0);
  double s = 8000.;

  double etamin,etamax;
  double sys_corr[fNData],lumi[fNData];//,tot_unc[fNData];
  double inmat[fNData][fNData];
  string line;

  for (int i = 0; i < fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);

      lstream >> etamin >> etamax;
      fKin1[i] = (etamin+etamax) * 0.5;       // <eta>
      fKin2[i] = MZ2;
      fKin3[i] = s;

      lstream >> fData[i];
      fData[i] *= convfac;
      // stat uncertainty
      lstream >> fStat[i];
      fStat[i] *= convfac;

      // Uncorrelated systematics
      lstream >> fSys[i][0].add;
      fSys[i][0].add *= convfac;
      fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
      fSys[i][0].type = MULT;
      fSys[i][0].name = "UNCORR";

      // Correlated and luminosity uncertainties
      lstream >> sys_corr[i] >> lumi[i];
      sys_corr[i] *= convfac;
      lumi[i] *= convfac;

      // Total uncertainty
      //tot_unc[i] = pow(sys_corr[i]*sys_corr[i]+fSys[i][0].add*fSys[i][0].add+fStat[i]*fStat[i],0.5);

      // Normalization uncertainty is defined both MULT (in %) and ADD (in abs value)
      fSys[i][fNSys-1].add   = lumi[i];
      fSys[i][fNSys-1].mult  = fSys[i][fNSys-1].add/fData[i]*1e2;
      fSys[i][fNSys-1].type = MULT;
      fSys[i][fNSys-1].name = "LHCBLUMI8TEV";
    }


  // Reading Covariance Matrix
  for (int i = 0; i < fNData; i++) {
    for (int j = 0; j < i+1; j++) {             // read only lower triangle
      f2 >> inmat[j][i];
    }
  }
  for (int i = 0; i < fNData; i++) {
    for (int j = i+1; j < fNData; j++) {
      inmat[j][i] = inmat[i][j];               // symmetrize
    }
  }

  //  Multiply by total systematic uncertainty
  double** covmat = new double*[fNData];
  for(int i = 0; i < fNData; i++)
  {
    covmat[i] = new double[fNData];
    for(int j = 0; j < fNData; j++) {
      covmat[i][j]=inmat[i][j]*sys_corr[i]*sys_corr[j];
      }
  }

  // Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
   {
     cerr << " in " << fSetName << endl;
     exit(-1);
   }

  for (int i = 0; i < fNData; i++)
    for (int l = 1; l < fNSys-1; l++)
    {
      fSys[i][l].add  = syscor[i][l-1];
      fSys[i][l].mult = fSys[i][l].add/fData[i]*1e2;
      fSys[i][l].type = MULT;
      fSys[i][l].name = "CORR";
    }

  for(int i = 0; i < fNData; i++) 
    {
      delete[] syscor[i];
      delete[] covmat[i];
    }
  delete[] syscor;
  delete[] covmat;
  
  f1.close();
  f2.close();
}
