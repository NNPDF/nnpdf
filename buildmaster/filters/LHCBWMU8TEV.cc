/* 
   This file implements the W production subset of the LHCBWZMU8TEV data set.
   This is required to separate CC DY from NC DY.
   Implemented by ERN June 2023.
*/

#include "LHCBWMU8TEV.h"

void LHCBWMU8TEVFilter::ReadData()
{
  // Opening files
  fstream fZ, fWp, fWm, fCorr;

  stringstream DataFileZ("");
  DataFileZ << dataPath() << "rawdata/LHCBWZMU8TEV/LHCBWZMU8TEV_zrap.data";
  fZ.open(DataFileZ.str().c_str(), ios::in);

  if (fZ.fail()) {
    cerr << "Error opening data file " << DataFileZ.str() << endl;
    exit(-1);
  }

  stringstream DataFileWp("");
  DataFileWp << dataPath() << "rawdata/LHCBWZMU8TEV/LHCBWZMU8TEV_wplrap.data";
  fWp.open(DataFileWp.str().c_str(), ios::in);

  if (fWp.fail()) {
    cerr << "Error opening data file " << DataFileWp.str() << endl;
    exit(-1);
  }

  stringstream DataFileWm("");
  DataFileWm << dataPath() << "rawdata/LHCBWZMU8TEV/LHCBWZMU8TEV_wmlrap.data";
  fWm.open(DataFileWm.str().c_str(), ios::in);

  if (fWm.fail()) {
    cerr << "Error opening data file " << DataFileWm.str() << endl;
  }

  stringstream DataFileCorr("");
  DataFileCorr << dataPath() << "rawdata/LHCBWZMU8TEV/LHCBWZMU8TEV_corrmat.data";
  fCorr.open(DataFileCorr.str().c_str(), ios::in);

  if (fCorr.fail()) {
    cerr << "Error opening data file " << DataFileCorr.str() << endl;
    exit(-1);
  }

  // Starting filter
  const int ndata_z  = 18;
  const int ndata_wp =  8;
  const int ndata_wm =  8;
  int ndata_tot = ndata_z + ndata_wp + ndata_wm;
  const double pb2fb = 1000.; // Must multiply from pb to fb
  double MW2 = pow(MW,2.0);
  double MZ2 = pow(MZ,2.0);
  double s = 8000;
  string line;

  double totsys_w[ndata_wp+ndata_wp];
  double totsys_z[ndata_z];
  double totsys[ndata_tot];
  double inmat[ndata_tot][ndata_tot];
  double etaavg, etamin, etamax;
  int idat = 0;

  double ddum;
  
  // Z
  getline(fZ,line);
  for (int i = 0; i < ndata_z; i++)
  {
    getline(fZ,line);
    istringstream lstream(line);

    lstream >> ddum >> ddum >> ddum >> ddum >> ddum >> totsys_z[i];
    //Chck if central value is needed or not
    totsys_z[i] *= pb2fb;
  }

  // W+
  getline(fWp,line);
  for (int i = 0; i < ndata_wp; i++)
  {
    getline(fWp,line);
    istringstream lstream(line);

    lstream >> etaavg >> etamin >> etamax;
    fKin1[idat+i] = etaavg;     // eta
    fKin2[idat+i] = MW2;        // W mass squared
    fKin3[idat+i] = s;          // sqrt(s)

    lstream >> fData[idat+i];
    fData[idat+i] *= pb2fb;
    lstream >> fStat[idat+i];
    fStat[idat+i] *= pb2fb;
    lstream >> totsys_w[idat+i];
    totsys_w[idat+i] *= pb2fb;

    lstream >> fSys[idat+i][ndata_tot].add;  // Beam uncertainty
    fSys[idat+i][ndata_tot].add *= pb2fb;
    fSys[idat+i][ndata_tot].mult = fSys[idat+i][ndata_tot].add/fData[idat+i]*1e2;
    fSys[idat+i][ndata_tot].type = MULT;
    fSys[idat+i][ndata_tot].name = "LHCBBEAM8TEV";

    lstream >> fSys[idat+i][ndata_tot+1].add;  // Lumi uncertainty
    fSys[idat+i][ndata_tot+1].add *= pb2fb;
    fSys[idat+i][ndata_tot+1].mult = fSys[idat+i][ndata_tot+1].add/fData[idat+i]*1e2;
    fSys[idat+i][ndata_tot+1].type = MULT;
    fSys[idat+i][ndata_tot+1].name = "LHCBLUMI8TEV";
  }
  idat+=ndata_wp;

  // W-
  getline(fWm,line);
  for (int i = 0; i < ndata_wp; i++)
  {
    getline(fWm,line);
    istringstream lstream(line);

    lstream >> etaavg >> etamin >> etamax;
    fKin1[idat+i] = etaavg;     // eta
    fKin2[idat+i] = MW2;        // W mass squared
    fKin3[idat+i] = s;          // sqrt(s)

    lstream >> fData[idat+i];
    fData[idat+i] *= pb2fb;
    lstream >> fStat[idat+i];
    fStat[idat+i] *= pb2fb;
    lstream >> totsys_w[idat+i];
    totsys_w[idat+i] *= pb2fb;

    lstream >> fSys[idat+i][ndata_tot].add;  // Beam uncertainty
    fSys[idat+i][ndata_tot].add *= pb2fb;
    fSys[idat+i][ndata_tot].mult = fSys[idat+i][ndata_tot].add/fData[idat+i]*1e2;
    fSys[idat+i][ndata_tot].type = MULT;
    fSys[idat+i][ndata_tot].name = "LHCBBEAM8TEV";

    lstream >> fSys[idat+i][ndata_tot+1].add;  // Lumi uncertainty
    fSys[idat+i][ndata_tot+1].add *= pb2fb;
    fSys[idat+i][ndata_tot+1].mult = fSys[idat+i][ndata_tot+1].add/fData[idat+i]*1e2;
    fSys[idat+i][ndata_tot+1].type = MULT;
    fSys[idat+i][ndata_tot+1].name = "LHCBLUMI8TEV";
  }
  idat+=ndata_wm;

  for(int i=0; i<ndata_z; i++)
    totsys[i] = totsys_z[i];
  for(int i=ndata_z; i<ndata_tot; i++)
    totsys[i] = totsys_w[i-ndata_z];

  // Reading Covariance Matrix
  for (int i = 0; i < ndata_tot; i++) {
    for (int j = 0; j < i+1; j++) {             // read only lower triangle
      fCorr >> inmat[j][i];
    }
  }
  for (int i = 0; i < ndata_tot; i++) {
    for (int j = i+1; j < ndata_tot; j++) {
      inmat[j][i] = inmat[i][j];               // symmetrize
    }
  }

  //  Multiply by total systematic uncertainty
  double** covmat = new double*[ndata_tot];
  for(int i = 0; i < ndata_tot; i++)
  {
    covmat[i] = new double[ndata_tot];
    for(int j = 0; j < ndata_tot; j++)
      {
	covmat[i][j]=inmat[i][j]*totsys[i]*totsys[j];
      }
  }

  // Generate artificial systematics
  double** syscor = new double*[ndata_tot];
  for(int i = 0; i < ndata_tot; i++)
    syscor[i] = new double[ndata_tot];

  if(!genArtSys(ndata_tot,covmat,syscor))
   {
     cerr << " in " << fSetName << endl;
     exit(-1);
   }

  for (int i = 0; i < fNData; i++)
    for (int l = 0; l < fNSys-2; l++)
    {
      fSys[i][l].add  = syscor[i+ndata_z][l];
      fSys[i][l].mult = fSys[i][l].add/fData[i]*1e2;
      fSys[i][l].type = MULT;
      ostringstream sysname;
      sysname << "LHCBWZMU8TEV_" << l;
      fSys[i][l].name = sysname.str();
    }

  for(int i = 0; i < ndata_tot; i++) 
    {
      delete[] syscor[i];
      delete[] covmat[i];
    }
  delete[] syscor;
  delete[] covmat;
  
  fZ.close();
  fWp.close();
  fWm.close();
  fCorr.close();
}
