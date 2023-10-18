/* 
   This file implements the W production subset of the LHCBWZMU7TEV data set.
   This is required to separate NC DY from CC DY.
   Implemented by ERN June 2023.
*/

#include "LHCBZMU7TEV.h"

void LHCBZMU7TEVFilter::ReadData()
{
  // Opening files
  fstream fZ, fWp, fWm, fCorr;

  stringstream DataFileZ("");
  DataFileZ << dataPath() << "rawdata/LHCBWZMU7TEV/LHCBWZMU7TEV_zrap.data";
  fZ.open(DataFileZ.str().c_str(), ios::in);

  if (fZ.fail()) {
    cerr << "Error opening data file " << DataFileZ.str() << endl;
    exit(-1);
  }

  stringstream DataFileWp("");
  DataFileWp << dataPath() << "rawdata/LHCBWZMU7TEV/LHCBWZMU7TEV_wplrap.data";
  fWp.open(DataFileWp.str().c_str(), ios::in);

  if (fWp.fail()) {
    cerr << "Error opening data file " << DataFileWp.str() << endl;
    exit(-1);
  }

  stringstream DataFileWm("");
  DataFileWm << dataPath() << "rawdata/LHCBWZMU7TEV/LHCBWZMU7TEV_wmlrap.data";
  fWm.open(DataFileWm.str().c_str(), ios::in);

  if (fWm.fail()) {
    cerr << "Error opening data file " << DataFileWm.str() << endl;
  }

  stringstream DataFileCorr("");
  DataFileCorr << dataPath() << "rawdata/LHCBWZMU7TEV/LHCBWZMU7TEV_corrmat.data";
  fCorr.open(DataFileCorr.str().c_str(), ios::in);

  if (fCorr.fail()) {
    cerr << "Error opening data file " << DataFileCorr.str() << endl;
    exit(-1);
  }

  // Starting filter
  const int ndata_z  = 17;
  const int ndata_wp =  8;
  const int ndata_wm =  8;
  int ndata_tot = ndata_z + ndata_wp + ndata_wm;
  const double pb2fb = 1000.; // Must multiply from pb to fb
  double MW2 = pow(MW,2.0);
  double MZ2 = pow(MZ,2.0);
  double s = 7000;
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

    lstream >> etaavg >> etamin >> etamax;
    fKin1[i] = etaavg;     // eta
    fKin2[i] = MZ2;        // Z mass squared
    fKin3[i] = s;          // sqrt(s)

    lstream >> fData[i];
    fData[i] *= pb2fb;
    lstream >> fStat[i];
    fStat[i] *= pb2fb;
    lstream >> totsys_z[i];
    totsys_z[i] *= pb2fb;

    lstream >> fSys[i][ndata_tot].add;  // Beam uncertainty
    fSys[i][ndata_tot].add *= pb2fb;
    fSys[i][ndata_tot].mult = fSys[idat+i][ndata_tot].add/fData[i]*1e2;
    fSys[i][ndata_tot].type = MULT;
    fSys[i][ndata_tot].name = "LHCBBEAM7TEV";

    lstream >> fSys[i][ndata_tot+1].add;  // Lumi uncertainty
    fSys[i][ndata_tot+1].add *= pb2fb;
    fSys[i][ndata_tot+1].mult = fSys[i][ndata_tot+1].add/fData[i]*1e2;
    fSys[i][ndata_tot+1].type = MULT;
    fSys[i][ndata_tot+1].name = "LHCBLUMI7TEV";
  }
  
  // W+
  getline(fWp,line);
  for (int i = 0; i < ndata_wp; i++)
  {
    getline(fWp,line);
    istringstream lstream(line);

    lstream >> ddum >> ddum >> ddum >> ddum >> ddum >> totsys_w[i+idat];
    //Chck if central value is needed or not
    totsys_w[i+idat] *= pb2fb;
  }
  idat+=ndata_wp;

  // W-
  getline(fWm,line);
  for (int i = 0; i < ndata_wp; i++)
  {
    getline(fWm,line);
    istringstream lstream(line);

    lstream >> ddum >> ddum >> ddum >> ddum >> ddum >> totsys_w[i+idat];
    //Chck if central value is needed or not
    totsys_w[i+idat] *= pb2fb;
  }
  idat+=ndata_wm;

  for(int i=0; i<ndata_z; i++)
    {
      totsys[i] = totsys_z[i];
    }
  for(int i=ndata_z; i<ndata_tot; i++)
    {
      totsys[i] = totsys_w[i-ndata_z];
    }

  for(int i=0; i< ndata_tot; i++)
    {
      cout << totsys[i] << endl;
    }
  
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
      fSys[i][l].add  = syscor[i][l];
      fSys[i][l].mult = fSys[i][l].add/fData[i]*1e2;
      fSys[i][l].type = MULT;
      ostringstream sysname;
      sysname << "LHCBWZMU7TEV_" << l;
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
