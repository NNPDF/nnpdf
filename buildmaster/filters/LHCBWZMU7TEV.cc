/********************* NEW DATA in NNPDF3.1 **************************
*
* LHCb 1fb^{-1}
*
* W production > mu nu_mu data from the LHCb experiment
*
* Final data from the LHCb preprint: 1408.4354
* Luminosity uncertainty is a 1.71% in
* all data points and is quoted separately from other sources
* of systematic uncertainty [1410.0149v2]
*/

#include "LHCb.h"

void LHCBWZMU7TEVFilter::ReadData()
{
  // Opening files
  fstream fZ, fWp, fWm, fCorr;

  stringstream DataFileZ("");
  DataFileZ << dataPath() << "rawdata/" << fSetName << "/LHCBWZMU7TEV_zrap.data";
  fZ.open(DataFileZ.str().c_str(), ios::in);

  if (fZ.fail()) {
    cerr << "Error opening data file " << DataFileZ.str() << endl;
    exit(-1);
  }

  stringstream DataFileWp("");
  DataFileWp << dataPath() << "rawdata/" << fSetName << "/LHCBWZMU7TEV_wplrap.data";
  fWp.open(DataFileWp.str().c_str(), ios::in);

  if (fWp.fail()) {
    cerr << "Error opening data file " << DataFileWp.str() << endl;
    exit(-1);
  }

  stringstream DataFileWm("");
  DataFileWm << dataPath() << "rawdata/" << fSetName << "/LHCBWZMU7TEV_wmlrap.data";
  fWm.open(DataFileWm.str().c_str(), ios::in);

  if (fWm.fail()) {
    cerr << "Error opening data file " << DataFileWm.str() << endl;
  }

  stringstream DataFileCorr("");
  DataFileCorr << dataPath() << "rawdata/" << fSetName << "/LHCBWZMU7TEV_corrmat.data";
  fCorr.open(DataFileCorr.str().c_str(), ios::in);

  if (fCorr.fail()) {
    cerr << "Error opening data file " << DataFileCorr.str() << endl;
    exit(-1);
  }

  // Starting filter
  const int ndata_z  = 17;
  const int ndata_wp =  8;
  const int ndata_wm =  8;
  //double binsize;  // Must multiply by binsize to match inclusive bin-by-bin data
  const double pb2fb = 1000.; // Must multiply from pb to fb
  double MW2 = pow(MW,2.0);
  double MZ2 = pow(MZ,2.0);
  double s = 7000;
  string line;

  double totsys[fNData];
  double inmat[fNData][fNData];
  double etaavg, etamin, etamax;
  int idat = 0;

  cout << "********** WARNING: Converting pb to fb to match ApplGrid output ********" << endl;

  // Z
  getline(fZ,line);
  for (int i = 0; i < ndata_z; i++)
  {
    getline(fZ,line);
    istringstream lstream(line);

    lstream >> etaavg >> etamin >> etamax;
    //binsize=1./(etamax-etamin);
    fKin1[idat+i] = etaavg;     // eta
    fKin2[idat+i] = MZ2;        // Z mass squared
    fKin3[idat+i] = s;          // sqrt(s)

    lstream >> fData[idat+i];
    fData[idat+i] *= pb2fb;
    //    fData[idat+i] *= binsize;
    lstream >> fStat[idat+i];
    fStat[idat+i] *= pb2fb;
    //    fStat[idat+i] *= binsize;
    lstream >> totsys[idat+i];
    totsys[idat+i] *= pb2fb;
    //    totsys[idat+i] *= binsize;

    lstream >> fSys[idat+i][fNData].add;  // Beam uncertainty
    fSys[idat+i][fNData].add *= pb2fb;
    fSys[idat+i][fNData].mult = fSys[idat+i][fNData].add/fData[idat+i]*1e2;
    fSys[idat+i][fNData].type = MULT;
    fSys[idat+i][fNData].name = "LHCBBEAM7TEV";

    lstream >> fSys[idat+i][fNData+1].add;  // Lumi uncertainty
    fSys[idat+i][fNData+1].add *= pb2fb;
    fSys[idat+i][fNData+1].mult = fSys[idat+i][fNData+1].add/fData[idat+i]*1e2;
    fSys[idat+i][fNData+1].type = MULT;
    fSys[idat+i][fNData+1].name = "LHCBLUMI7TEV";
  }
  idat+=ndata_z;

  // W+
  getline(fWp,line);
  for (int i = 0; i < ndata_wp; i++)
  {
    getline(fWp,line);
    istringstream lstream(line);

    lstream >> etaavg >> etamin >> etamax;
    //binsize=1./(etamax-etamin);
    fKin1[idat+i] = etaavg;     // eta
    fKin2[idat+i] = MW2;        // W mass squared
    fKin3[idat+i] = s;          // sqrt(s)

    lstream >> fData[idat+i];
    fData[idat+i] *= pb2fb;
    //    fData[idat+i] *= binsize;
    lstream >> fStat[idat+i];
    fStat[idat+i] *= pb2fb;
    //    fStat[idat+i] *= binsize;
    lstream >> totsys[idat+i];
    totsys[idat+i] *= pb2fb;
    //    totsys[idat+i] *= binsize;

    lstream >> fSys[idat+i][fNData].add;  // Beam uncertainty
    fSys[idat+i][fNData].add *= pb2fb;
    fSys[idat+i][fNData].mult = fSys[idat+i][fNData].add/fData[idat+i]*1e2;
    fSys[idat+i][fNData].type = MULT;
    fSys[idat+i][fNData].name = "LHCBBEAM7TEV";

    lstream >> fSys[idat+i][fNData+1].add;  // Lumi uncertainty
    fSys[idat+i][fNData+1].add *= pb2fb;
    fSys[idat+i][fNData+1].mult = fSys[idat+i][fNData+1].add/fData[idat+i]*1e2;
    fSys[idat+i][fNData+1].type = MULT;
    fSys[idat+i][fNData+1].name = "LHCBLUMI7TEV";
  }
  idat+=ndata_wp;

  // W-
  getline(fWm,line);
  for (int i = 0; i < ndata_wp; i++)
  {
    getline(fWm,line);
    istringstream lstream(line);

    lstream >> etaavg >> etamin >> etamax;
    //binsize=1./(etamax-etamin);
    fKin1[idat+i] = etaavg;     // eta
    fKin2[idat+i] = MW2;        // W mass squared
    fKin3[idat+i] = s;          // sqrt(s)

    lstream >> fData[idat+i];
    fData[idat+i] *= pb2fb;
    //    fData[idat+i] *= binsize;
    lstream >> fStat[idat+i];
    fStat[idat+i] *= pb2fb;
    //    fStat[idat+i] *= binsize;
    lstream >> totsys[idat+i];
    totsys[idat+i] *= pb2fb;
    //    totsys[idat+i] *= binsize;

    lstream >> fSys[idat+i][fNData].add;  // Beam uncertainty
    fSys[idat+i][fNData].add *= pb2fb;
    fSys[idat+i][fNData].mult = fSys[idat+i][fNData].add/fData[idat+i]*1e2;
    fSys[idat+i][fNData].type = MULT;
    fSys[idat+i][fNData].name = "LHCBBEAM7TEV";

    lstream >> fSys[idat+i][fNData+1].add;  // Lumi uncertainty
    fSys[idat+i][fNData+1].add *= pb2fb;
    fSys[idat+i][fNData+1].mult = fSys[idat+i][fNData+1].add/fData[idat+i]*1e2;
    fSys[idat+i][fNData+1].type = MULT;
    fSys[idat+i][fNData+1].name = "LHCBLUMI7TEV";
  }
  idat+=ndata_wm;

  for(int i=0; i<fNData; i++)
    {
      cout << totsys[i] << endl;
    }
  
  // Reading Covariance Matrix
  for (int i = 0; i < fNData; i++) {
    for (int j = 0; j < i+1; j++) {             // read only lower triangle
      fCorr >> inmat[j][i];
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
      covmat[i][j]=inmat[i][j]*totsys[i]*totsys[j];
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
    for (int l = 0; l < fNSys-2; l++)
    {
      fSys[i][l].add  = syscor[i][l];
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
  
  fZ.close();
  fWp.close();
  fWm.close();
  fCorr.close();
}
