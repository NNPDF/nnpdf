/********************* NEW DATA in NNPDF3.1 **************************
*
* LHCb 2fb^{-1}, 8TeV
*
* W/Z LHCb experiment
*
*/

#include "LHCb.h"

void LHCBWZMU8TEVFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3, f4;

  stringstream datafile1("");
  datafile1 << dataPath() << "rawdata/" << fSetName << "/LHCBWZMU8TEV_wplrap.data";
  f1.open(datafile1.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile1.str() << endl;
    exit(-1);
  }

  stringstream datafile2("");
  datafile2 << dataPath() << "rawdata/" << fSetName << "/LHCBWZMU8TEV_wmlrap.data";
  f2.open(datafile2.str().c_str(), ios::in);

  if (f2.fail()) {
    cerr << "Error opening data file " << datafile2.str() << endl;
    exit(-1);
  }

  stringstream datafile3("");
  datafile2 << dataPath() << "rawdata/" << fSetName << "/LHCBWZMU8TEV_zrap.data";
  f2.open(datafile2.str().c_str(), ios::in);

  if (f2.fail()) {
    cerr << "Error opening data file " << datafile2.str() << endl;
    exit(-1);
  }

  stringstream datafile4("");
  datafile3 << dataPath() << "rawdata/" << fSetName << "/LHCBWZMU8TEV_covmat.data";
  f3.open(datafile3.str().c_str(), ios::in);

  if (f3.fail()) {
    cerr << "Error opening data file " << datafile3.str() << endl;
    exit(-1);
  }

  // Starting filter
  const int ndata_z  = 18;
  const int ndata_wp =  8;
  const int ndata_wm =  8;
  double binsize;  // Must multiply by binsize to match inclusive bin-by-bin data
  const double convfac = 1000.; // Must multiply from pb to fb
  double MW2 = pow(MW,2.0);
  double MZ2 = pow(MZ,2.0);
  double s   = 8000;
  string line;

  double totsys[fNData];
  double etaavg, etamin, etamax;
  int idat = 0;

  cout << "********** WARNING: Converting pb to fb to match ApplGrid output ********" << endl;

  // Z
  getline(f1,line);
  for (int i = 0; i < ndata_z; i++)
  {
    getline(f1,line);
    istringstream lstream(line);

    lstream >> etaavg >> etamin >> etamax;
    binsize=1./(etamax-etamin);
    fKin1[idat+i] = etaavg;     // eta
    fKin2[idat+i] = MZ2;        // Z mass squared
    fKin3[idat+i] = s;          // sqrt(s)

    lstream >> fData[idat+i];
    fData[idat+i] *= convfac;
    //    fData[idat+i] *= binsize;
    lstream >> fStat[idat+i];
    fStat[idat+i] *= convfac;
    //    fStat[idat+i] *= binsize;
    lstream >> totsys[idat+i];
    totsys[idat+i] *= convfac;
    //    totsys[idat+i] *= binsize;

    lstream >> fSys[idat+i][fNData].add;  // Beam uncertainty
    fSys[idat+i][fNData].mult = fSys[idat+i][0].add/fData[i]*1e2;
    fSys[idat+i][fNData].type = MULT;
    fSys[idat+i][fNData].name = "LHCBBEAM11";

    lstream >> fSys[idat+i][fNData+1].add;  // Lumi uncertainty
    fSys[idat+i][fNData+1].mult = fSys[idat+i][0].add/fData[i]*1e2;
    fSys[idat+i][fNData+1].type = MULT;
    fSys[idat+i][fNData+1].name = "LHCBLUMI11";
  }
  idat+=ndata_z;

  // W+
  getline(f2,line);
  for (int i = 0; i < ndata_wp; i++)
  {
    getline(f1,line);
    istringstream lstream(line);

    lstream >> etaavg >> etamin >> etamax;
    binsize=1./(etamax-etamin);
    fKin1[idat+i] = etaavg;     // eta
    fKin2[idat+i] = MZ2;        // Z mass squared
    fKin3[idat+i] = s;          // sqrt(s)

    lstream >> fData[idat+i];
    fData[idat+i] *= convfac;
    //    fData[idat+i] *= binsize;
    lstream >> fStat[idat+i];
    fStat[idat+i] *= convfac;
    //    fStat[idat+i] *= binsize;
    lstream >> totsys[idat+i];
    totsys[idat+i] *= convfac;
    //    totsys[idat+i] *= binsize;

    lstream >> fSys[idat+i][fNData].add;  // Beam uncertainty
    fSys[idat+i][fNData].mult = fSys[idat+i][fNData].add/fData[i]*1e2;
    fSys[idat+i][fNData].type = MULT;
    fSys[idat+i][fNData].name = "LHCBBEAM11";

    lstream >> fSys[idat+i][fNData+1].add;  // Lumi uncertainty
    fSys[idat+i][fNData+1].mult = fSys[idat+i][fNData+1].add/fData[i]*1e2;
    fSys[idat+i][fNData+1].type = MULT;
    fSys[idat+i][fNData+1].name = "LHCBLUMI11";
  }
  idat+=ndata_wp;

  // W-
  getline(f3,line);
  for (int i = 0; i < ndata_wp; i++)
  {
    getline(f3,line);
    istringstream lstream(line);

    lstream >> etaavg >> etamin >> etamax;
    binsize=1./(etamax-etamin);
    fKin1[idat+i] = etaavg;     // eta
    fKin2[idat+i] = MZ2;        // Z mass squared
    fKin3[idat+i] = s;          // sqrt(s)

    lstream >> fData[idat+i];
    fData[idat+i] *= convfac;
    //    fData[idat+i] *= binsize;
    lstream >> fStat[idat+i];
    fStat[idat+i] *= convfac;
    //    fStat[idat+i] *= binsize;
    lstream >> totsys[idat+i];
    totsys[idat+i] *= convfac;
    //    totsys[idat+i] *= binsize;

    lstream >> fSys[idat+i][fNData].add;  // Beam uncertainty
    fSys[idat+i][fNData].mult = fSys[idat+i][fNData].add/fData[i]*1e2;
    fSys[idat+i][fNData].type = MULT;
    fSys[idat+i][fNData].name = "LHCBBEAM11";

    lstream >> fSys[idat+i][fNData+1].add;  // Lumi uncertainty
    fSys[idat+i][fNData+1].mult = fSys[idat+i][fNData+1].add/fData[i]*1e2;
    fSys[idat+i][fNData+1].type = MULT;
    fSys[idat+i][fNData+1].name = "LHCBLUMI11";
  }
  idat+=ndata_wm;

  //Reading correlation matrix
  //Format of table means we need to read all 15 points (W+,W- and Z)
  double inmat[fNData][fNData];
  for (int i = 0; i < fNData; i++)
  {
    getline(f4,line);
    istringstream lstream(line);
    for (int j = 0; j < i+1; j++)      //Only lower diagonal in file
    {
      lstream >> inmat[i][j];
      inmat[j][i]=inmat[i][j];         //symmetrize
    }
  }

  //Put corrmat entries into correct order
  double** covmat = new double*[fNData];
  for(int i = 0; i < fNData-1; i++)
  {
    covmat[i] = new double[fNData];
    for(int j = 0; j < fNData-1; j++)
      covmat[i][j]=inmat[(2*i)%15][(2*j)%15]*totsys[i]*totsys[j];
    covmat[i][15]=inmat[(2*i)%15][15]*totsys[i]*totsys[15];
  }
  covmat[15] = new double[fNData];
  for(int j = 0; j < fNData-1;j++)
  	covmat[15][j]=inmat[15][(2*j)%15]*totsys[15]*totsys[j];
  covmat[15][15]=inmat[15][15]*totsys[15]*totsys[15];

  // Now generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
   {
     cerr << " in " << fSetName << endl;
     exit(-1);
   }

  for (int i = 0; i < fNData; i++)
    for (int l = 0; l < fNData; l++)
    {
      fSys[i][l+1].add = syscor[i][l];
      fSys[i][l+1].mult = fSys[i][l+1].add*100/fData[i];
      fSys[i][l+1].type = ADD;
      fSys[i][l+1].name = "CORR";
    }
  f1.close();
  f2.close();
  f3.close();
}
