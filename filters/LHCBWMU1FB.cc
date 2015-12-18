/********************* NEW DATA in NNPDF3.1 **************************
*
* LHCb 1fb^{-1}
*
* W production > mu nu_mu data from the LHCb experiment
*
* Final data from the LHCb preprint: 1408.4354
* Luminosity uncertainty is a 1.71% in
* all data points and is quoted separately from other sources
* of systematic uncertainty
*/

#include "LHCb.h"

void LHCBWMU1FBFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/lhcb_wplrap.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream datafile2("");
  datafile2 << dataPath() << "rawdata/"
  << fSetName << "/lhcb_wmlrap.data";
  f2.open(datafile2.str().c_str(), ios::in);

  if (f2.fail()) {
    cerr << "Error opening data file " << datafile2.str() << endl;
    exit(-1);
  }

  stringstream datafile3("");
  datafile3 << dataPath() << "rawdata/"
   << fSetName << "/lhcb_covmat.data";
  f3.open(datafile3.str().c_str(), ios::in);

  if (f3.fail()) {
    cerr << "Error opening data file " << datafile3.str() << endl;
    exit(-1);
  }

  // Starting filter
  const int ndata_wp = 8;
  const int ndata_wm = 8;
  double binsize;  // Must multiply by binsize to match inclusive bin-by-bin data
  const double convfac = 1000.; // Must multiply from pb to fb
  double MW2 = pow(MW,2.0);
  double s = 7000;
  string line;

  double totsys[fNData];
  double etamin, etamax;
  int idat = 0;

  cout << "********** WARNING: Converting pb to fb to match ApplGrid output ********" << endl;

  // W+
  for (int i = 0; i < ndata_wp; i++)
  {
    getline(f1,line);
    istringstream lstream(line);

    lstream >> etamin >> etamax;
    binsize=1./(etamax-etamin);
    fKin1[idat+i] = (etamax + etamin)*0.5;     //eta
    fKin2[idat+i] = MW2;                       //mass W squared
    fKin3[idat+i] = s;                         //sqrt(s)

    lstream >> fData[idat+i];
    fData[idat+i] *= convfac;
    fData[idat+i] *= binsize;
    lstream >> fStat[idat+i];
    fStat[idat+i] *= convfac;
    fStat[idat+i] *= binsize;
    lstream >> totsys[idat+i];
    totsys[idat+i] *= convfac;
    totsys[idat+i] *= binsize;

    fSys[idat+i][0].mult = 1.71;                     //luminosity uncertainty of 1.71%
    fSys[idat+i][0].add = fSys[idat+i][0].mult*fData[i]*1e-2;
    fSys[idat+i][0].type = MULT;
    fSys[idat+i][0].name = "CORR";
  }
  idat+=ndata_wp;

  // W-
  for (int i = 0; i < ndata_wm; i++)
  {
    getline(f2,line);
    istringstream lstream(line);

    lstream >> etamin >> etamax;
    binsize=1./(etamax-etamin);
    fKin1[idat+i] = (etamax + etamin)*0.5;     //eta
    fKin2[idat+i] = MW2;                       //mass W squared
    fKin3[idat+i] = s;                         //sqrt(s)

    lstream >> fData[idat+i];
    fData[idat+i] *= convfac;
    fData[idat+i] *= binsize;
    lstream >> fStat[idat+i];
    fStat[idat+i] *= convfac;
    fStat[idat+i] *= binsize;
    lstream >> totsys[idat+i];
    totsys[idat+i] *= convfac;
    totsys[idat+i] *= binsize;

    fSys[idat+i][0].mult = 1.71;                     //luminosity uncertainty of 1.71%
    fSys[idat+i][0].add = fSys[idat+i][0].mult*fData[i]*1e-2;
    fSys[idat+i][0].type = MULT;
    fSys[idat+i][0].name = "CORR";
  }
  idat+=ndata_wm;

  //Reading correlation matrix
  //Format of table means we need to read all 15 points (W+,W- and Z)
  double inmat[16][16];
  for (int i = 0; i < 16; i++)
  {
    getline(f3,line);
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
