/*
  Z boson production at the LHCb experiment in proton proton collision 
  at 13 TeV, The analysis uses a dataset corresponding to an integrated 
  luminosity of 294 +- 11 pb^(-1). Event considered are those where Z decays 
  either to a dimuon or a dielecron final state Z -> mu mubar and Z -> e ebar 
  Both normalized and unnormalized differential distributions are implemented 
  for both channels.

  Data from paper [1607.06495] 
  https://arxiv.org/abs/1607.06495
*/

#include "LHCb.h"

//1) Distribution differential in Z boson rapidity in di-electron channel

void LHCB_Z_13TEV_DIELECTRONFilter::ReadData()
{

  // Opening files
  fstream rZ, rCorr;

  // rapidity distribution
  stringstream DataFileZ("");
  DataFileZ << dataPath() << "rawdata/" << fSetName 
	    << "/LHCBZ13TEV_zrap.data";
  rZ.open(DataFileZ.str().c_str(), ios::in);

  if (rZ.fail()) {
    cerr << "Error opening data file " << DataFileZ.str() << endl;
    exit(-1);
  }

  // correlation matrix
  stringstream DataFileCorr("");
  DataFileCorr << dataPath() << "rawdata/" << fSetName 
	       << "/LHCBZ13TEV_corrmat.data";
  rCorr.open(DataFileCorr.str().c_str(), ios::in);

  if (rCorr.fail()) {
    cerr << "Error opening data file " << DataFileCorr.str() << endl;
    exit(-1);
  }

  // Starting filter
  const int ndata_z  = 17;
  const double pb2fb = 1000.;  // Must multiply from pb to fb  
  double MZ2 = pow(MZ,2.0);
  double s = 13000;
  string line;

  std::vector<double> totsys(fNData);        
  double etaavg, etamin, etamax;

  // Z rapidity distribution
  for (int i = 0; i < 3; i++)
    getline(rZ,line);
  
  for (int i = 0; i < ndata_z; i++)
  {
    getline(rZ,line);                  
    istringstream lstream(line); 

    lstream >> etaavg >> etamin >> etamax;
    fKin1[i] = etaavg;          // eta, central value of the bin
    fKin2[i] = MZ2;             // Z mass squared
    fKin3[i] = s;               // sqrt(s)

    lstream >> fData[i];        
    fData[i] *= pb2fb;          
    lstream >> fStat[i];        
    fStat[i] *= pb2fb;
    lstream >> totsys[i];       
    totsys[i] *= pb2fb;

    lstream >> fSys[i][fNSys-1].add;  // Lumi uncertainty
    fSys[i][fNSys-1].add *= pb2fb;
    fSys[i][fNSys-1].mult = fSys[i][fNSys-1].add/fData[i]*1e2;
    fSys[i][fNSys-1].type = MULT;
    fSys[i][fNSys-1].name = "LHCBLUMI13TEV";
  }
  
  // Defining covariance matrix
  double** covmat = new double*[fNData];
  for (int i = 0; i < fNData; i++) 
    covmat[i] = new double[fNData];
 
  // Reading Covariance Matrix
  for (int i = 0; i < fNData; i++)
    {
      for (int j = 0; j < i+1; j++) 
	{
	  rCorr >> covmat[j][i];
	  covmat[j][i]=covmat[j][i]*totsys[j]*totsys[i];
	  covmat[i][j]=covmat[j][i];
	}
    }
  
  // Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];
  
  if(!genArtSys(fNData,covmat,syscor))
    {
      cerr << " in " << fSetName 
	   << " : cannot generate artificial systematics" << endl;
      exit(-1);
    }
  
  // Copy the artificial systematics in the fSys matrix
  for (int i = 0; i < fNData; i++)
    for (int l = 0; l < fNSys-1; l++)
    {
      fSys[i][l].add  = syscor[i][l];
      fSys[i][l].mult = fSys[i][l].add/fData[i]*1e2;
      fSys[i][l].type = MULT;  
      fSys[i][l].name = "CORR";
    }
  
  rZ.close();
  rCorr.close();
  
  for(int i = 0; i < fNData; i++) 
    delete[] covmat[i];
  delete[] covmat;
  
  for(int i = 0; i < fNData; i++) 
    delete[] syscor[i];
  delete[] syscor;
  
}

//2) Distribution differential in Z boson rapidity in di-muon channel

void LHCB_Z_13TEV_DIMUONFilter::ReadData()
{

  // Opening files
  fstream rZ, rCorr;

  // rapidity distribution
  stringstream DataFileZ("");
  DataFileZ << dataPath() << "rawdata/" << fSetName 
	    << "/LHCBZ13TEV_zrap_dimuon.data";
  rZ.open(DataFileZ.str().c_str(), ios::in);

  if (rZ.fail()) {
    cerr << "Error opening data file " << DataFileZ.str() << endl;
    exit(-1);
  }

  // correlation matrix
  stringstream DataFileCorr("");
  DataFileCorr << dataPath() << "rawdata/" << fSetName 
	       << "/LHCBZ13TEV_corrmat_dimuon.data";
  rCorr.open(DataFileCorr.str().c_str(), ios::in);

  if (rCorr.fail()) {
    cerr << "Error opening data file " << DataFileCorr.str() << endl;
    exit(-1);
  }

  // Starting filter
  const int ndata_z  = 18;
  const double pb2fb = 1000.;   
  double MZ2 = pow(MZ,2.0);
  double s = 13000;
  string line;

  std::vector<double> totsys(fNData);        
  double etaavg, etamin, etamax;

  // Z rapidity distribution
  for (int i = 0; i < 3; i++)
    getline(rZ,line);
  
  for (int i = 0; i < ndata_z; i++)
    {
      getline(rZ,line);        
      istringstream lstream(line); 
      
      lstream >> etaavg >> etamin >> etamax;
      fKin1[i] = etaavg;          // eta, central value of the bin
      fKin2[i] = MZ2;             // Z mass squared
      fKin3[i] = s;               // sqrt(s)
      
      lstream >> fData[i];        
      fData[i] *= pb2fb;          
      lstream >> fStat[i];        
      fStat[i] *= pb2fb;
      lstream >> totsys[i];       
      totsys[i] *= pb2fb;
      
      lstream >> fSys[i][fNSys-1].add;  // Lumi uncertainty
      fSys[i][fNSys-1].add *= pb2fb;
      fSys[i][fNSys-1].mult = fSys[i][fNSys-1].add/fData[i]*1e2;
      fSys[i][fNSys-1].type = MULT;
      fSys[i][fNSys-1].name = "LHCBLUMI13TEV";
    }
  
  // Defining covariance matrix
  double** covmat = new double*[fNData];
  for (int i = 0; i < fNData; i++) 
    covmat[i] = new double[fNData];
  
  // Reading Covariance Matrix
  for (int i = 0; i < fNData; i++)
    { 
      for (int j = 0; j < i+1; j++) 
	{
	  rCorr >> covmat[j][i];
	  covmat[j][i]=covmat[j][i]*totsys[j]*totsys[i];
	  covmat[i][j]=covmat[j][i];
	}
    }
  
  // Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];
  
  if(!genArtSys(fNData,covmat,syscor))
   {
     cerr << " in " << fSetName 
	  << " : cannot generate artificial systematics" << endl;
     exit(-1);
   }
  
  // Copy the artificial systematics in the fSys matrix
  for (int i = 0; i < fNData; i++)
    for (int l = 0; l < fNSys-1; l++)
      {
	fSys[i][l].add  = syscor[i][l];
	fSys[i][l].mult = fSys[i][l].add/fData[i]*1e2;
	fSys[i][l].type = MULT;  
	fSys[i][l].name = "CORR";
      }
  
  rZ.close();
  rCorr.close();
  
  for(int i = 0; i < fNData; i++) 
    delete[] covmat[i];
  delete[] covmat;
  
  for(int i = 0; i < fNData; i++) 
    delete[] syscor[i];
  delete[] syscor;
  
}

