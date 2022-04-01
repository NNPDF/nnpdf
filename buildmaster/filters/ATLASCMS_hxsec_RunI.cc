/*
References: 
1) For the first 20 signal strenghts
   [arXiv:1606.02266]
   Measurements of the Higgs boson production and decay rates and constraints 
   on its couplings from a combined ATLAS and CMS analysis of the LHC pp 
   collision data at √s = 7 and 8 TeV
   JHEP 1608 (2016) 045
   https://hepdata.net/record/ins1468068
2) For the last 2 signal strengths
   [arXiv:1507.04548]
   Measurements of the Higgs boson production and decay rates and coupling 
   strengths using pp collision data at √s = 7 and 8 TeV in the ATLAS experiment
   Eur.Phys.J. C76 (2016) no.1, 6

Best fit values of $\sigma_i \cdot \mathrm{B}^f$ relative to their SM 
prediction for each specific channel $i \to H\to f$, as obtained from the 
generic parameterisation with 23 parameters for the combination of the ATLAS 
and CMS measurements, using the $\sqrt{s}$=7 and 8 TeV data. 
The order of the signal strengths is as follows
ggF gamma gamma
    Z Z
    W W
    tau tau
VBF gamma gamma
    Z Z
    W W
    tau tau
Wh  gamma gamma
    W W
    tau tau
    b b
Zh  gamma gamma
    W W
    tau tau
    b b
tth gamma gamma
    W W
    tau tau
    b b
pp  Z gamm
    mu mu
The results are implemented their breakdown into statistical and systematic 
components. The missing values are either not measured with a meaningful 
precision and therefore not quoted, in the case of the $H\to ZZ$ decay channel 
for the $WH$, $ZH$, and $ttH$ production processes, or not measured at all and 
therefore fixed to their corresponding SM predictions, in the case of the 
$H\to bb$ decay mode for the $gg\mathrm{F}$ and VBF production processes.
*/
#include "ATLASCMS_hxsec_RunI.h"

void ATLASCMS_hxsec_RunIFilter::ReadData()
{
  fstream f1;
  fstream f2;

  //Central values, statistical and systematic ucnertainties
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/ATLASCMS_hxsec_RunI/HEPData-ins1468068-v1-Table_2.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Correlations between the 20 data points
  stringstream datafile_corr("");
  datafile_corr << dataPath()
		<< "rawdata/ATLASCMS_hxsec_RunI/HEPData-ins1468068-v1-Table_18.csv";
  f2.open(datafile_corr.str().c_str(), ios::in);

  if (f2.fail())
    {
      cerr << "Error opening data file " << datafile_corr.str() << endl;
      exit(-1);
    }

  //Read central values, statistical and systematic ucnertainties
  string line;
  for(int i=0; i<12; i++)
    {
      getline(f1,line);
    }
  for(int i=0; i<11; i++)
    {
      getline(f2,line);
    }
  
  string sdum;
  char comma;
  double* Sys = new double[fNData];
  double** corrmat = new double*[fNData-2];
  double** syscor  = new double*[fNData];
  double statR, statL, systR, systL;
  double stmp, dtmp;
  double shift;
  double data_sh;
  
  for(int i=0; i<fNData-2; i++)
    {
      shift=0.;
      getline(f1,line);
      istringstream lstream(line);
      fKin1[i] = 0.;
      fKin2[i] = 0.;
      fKin3[i] = 0.;
      lstream >> sdum  
	      >> comma >> fData[i]
	      >> comma >> statR 
	      >> comma >> statL
	      >> comma >> systR
	      >> comma >> systL;

      statR = statR/fData[i]*100;
      statL = statL/fData[i]*100;
      systR = systR/fData[i]*100;
      systL = systL/fData[i]*100;

      symmetriseErrors(statR,statL,&stmp,&dtmp);
      data_sh = fData[i]*(1.0 + dtmp*0.01);
      fStat[i]=stmp/100.*data_sh;
      shift += dtmp;

      symmetriseErrors(systR,systL,&stmp,&dtmp);
      data_sh = fData[i]*(1.0 + dtmp*0.01);
      Sys[i]=stmp/100*data_sh;
      shift += dtmp;
      
      fData[i]*=(1.0 + shift*0.01);

      corrmat[i] = new double[fNData-2];

      for(int j=0; j<fNData-2; j++)
	{
	  getline(f2,line);
	  istringstream kstream(line);
	  kstream >> sdum >> comma >> corrmat[i][j];
	}
    }

  //Add the missing two data points and generate artificial systematics
  fKin1[20]=0.;
  fKin2[20]=0.;
  fKin3[20]=0.;
  fData[20]=2.7;
  fStat[20]=0.0;
  Sys[20]=4.6;

  fKin1[21]=0.;
  fKin2[21]=0.;
  fKin3[21]=0.;
  fData[21]=-0.7;
  fStat[21]=0.0;
  Sys[21]=3.7;

  double** corrmat_cmp = new double*[fNData];

  for(int i=0; i<fNData; i++)
    {
      corrmat_cmp[i] = new double[fNData];
      syscor[i]  = new double[fNData];
      for(int j=0; j<fNData; j++)
	{
	  if(i<fNData-2 && j<fNData-2)
	    {
	      corrmat_cmp[i][j] = corrmat[i][j];
	    }
	  else
	    {
	      if(i==j)
		corrmat_cmp[i][j]=1.;
	      else
		corrmat_cmp[i][j]=0.;
	    }
	  corrmat_cmp[i][j] *= Sys[i]*Sys[j];
	}
    }

  if(!genArtSys(fNData,corrmat_cmp,syscor))
    {
      throw runtime_error("Couldn't generate artificial systematics for " + fSetName);
    }

  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNSys; j++)
	{
	  fSys[i][j].add  = syscor[i][j];
	  fSys[i][j].mult = fSys[i][j].add*1e2/fData[i];
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = "CORR";
	}
    }

  for(int i=0; i<fNData; i++)
    {
      delete[] syscor[i];
      delete[] corrmat_cmp[i];
    }
  
  for(int i=0; i<fNData-2; i++)
    {
      delete[] corrmat[i];;
    }
  
  delete[] corrmat;
  delete[] syscor;
  delete[] corrmat_cmp;
  delete[] Sys;
  
  f1.close();
  f2.close();

}
