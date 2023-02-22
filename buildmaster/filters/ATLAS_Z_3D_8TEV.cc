/*
Name_exp  : ATLAS_Z_3D_8TEV
Reference : Measurement of the triple-differential cross section for the 
            Drell–Yan Z at 8 TeV with ATLAS
ArXiv     : arXiv:1710.05167
Published : JHEP12(2017)059
Hepdata   : https://www.hepdata.net/record/ins1630886

The data are presented in bins of invariant mass, absolute dilepton rapidity, 
|yll|, and the angular variable cosθ* between the outgoing lepton and the 
incoming quark in the Collins–Soper frame. There are two data sets: one in
which the measurements are performed in the range |yll|<2.4 in the muon 
and electron channel (and are subsequently combined); and one in which
the measurements are extended to |yll|<3.6 only in the electron channel. 
Otherwise, the measurement is performed in seven bins of mll from 46 GeV to 
200 GeV with edges set at 66, 80, 91, 102, 116, and 150 GeV; 12 equidistant 
bins of |yll| from 0 to 2.4; and bins of cosθ∗ from −1 to +1, separated at 
−0.7, −0.4, 0.0, +0.4, +0.7 giving 6 bins. In total, 504 measurement bins are 
used for the central rapidity electron and muon channel measurements.
For the high rapidity electron channel the measurement is restricted to the 5 
invariant mass bins in the region 66 < mll < 150 GeV. The |yll| region 
measured in this channel ranges from 1.2 to 3.6 in 5 bins with boundaries 
at 1.6, 2.0, 2.4, 2.8. The cos θ∗ binning is identical to the binning of the 
central analyses. A total of 150 measurement bins is used in this channel.

Breakdown of uncertainties:
- statistical uncertainty (1)
- correlated systematic uncertainties (331)
- uncorrelated systematic uncertainties (1)
- luminosity uncertainty (1)
*/

#include "ATLAS_Z_3D_8TEV.h"

//Central rapidity measurement (electron and muon combined)
void ATLAS_Z_3D_EMU_CRAP_8TEVFilter::ReadData()
{
  // Opening files
  fstream f;
  
  const double convfact = 1000.; //Conversion factor from pb to fb
  const int ncorrsys = 331;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/ATLAS_Z_3D_8TEV/Table4.csv";
  f.open(datafile.str().c_str(), ios::in);
  
  if (f.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  string line;
  double ddum;
  double sys;
  char comma;
  
  //Skip 15 header lines
  for (int i = 0; i < 15; i++)
    {
      getline(f, line);
    }
  
  for (int i = 0; i < fNData; i++)
    {
      getline(f, line);
      istringstream lstream(line);
      double ctmin, ctmax;
      double mlmin, mlmax;
      double binfactor;

      lstream >> fKin3[i] >> comma //costhetastar
	      >> ctmin    >> comma
	      >> ctmax    >> comma
	      >> fKin1[i] >> comma //dilepton rapidity
	      >> ddum     >> comma
	      >> ddum     >> comma
	      >> fKin2[i] >> comma //dilepton mass [GeV]
	      >> mlmin    >> comma
	      >> mlmax    >> comma
	      >> fData[i] >> comma //3D cross section [pb]
	      >> fStat[i] >> comma //statistical uncertainty
	      >> ddum;
      
      //fKin2[i] = mass*mass;
      //fKin3[i] = 8000.;            //c.m. energy [GeV]
      
      //Conversion to [fb] and normalistaion to bin width
      binfactor = (mlmax-mlmin)*0.4;
      fData[i] *= convfact*binfactor;        
      fStat[i] *= convfact*binfactor;        

      //Correlated uncertainty
      for (int isys = 0; isys < ncorrsys; isys++)
        {
	  lstream >> comma>> sys >> comma >> ddum;
	  fSys[i][isys].add = sys * convfact * binfactor;
	  fSys[i][isys].mult = fSys[i][isys].add * 1e2 / fData[i];
	  fSys[i][isys].type = MULT;
	  fSys[i][isys].name = "CORR";
        }
      
      //Uncorrelated uncertainty
      lstream >> comma >> sys >> comma >> ddum;
      fSys[i][ncorrsys].add = sys * convfact * binfactor;
      fSys[i][ncorrsys].mult = fSys[i][ncorrsys].add * 1e2 / fData[i];
      fSys[i][ncorrsys].type = MULT;
      fSys[i][ncorrsys].name = "UNCORR";

      //Luminosity uncertainty
      fSys[i][ncorrsys + 1].mult = 1.8; // in percent
      fSys[i][ncorrsys + 1].add = fSys[i][ncorrsys + 1].mult * fData[i] / 100;
      fSys[i][ncorrsys + 1].type = MULT;
      fSys[i][ncorrsys + 1].name = "ATLASLUMI12";
    }

    f.close();
}

//Forward rapidity measurement (electron only)
void ATLAS_Z_3D_ELE_HRAP_8TEVFilter::ReadData()
{
  // Opening files
  fstream f;
  
  const double convfact = 1000.; //Conversion factor from pb to fb
  const int ncorrsys = 331;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/ATLAS_Z_3D_8TEV/Table4.csv";
  f.open(datafile.str().c_str(), ios::in);
  
  if (f.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  string line;
  double ddum;
  double sys;
  char comma;
  
  //Skip 15 header lines
  for (int i = 0; i < 15; i++)
    {
      getline(f, line);
    }

  //Skip 458 data points for central rapidity
  for(int i = 0; i<458; i++)
    {
      getline(f, line);
    }
  
  for (int i = 0; i < fNData; i++)
    {
      getline(f, line);
      istringstream lstream(line);
      double ctmin, ctmax;
      double mlmin, mlmax;
      double binfactor;

      lstream >> fKin3[i] >> comma //costhetastar
	      >> ctmin    >> comma
	      >> ctmax    >> comma
	      >> fKin1[i] >> comma //dilepton rapidity
	      >> ddum     >> comma
	      >> ddum     >> comma
	      >> fKin2[i] >> comma //dilepton mass [GeV]
	      >> mlmin    >> comma
	      >> mlmax    >> comma
	      >> fData[i] >> comma //3D cross section [pb]
	      >> fStat[i] >> comma //statistical uncertainty
	      >> ddum;
      
      //fKin2[i] = mass*mass;
      //fKin3[i] = 8000.;            //c.m. energy [GeV]
      
      //Conversion to [fb]
      binfactor = (mlmax-mlmin)*0.4;
      fData[i] = fData[i] * convfact * binfactor;        
      fStat[i] = fStat[i] * convfact * binfactor;        

      //Correlated uncertainty
      for (int isys = 0; isys < ncorrsys; isys++)
        {
	  lstream >> comma>> sys >> comma >> ddum;
	  fSys[i][isys].add = sys * convfact * binfactor;
	  fSys[i][isys].mult = fSys[i][isys].add * 1e2 / fData[i];
	  fSys[i][isys].type = MULT;
	  fSys[i][isys].name = "CORR";
        }
      
      //Uncorrelated uncertainty
      lstream >> comma >> sys >> comma >> ddum;
      fSys[i][ncorrsys].add = sys * convfact * binfactor;
      fSys[i][ncorrsys].mult = fSys[i][ncorrsys].add * 1e2 / fData[i];
      fSys[i][ncorrsys].type = MULT;
      fSys[i][ncorrsys].name = "UNCORR";

      //Luminosity uncertainty
      fSys[i][ncorrsys + 1].mult = 1.8; // in percent
      fSys[i][ncorrsys + 1].add = fSys[i][ncorrsys + 1].mult * fData[i] / 100;
      fSys[i][ncorrsys + 1].type = MULT;
      fSys[i][ncorrsys + 1].name = "ATLASLUMI12";
    }

    f.close();
}
