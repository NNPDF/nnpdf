/*
Name_exp : ATLAS_WMU_8TEV
Reference: Measurement of differential cross sections and 𝑊+/𝑊−
           cross-section ratios for 𝑊 boson production in association 
           with jets at 𝑠√=8 TeV with the ATLAS detector
ArXiv    : arxiv:1904.05631
Published: Eur.Phys.J.C 79 (2019) 9, 760
Hepdata  : https://www.hepdata.net/record/ins1729240

Measurements of the W+→μ+ν and W−→μ−ν¯ cross-sections and the associated 
charge asymmetry as a function of the absolute pseudorapidity of the decay muon.
The data were collected in proton-proton collisions at a centre-of-mass energy 
of 8 TeV with the ATLAS experiment at the LHC and correspond to a total 
integrated luminosity of 20.2 fb−1. The integrated cross-sections for W+→μ+ν  
and W−→μ−ν¯ production are measured in a fiducial phase space defined 
at the particle level by requiring the muon transverse momentum to be greater 
than 25 GeV  and the neutrino transverse momentum to be greater than 25 GeV. 
The absolute muon pseudorapidity is required to be less than 2.4. 
The W boson transverse mass must be at least 40 GeV.
The differential cross-sections and charge asymmetry are measured in the same 
fiducial phase space as for the integrated measurement. 
These are measured in 11 bins of absolute muon pseudorapidity between 0 and 
2.4 with bin edges at 0, 0.21, 0.42, 0.63, 0.84, 1.05, 1.37, 1.52, 1.74, 1.95, 
2.18, and 2.4. The carge asymmetry is not implemented. Systematic uncertainties
are correlated between W+ and W-. An extra systematic uncertainty due to the
ATLAS luminosity at 8 TeV is implemented.
*/

#include "ATLAS_WMU_8TEV.h"

//1)W+ distribution

void ATLAS_WP_MU_8TEVFilter::ReadData()
{
  fstream f1;

  //Central values and uncertainties
  stringstream datafile_sys("");
  datafile_sys << dataPath()
	   << "rawdata/ATLAS_WMU_8TEV/HEPData-ins1729240-v1-Table_3.csv";
  f1.open(datafile_sys.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile_sys.str() << endl;
      exit(-1);
    }


  //Read central value
  string line;
  for(int i=0; i<12; i++)
    {
      getline(f1,line);
    }

  double ddum;
  char comma;

  cout << fNData << endl;
  
  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      lstream >> fKin1[i]        >> comma
	      >> ddum            >> comma
	      >> ddum            >> comma
	      >> fData[i]        >> comma //in [pb]
	      >> fStat[i]        >> comma //in percentage
	      >> ddum;

      fData[i] *= 1000.; //Convert to [fb]
      fKin2[i] = MW*MW;
      fKin3[i] = 8000.; // GeV
      fStat[i] = fStat[i]/100. * fData[i]; //convert to absolute value
      
      //Systematic uncertainties
      for(int j=0; j<fNSys-1; j++)
	{
	  lstream >> comma >> fSys[i][j].mult
		  >> comma >> ddum;
	  fSys[i][j].add = fSys[i][j].mult/100. * fData[i];
	  fSys[i][j].type = MULT;
	}
      
      fSys[0][j].name = "ATLASWMU_ETmiss";
      fSys[1][j].name = "ATLASWMU_MuonReco";
      fSys[2][j].name = "Background";
      fSys[3][j].name = "UNCORR";
      fSys[4][j].name = "Modelling";
	
      //Luminosity uncertainty
      fSys[i][5].mult = 1.9;
      fSys[i][5].add = fSys[i][5].mult/100. * fData[i];
      fSys[i][5].type = MULT;
      fSys[i][5].name = "ATLASLUMI12";
      
    }
  
  f1.close();
  
} 

//2)W- distribution

void ATLAS_WM_MU_8TEVFilter::ReadData()
{
  fstream f1;

  //Central values and uncertainties
  stringstream datafile_sys("");
  datafile_sys << dataPath()
	   << "rawdata/ATLAS_WMU_8TEV/HEPData-ins1729240-v1-Table_3.csv";
  f1.open(datafile_sys.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile_sys.str() << endl;
      exit(-1);
    }

  //Read central value
  string line;
  for(int i=0; i<28; i++)
    {
      getline(f1,line);
    }

  double ddum;
  char comma;

  cout << fNData << endl;
  
  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      lstream >> fKin1[i]        >> comma
	      >> ddum            >> comma
	      >> ddum            >> comma
	      >> fData[i]        >> comma //[pb]
	      >> fStat[i]        >> comma //in percentage
	      >> ddum;
      
      fData[i] *= 1000.; //Convert to [fb]  
      fKin2[i] = MW*MW;
      fKin3[i] = 8000.; // GeV
      fStat[i] = fStat[i]/100. * fData[i]; //convert to absolute value
      
      //Systematic uncertainties
      for(int j=0; j<fNSys-1; j++)
	{
	  lstream >> comma >> fSys[i][j].mult
		  >> comma >> ddum;
	  fSys[i][j].add = fSys[i][j].mult/100. * fData[i];
	  fSys[i][j].type = MULT;
	}

      fSys[0][j].name = "ATLASWMU_ETmiss";
      fSys[1][j].name = "ATLASWMU_MuonReco";
      fSys[2][j].name = "Background";
      fSys[3][j].name = "UNCORR";
      fSys[4][j].name = "Modelling";
      
      //Luminosity uncertainty
      fSys[i][5].mult = 1.9;
      fSys[i][5].add = fSys[i][5].mult/100. * fData[i];
      fSys[i][5].type = MULT;
      fSys[i][5].name = "ATLASLUMI12";
      
    }
  
  f1.close();
  
} 









































