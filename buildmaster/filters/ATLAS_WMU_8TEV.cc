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
2.18, and 2.4. The charge asymmetry is not implemented. Systematic uncertainties
are correlated between W+ and W-. An extra systematic uncertainty due to the
ATLAS luminosity at 8 TeV is implemented.
*/

#include "ATLAS_WMU_8TEV.h"

//Combined W+ and W- distributions

void ATLAS_WMU_8TEVFilter::ReadData()
{
  fstream f1;

  //Central values and uncertainties: W+
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/ATLAS_WMU_8TEV/HEPData-ins1729240-v1-Table_3.csv";
  f1.open(datafile.str().c_str(), ios::in); 

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
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

  for(int i=0; i<fNData/2; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      lstream >> fKin1[i]        >> comma
	      >> ddum            >> comma
	      >> ddum            >> comma
	      >> fData[i];       //in [pb]
	
      fData[i] *= 1000.; //Convert to [fb]
      fKin2[i] = MW*MW;
      fKin3[i] = 8000.;  //GeV

      //Statistical uncertainty
      lstream >> comma >> fStat[i]
	      >> comma >> ddum;
      fStat[i] = fStat[i]/100. * fData[i]; //convert to absolute value
      
      //Systematic uncertainties
      for(int j=0; j<fNSys-1; j++)
	{
	  lstream >> comma >> fSys[i][j].mult
		  >> comma >> ddum;
	  
	  fSys[i][j].add = fSys[i][j].mult/100. * fData[i];
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = "CORR";
	}

      //Luminosity uncertainty
      fSys[i][5].mult = 1.9;
      fSys[i][5].add = fSys[i][5].mult/100. * fData[i];
      fSys[i][5].type = MULT;
      fSys[i][5].name = "ATLASLUMI12";

      //Uncorrelated uncertainties
      fSys[i][1].name = "UNCORR";
      fSys[i][3].name = "UNCORR";
    }

  for(int i=0; i<5; i++)
    {
      getline(f1,line);
    }
  
  for(int i=fNData/2; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      lstream >> fKin1[i]        >> comma
	      >> ddum            >> comma
	      >> ddum            >> comma
	      >> fData[i];       //in [pb]
	
      fData[i] *= 1000.; //Convert to [fb]
      fKin2[i] = MW*MW;
      fKin3[i] = 8000.;  //GeV

      //Statistical uncertainty
      lstream >> comma >> fStat[i]
	      >> comma >> ddum;
      fStat[i] = fStat[i]/100. * fData[i]; //convert to absolute value
      
      //Systematic uncertainties
      for(int j=0; j<fNSys-1; j++)
	{
	  lstream >> comma >> fSys[i][j].mult
		  >> comma >> ddum;
	  
	  fSys[i][j].add = fSys[i][j].mult/100. * fData[i];
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = "CORR";
	}

      //Luminosity uncertainty
      fSys[i][5].mult = 1.9;
      fSys[i][5].add = fSys[i][5].mult/100. * fData[i];
      fSys[i][5].type = MULT;
      fSys[i][5].name = "ATLASLUMI12";

      //Uncorrelated uncertainties
      fSys[i][1].name = "UNCORR";
      fSys[i][3].name = "UNCORR";
    }
  
  f1.close();
}


/*
void ATLAS_WMU_8TEVFilter::ReadData()
{
  fstream f1, f2;

  //Central values and uncertainties: W+
  stringstream datafile_WP("");
  datafile_WP << dataPath()
	   << "rawdata/ATLAS_WMU_8TEV/HEPData-ins1729240-v1-Table_5a.csv";
  f1.open(datafile_WP.str().c_str(), ios::in); 

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile_WP.str() << endl;
      exit(-1);
    }

  //Central values and uncertainties: W-
  stringstream datafile_WM("");
  datafile_WM << dataPath()
	   << "rawdata/ATLAS_WMU_8TEV/HEPData-ins1729240-v1-Table_5b.csv";
  f2.open(datafile_WM.str().c_str(), ios::in); 

  if (f2.fail())
    {
      cerr << "Error opening data file " << datafile_WM.str() << endl;
      exit(-1);
    }
  
  //Read central value
  string line;
  for(int i=0; i<12; i++)
    {
      getline(f1,line);
      getline(f2,line);
    }
  
  double ddum;
  char comma;

  for(int i=0; i<fNData/2; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      lstream >> fKin1[i]        >> comma
	      >> ddum            >> comma
	      >> ddum            >> comma
	      >> fData[i];       //in [pb]
	
      fData[i] *= 1000.; //Convert to [fb]
      fKin2[i] = MW*MW;
      fKin3[i] = 8000.;  //GeV

      //Systematic uncertainties
      for(int j=0; j<fNSys-1; j++)
	{
	  lstream >> comma >> fSys[i][j].mult
		  >> comma >> ddum;
	  
	  fSys[i][j].add = fSys[i][j].mult/100. * fData[i];
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = "CORR";
	}

      //Statistical uncertainty
      lstream >> comma >> fStat[i];
      fStat[i] = fStat[i]/100. * fData[i]; //convert to absolute value

      //Luminosity uncertainty
      fSys[i][47].mult = 1.9;
      fSys[i][47].add = fSys[i][47].mult/100. * fData[i];
      fSys[i][47].type = MULT;
      fSys[i][47].name = "ATLASLUMI12";

      //Uncorrelated uncertainties
      fSys[i][6].name  = "UNCORR";
      fSys[i][43].name = "UNCORR";
      fSys[i][44].name = "UNCORR";
      fSys[i][46].name = "UNCORR";
    }
  
  for(int i=fNData/2; i<fNData; i++)
    {
      getline(f2,line);
      istringstream lstream(line);
      lstream >> fKin1[i]        >> comma
	      >> ddum            >> comma
	      >> ddum            >> comma
	      >> fData[i];       //in [pb]
      
      fData[i] *= 1000.; //Convert to [fb]
      fKin2[i] = MW*MW;
      fKin3[i] = 8000.;  //GeV
      
      //Systematic uncertainties
      for(int j=0; j<fNSys-1; j++)
	{
	  lstream >> comma >> fSys[i][j].mult
		  >> comma >> ddum;
	  fSys[i][j].add = fSys[i][j].mult/100. * fData[i];
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = "CORR";
	}
      
      //Statistical uncertainty
      lstream >> comma >> fStat[i];
      fStat[i] = fStat[i]/100. * fData[i]; //convert to absolute value
      
      //Luminosity uncertainty
      fSys[i][47].mult = 1.9;
      fSys[i][47].add = fSys[i][47].mult/100. * fData[i];
      fSys[i][47].type = MULT;
      fSys[i][47].name = "ATLASLUMI12";

      //Uncorrelated uncertainties
      fSys[i][6].name  = "UNCORR";
      fSys[i][43].name = "UNCORR";
      fSys[i][44].name = "UNCORR";
      fSys[i][46].name = "UNCORR";
    }
  
  f1.close();
  f2.close();
  
}

*/
