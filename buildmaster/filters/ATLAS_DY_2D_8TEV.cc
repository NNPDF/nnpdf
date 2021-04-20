//Some description is in order here


#include "ATLAS_DY_2D_8TEV.h"

void ATLAS_DY_2D_8TEVFilter::ReadData()
{
  fstream f1;

  //Full breakdown of systematic ucnertainties
  stringstream datafile("");
  datafile << dataPath()
	  << "rawdata/ATLAS_DY_2D_8TEV/HEPData-ins1630886-v3-Table_5.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read central value and uncertainties
  string line;
  for(int i=0; i<15; i++)
    {
      getline(f1,line);     
    }

  double ddum;
  char comma;

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      lstream >> ddum >> comma
	      >> ddum >> comma
	      >> ddum >> comma
	      >> fKin1[i] >> comma
	      >> ddum >> comma
	      >> ddum >> comma
	      >> fKin2[i] >> comma
	      >> ddum >> comma
	      >> ddum >> comma
	      >> fData[i] >> comma  //[pb/GeV]
	      >> fStat[i] >> comma
	      >> ddum;

      fKin2[i] = fKin2[i]*fKin2[i]; //[GeV2]
      fKin3[i] = 8000.;             //[GeV]

      //Convert to [fb]
      fData[i] *= 1000.;
      fStat[i] *= 1000.;

      //Systematic ucnertainties
      for(int j=0; j<fNSys-1; j++)
	{
	  lstream >> comma >> fSys[i][j].add
		  >> comma >> ddum;

	  fSys[i][j].add *= 1000.; //Convert to [fb]
	  
	  fSys[i][j].mult = fSys[i][j].add/fData[i]*100.;
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = "CORR";
	}

      //Luminosity uncertainty
      fSys[i][276].mult = 1.8;
      fSys[i][276].add = fSys[i][276].mult/100. * fData[i];
      fSys[i][276].type = MULT;
      fSys[i][276].name = "ATLASLUMI12";

      //Uncorrelated uncertainties
      fSys[i][275].name = "UNCORR";      
    }

  f1.close();
}



