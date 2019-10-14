/*
write here a description of the data set
*/

#include "ATLAS_WJET_8TEV.h"

//Distribution differential in PT

void ATLAS_WP_JET_8TEV_PTFilter::ReadData()
{
  fstream f1;
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/ATLAS_WJET_8TEV/HEPData-ins1635273-v1-Table_13.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read results
  string line;
  for(int i=0; i<16; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      double ddum, sysL[fNSys/2], sysR[fNSys/2], lumi;
      char comma;
      lstream >> fKin1[i] >> comma
	      >> ddum     >> comma
	      >> ddum     >> comma
	      >> fData[i] >> comma
	      >> fStat[i] >> comma
	      >> ddum     >> comma;

      fKin2[i] = 0;
      fKin3[i] = 8000; //GeV

      for (int j=0; j<(fNSys+1)/2-4; j++)
	{
	  lstream >> sysL[j] >> comma >> sysR[j] >> comma;

	  cout << sysL[j] << "   " << sysR[j] << endl;

	  sysL[j] = sysL[j]/sqrt(2);
	  sysR[j] = sysR[j]/sqrt(2);

	  if(sysL[j]<0 && sysR[j]<0)
	    {
	      cout << "same sign" << endl;
	    }
	  if(sysL[j]>0 && sysR[j]>0)
	    {
	      cout << "same sign" << endl;
	    }

	  fSys[i][2*j].add = sysL[j];
	  fSys[i][2*j].mult = fSys[i][2*j].add*1e2/fData[i];
	  fSys[i][2*j].type = MULT;
	  fSys[i][2*j].name = "ATLASWJ";
	  
	  fSys[i][2*j+1].add = sysR[j];
	  fSys[i][2*j+1].mult = fSys[i][2*j+1].add*1e2/fData[i];
	  fSys[i][2*j+1].type = MULT;
	  fSys[i][2*j+1].name = "ATLASWJ";

	}

      lstream >> lumi     >> comma
	      >> ddum     >> comma
	      >> sysL[51] >> comma 
	      >> sysR[51] >> comma
	      >> sysL[52] >> comma 
	      >> sysR[52] >> comma
	      >> sysL[53] >> comma 
	      >> sysR[53] >> comma;
      
      fSys[i][100].add = lumi;
      fSys[i][100].mult = fSys[i][100].add*1e2/fData[i];
      fSys[i][100].type = MULT;
      fSys[i][100].name = "ATLASLUMI12"; 

      fSys[i][101].add = sysL[51]/sqrt(2);
      fSys[i][102].add = sysR[51]/sqrt(2);
      fSys[i][103].add = sysL[52]/sqrt(2);
      fSys[i][104].add = sysR[52]/sqrt(2);
      fSys[i][105].add = sysL[53]/sqrt(2);
      fSys[i][106].add = sysR[53]/sqrt(2);
      
      for(int j=101; j<fNSys; j++)
	{
	  fSys[i][j].mult = fSys[i][j].add*1e2/fData[i];
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = "UNCORR"; 
	}

    }

  f1.close();

} 

