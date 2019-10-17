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
	  lstream >> sysR[j] >> comma >> sysL[j] >> comma;

	  sysR[j] /= sqrt(2.);
	  sysL[j] /= sqrt(2.);

	  double tmp1, tmp2;
	  tmp1=sysR[j];
	  tmp2=sysL[j];
	  /*	  
	  //Case 1: sysL and sysR are both negative
	  if(tmp1<0 && tmp2<0)
	    {
	      if(tmp2<tmp1 || tmp2==tmp1)
		{
		  sysR[j] = 0.0;
		  sysL[j] = tmp2;
		}
	      if(tmp2>tmp1)
		{
		  sysR[j] = 0.0;
		  sysL[j] = tmp1;
		}
	    }

	  //Case 2: sysL and sysR are both positive
	  if(tmp1>0.0 && tmp2>0.0)
	    {
	      if(tmp1>tmp2 || tmp1==tmp2)
		{
		  sysR[j] = tmp1;
		  sysL[j] = 0.0;
		}
	      if(tmp1<tmp2)
		{
		  sysR[j] = tmp2;
		  sysL[j] = 0.0;
		}
	    }
	  
	  //Case3: sys1 is negative and sys2 is positive
	  if(tmp1<0.0 && tmp2>0.0)
	    {
	      sys1[idat][isys] = tmp2;
	      sys2[idat][isys] = tmp1;
	    }
	  */	  

	  fSys[i][2*j].add = sysR[j];
	  fSys[i][2*j].mult = fSys[i][2*j].add*1e2/fData[i];
	  fSys[i][2*j].type = MULT;
	  fSys[i][2*j].name = "ATLASWJ";
	  //fSys[i][2*j].name = "CORR";
	  
	  fSys[i][2*j+1].add = sysL[j];
	  fSys[i][2*j+1].mult = fSys[i][2*j+1].add*1e2/fData[i];
	  fSys[i][2*j+1].type = MULT;
	  fSys[i][2*j+1].name = "ATLASWJ";
	  //fSys[i][2*j+1].name = "CORR";

	}

      lstream >> lumi     >> comma
	      >> ddum     >> comma
	      >> sysR[51] >> comma 
	      >> sysL[51] >> comma
	      >> sysR[52] >> comma 
	      >> sysL[52] >> comma
	      >> sysR[53] >> comma 
	      >> sysL[53] >> comma;
      
      fSys[i][100].add = lumi;
      fSys[i][100].mult = fSys[i][100].add*1e2/fData[i];
      fSys[i][100].type = MULT;
      fSys[i][100].name = "ATLASLUMI12"; 

      fSys[i][101].add = sysR[51]/sqrt(2);
      fSys[i][102].add = sysL[51]/sqrt(2);
      fSys[i][103].add = sysR[52]/sqrt(2);
      fSys[i][104].add = sysL[52]/sqrt(2);
      fSys[i][105].add = sysR[53]/sqrt(2);
      fSys[i][106].add = sysL[53]/sqrt(2);
      
      for(int j=101; j<fNSys; j++)
	{
	  fSys[i][j].mult = fSys[i][j].add*1e2/fData[i];
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = "UNCORR"; 
	}

      fSys[i][96].name  = "UNCORR";
      fSys[i][97].name  = "UNCORR";

    }

  f1.close();

} 

