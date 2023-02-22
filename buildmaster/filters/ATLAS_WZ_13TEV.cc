/*
Reference:
   arXiv:     [1902.05759]
   hepdata:   https://www.hepdata.net/record/ins1720438
   published: Eur.Phys.J. C79 (2019) 535, 2019 
Description:
   Measurements of W± Z production cross sections in pp collisions at a centre-
   of-mass energy of 13 TeV. The data were collected in 2015 and 2016 by the 
   ATLAS experiment at the Large Hadron Collider, and correspond to an 
   integrated luminosity of 36.1 fb-1 . The W± Z candidate events are 
   reconstructed using leptonic decay modes of the gauge bosons into electrons 
   and muons. Differential distributions in the following kinematic variables
   are implemented: pTZ, pTW, mT(W,Z), and phi(W,Z). The total fiducial cross 
   section is also implemented. The implementation is based on 
   Tabs. 8, 10, 12, 14 and 1 available on HepData.
*/

#include "ATLAS_WZ_13TEV.h"

//ATLAS_WZ_13TEV_pTZ: combined pT spectrum
void ATLAS_WZ_13TEV_pTZFilter::ReadData()
{
  fstream f1;

  //Central values and breakdown of uncertainties
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/ATLAS_WZ_13TEV/Table8.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Read central values and uncertainties
  string line;
  for(int i=0; i<17; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<fNData; i++)
    {
      double ddum;
      double stat;
      char comma;
      getline(f1,line);
      istringstream lstream(line);

      fKin2[i] = 0.;
      fKin3[i] = 13000.;    //c.m. energy
      lstream >> fKin1[i]        >> comma
	      >> ddum            >> comma
	      >> ddum            >> comma
	      >> fData[i]        >> comma
	      >> stat            >> comma;

      fStat[i] = stat/100.*fData[i];

      for(int j=0; j<fNSys; j++)
	{
	  lstream >> comma >> fSys[i][j].mult >> comma
		  >> comma >> ddum            >> comma;
	  fSys[i][j].add = fSys[i][j].mult*fData[i]/100.;
	}
      
      //Unocrrelated systematic uncertainty
      fSys[i][0].type = ADD;
      fSys[i][0].name = "UNCORR";
      
      //Unfolding uncertainty
      fSys[i][1].type = ADD;
      fSys[i][1].name = "CORR";

      //Electron
      fSys[i][2].type = MULT;
      fSys[i][2].name = "CORR";

      //Muon
      fSys[i][3].type = MULT;
      fSys[i][3].name = "CORR";

      //Jets
      fSys[i][4].type = MULT;
      fSys[i][4].name = "CORR";

      //Background1
      fSys[i][5].type = MULT;
      fSys[i][5].name = "CORR";

      //Background2
      fSys[i][6].type = MULT;
      fSys[i][6].name = "CORR";

      //Pileup
      fSys[i][7].type = MULT;
      fSys[i][7].name = "CORR";

      //Luminosity
      fSys[i][8].type = MULT;
      fSys[i][8].name = "ATLASLUMI13";

    }

  f1.close();

}

//ATLAS_WZ_13TEV_pTW: combined pT spectrum
void ATLAS_WZ_13TEV_pTWFilter::ReadData()
{
  fstream f1;

  //Central values and breakdown of uncertainties
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/ATLAS_WZ_13TEV/Table10.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Read central values and uncertainties
  string line;
  for(int i=0; i<17; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<fNData; i++)
    {
      double ddum;
      double stat;
      char comma;
      getline(f1,line);
      istringstream lstream(line);

      fKin2[i] = 0.;
      fKin3[i] = 13000.;    //c.m. energy
      lstream >> fKin1[i]        >> comma
	      >> ddum            >> comma
	      >> ddum            >> comma
	      >> fData[i]        >> comma
	      >> stat            >> comma;

      fStat[i] = stat/100.*fData[i];

      for(int j=0; j<fNSys; j++)
	{
	  lstream >> comma >> fSys[i][j].mult >> comma
		  >> comma >> ddum            >> comma;
	  fSys[i][j].add = fSys[i][j].mult*fData[i]/100.;
	}
      
      //Unocrrelated systematic uncertainty
      fSys[i][0].type = ADD;
      fSys[i][0].name = "UNCORR";
      
      //Unfolding uncertainty
      fSys[i][1].type = ADD;
      fSys[i][1].name = "CORR";

      //Electron
      fSys[i][2].type = MULT;
      fSys[i][2].name = "CORR";

      //Muon
      fSys[i][3].type = MULT;
      fSys[i][3].name = "CORR";

      //Jets
      fSys[i][4].type = MULT;
      fSys[i][4].name = "CORR";

      //Background1
      fSys[i][5].type = MULT;
      fSys[i][5].name = "CORR";

      //Background2
      fSys[i][6].type = MULT;
      fSys[i][6].name = "CORR";

      //Pileup
      fSys[i][7].type = MULT;
      fSys[i][7].name = "CORR";

      //Luminosity
      fSys[i][8].type = MULT;
      fSys[i][8].name = "ATLASLUMI13";

    }

  f1.close();

}

//ATLAS_WZ_13TEV_mTWZ: combined pT spectrum
void ATLAS_WZ_13TEV_mTWZFilter::ReadData()
{
  fstream f1;

  //Central values and breakdown of uncertainties
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/ATLAS_WZ_13TEV/Table12.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Read central values and uncertainties
  string line;
  for(int i=0; i<17; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<fNData; i++)
    {
      double ddum;
      double stat;
      char comma;
      getline(f1,line);
      istringstream lstream(line);

      fKin2[i] = 0.;
      fKin3[i] = 13000.;    //c.m. energy
      lstream >> fKin1[i]        >> comma
	      >> ddum            >> comma
	      >> ddum            >> comma
	      >> fData[i]        >> comma
	      >> stat            >> comma;

      fStat[i] = stat/100.*fData[i];

      for(int j=0; j<fNSys; j++)
	{
	  lstream >> comma >> fSys[i][j].mult >> comma
		  >> comma >> ddum            >> comma;
	  fSys[i][j].add = fSys[i][j].mult*fData[i]/100.;
	}
      
      //Unocrrelated systematic uncertainty
      fSys[i][0].type = ADD;
      fSys[i][0].name = "UNCORR";
      
      //Unfolding uncertainty
      fSys[i][1].type = ADD;
      fSys[i][1].name = "CORR";

      //Electron
      fSys[i][2].type = MULT;
      fSys[i][2].name = "CORR";

      //Muon
      fSys[i][3].type = MULT;
      fSys[i][3].name = "CORR";

      //Jets
      fSys[i][4].type = MULT;
      fSys[i][4].name = "CORR";

      //Background1
      fSys[i][5].type = MULT;
      fSys[i][5].name = "CORR";

      //Background2
      fSys[i][6].type = MULT;
      fSys[i][6].name = "CORR";

      //Pileup
      fSys[i][7].type = MULT;
      fSys[i][7].name = "CORR";

      //Luminosity
      fSys[i][8].type = MULT;
      fSys[i][8].name = "ATLASLUMI13";

    }

  f1.close();

}

//ATLAS_WZ_13TEV_phiWZ: combined pT spectrum
void ATLAS_WZ_13TEV_phiWZFilter::ReadData()
{
  fstream f1;

  //Central values and breakdown of uncertainties
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/ATLAS_WZ_13TEV/Table14.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Read central values and uncertainties
  string line;
  for(int i=0; i<17; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<fNData; i++)
    {
      double ddum;
      double stat;
      char comma;
      getline(f1,line);
      istringstream lstream(line);

      fKin2[i] = 0.;
      fKin3[i] = 13000.;    //c.m. energy
      lstream >> fKin1[i]        >> comma
	      >> ddum            >> comma
	      >> ddum            >> comma
	      >> fData[i]        >> comma
	      >> stat            >> comma;

      fStat[i] = stat/100.*fData[i];

      for(int j=0; j<fNSys; j++)
	{
	  lstream >> comma >> fSys[i][j].mult >> comma
		  >> comma >> ddum            >> comma;
	  fSys[i][j].add = fSys[i][j].mult*fData[i]/100.;
	}
      
      //Unocrrelated systematic uncertainty
      fSys[i][0].type = ADD;
      fSys[i][0].name = "UNCORR";
      
      //Unfolding uncertainty
      fSys[i][1].type = ADD;
      fSys[i][1].name = "CORR";

      //Electron
      fSys[i][2].type = MULT;
      fSys[i][2].name = "CORR";

      //Muon
      fSys[i][3].type = MULT;
      fSys[i][3].name = "CORR";

      //Jets
      fSys[i][4].type = MULT;
      fSys[i][4].name = "CORR";

      //Background1
      fSys[i][5].type = MULT;
      fSys[i][5].name = "CORR";

      //Background2
      fSys[i][6].type = MULT;
      fSys[i][6].name = "CORR";

      //Pileup
      fSys[i][7].type = MULT;
      fSys[i][7].name = "CORR";

      //Luminosity
      fSys[i][8].type = MULT;
      fSys[i][8].name = "ATLASLUMI13";

    }

  f1.close();

}

//ATLAS_WZ_13TEV_totWZ: combined pT spectrum
void ATLAS_WZ_13TEV_totWZFilter::ReadData()
{
  fstream f1;

  //Central values and breakdown of uncertainties
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/ATLAS_WZ_13TEV/Table1.csv";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Read central values and uncertainties
  string line;
  for(int i=0; i<38; i++)
    {
      getline(f1,line);
    }

  for(int i=0; i<fNData; i++)
    {
      double ddum;
      char comma;
      getline(f1,line);
      istringstream lstream(line);

      fKin1[i] = 0.;
      fKin2[i] = 0.;
      fKin3[i] = 13000.;    //c.m. energy
      lstream >> ddum           >> comma
	      >> fData[i]       >> comma
	      >> fStat[i]       >> comma
	      >> ddum           >> comma
	      >> fSys[i][0].add >> comma
	      >> ddum           >> comma
	      >> fSys[i][1].add >> comma
	      >> ddum           >> comma
	      >> fSys[i][2].add >> comma
	      >> ddum;

      for(int j=0; j<fNSys; j++)
	{
	  fSys[i][j].mult = fSys[i][j].add/fData[i]*100.;
	}
      
      //Total systematic uncertainty
      fSys[i][0].type = MULT;
      fSys[i][0].name = "CORR";
      
      //Modelling uncertainty
      fSys[i][1].type = MULT;
      fSys[i][1].name = "CORR";

      //Luminosity uncertainty
      fSys[i][2].type = MULT;
      fSys[i][2].name = "ATLASLUMI13";

    }

  f1.close();

}
