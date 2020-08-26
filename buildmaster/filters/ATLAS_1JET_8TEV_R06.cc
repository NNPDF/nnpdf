/*
Name_exp  : ATLAS_1JET_8TEV_R06
Reference : Measurement of the inclusive jet cross-sections in
            proton–proton collisions at √s = 8 TeV with the ATLAS detector 
ArXiv     : arXiv:1706.03192
Published : JHEP09(2017)020
Hepdata   : https://www.hepdata.net/record/ins1604271

Inclusive jet production cross-sections are measured in proton–proton collisions at a centre of-mass
energy of √s = 8 TeV recorded by the ATLAS experiment at the Large Hadron Collider
at CERN. The total integrated luminosity of the analysed data set amounts to 20.2 fb−1.
Double-differential cross-sections are measured for jets defined by the anti-kt
jet clustering algorithm with radius parameter of R = 0.6 and are presented as a function
of the jet transverse momentum, in the range between 70GeV and 2.5TeV and in six bins of
the absolute jet rapidity, between 0 and 3.0.



bin 1:   0 < |y| < 0.5
========================
points    34 
real sys  329x2 + lumi 


bin 2:   0.5 < |y| < 1.0
========================
points    33 
real sys  329x2 + lumi


bin 3:   1.0 < |y| < 1.5
========================
points    32 
real sys  329x2 + lumi 


bin 4:   1.5 < |y| < 2.0
========================
points    30 
real sys  329x2 + lumi 


bin 5:   2.0 < |y| < 2.5
========================
points    24 
real sys  329x2 + lumi 


bin 6:   2.5 < |y| < 3.0
========================
points    18 
real sys  329x2 + lumi 


tot points         171
tot sys per point  659

Implemented by TG July 2019. 
*/


#include "ATLAS_1JET_8TEV_R06.h"

void ATLAS_1JET_8TEV_R06Filter::ReadData()
{

  //opening files
  fstream rS, rCorr;

  //bins specification
  int nbins = 6;
  std::vector<double> y = {0.25, 0.75, 1.25, 1.75, 2.25, 2.75};
  std::vector<int> ndata  = {34, 33, 32, 30, 24, 18};    
  std::vector<double> stat(fNData);   

  int n = 0;                                //count total number of datapoints

  for(int bin=0; bin < nbins; bin++ )
  {

    string data_file = "/bin_" + to_string(bin+1) + ".dat";
    stringstream DataFile("");
    DataFile << dataPath() << "rawdata/" << fSetName << data_file;
    rS.open(DataFile.str().c_str(), ios::in);
    if (rS.fail()) {
      cerr << "Error opening data file " << DataFile.str() << endl;
      exit(-1);
    }
  
    string line;
    double pt, dum;    

    for (int i = n; i < n + ndata[bin]; i++)
    {

      getline(rS,line);                     
      istringstream lstream(line);          

      lstream >> pt;
      lstream >> dum >> dum;
      
      fKin1[i] = y[bin];           // y, central value of the bin
      fKin2[i] = pow(pt,2.);       // pt2, central value of the bin
      fKin3[i] = 8000;             // sqrt(s)

      lstream >> fData[i];         // cross section [fb/GeV]

      //Read uncertainties
      double sys1, sys2;

      //1) Stat
      lstream >> sys1 >> sys2;
      fStat[i] = sys1;               

      //2) Systematics
      //note: sys1 and sys2 can be: both positive, both negative, negative and positive, positive and negative
      for (int k=0; k<329; k++)
	{

	  lstream >> sys1 >> sys2;
	  	  
	  //sort out uncertainties
	  double tmp1, tmp2;
	  tmp1 = sys1;
	  tmp2 = sys2;

	  //case 1: sys1 and sys2 are both negative
	  if(tmp1<0.0 && tmp2<0.0)
	    {
	      if(tmp2<tmp1)
		{
		  sys1 = 0.0;
		  sys2 = tmp2;
		}
	      if(tmp2>tmp1)
		{
		  sys1 = 0.0;
		  sys2 = tmp1;
		}
	    }
	  
	  //Case 2: sys1 and sys2 are both positive
	  if(tmp1>0.0 && tmp2>0.0)
	    {
	      if(tmp1>tmp2)
		{
		  sys1 = tmp1;
		  sys2 = 0.0;
		}
	      if(tmp1<tmp2)
		{
		  sys1 = tmp2;
		  sys2 = 0.0;
		}
	    }	  

	  fSys[i][2*k].type = MULT;   
	  fSys[i][2*k].name = "CORR";
	  fSys[i][2*k].add = sys1/sqrt(2.);
	  fSys[i][2*k].mult  = fSys[i][2*k].add*100/fData[i];

	  fSys[i][2*k+1].type = MULT;   
	  fSys[i][2*k+1].name = "CORR";
	  fSys[i][2*k+1].add = sys2/sqrt(2.);
	  fSys[i][2*k+1].mult  = fSys[i][2*k+1].add*100/fData[i];
	}

      //5) Luminosity uncertainty
      lstream >> sys1 >> sys2;
      fSys[i][658].type = MULT;   
      fSys[i][658].name = "ATLASLUMI12";   
      fSys[i][658].add = sys1;
      fSys[i][658].mult  = fSys[i][658].add*100/fData[i];

      //Decorrelate uncertainties
      fSys[i][528].name = "UNCORR";
      fSys[i][529].name = "UNCORR";
      fSys[i][580].name = "UNCORR";
      fSys[i][581].name = "UNCORR";
      fSys[i][610].name = "UNCORR";
      fSys[i][611].name = "UNCORR";
      
      /*
      cout << fSys[i][528].add*sqrt(2.) << "   " << fSys[i][529].add*sqrt(2.) << "   "
	   << fSys[i][580].add*sqrt(2.) << "   " << fSys[i][581].add*sqrt(2.) << "   "
	   << fSys[i][610].add*sqrt(2.) << "   " << fSys[i][611].add*sqrt(2.) << endl;
      */
      
    }
    
    rS.close();
    n += ndata[bin];

  }
  
}

