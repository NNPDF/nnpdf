/*
Name_exp  : ATLAS_1JET_8TEV_R06_DEC
Reference : Measurement of the inclusive jet cross-sections in
            proton–proton collisions at √s = 8 TeV with the ATLAS detector 
ArXiv     : arXiv:1706.03192
Published : JHEP09(2017)020
Hepdata   : https://www.hepdata.net/record/ins1604271

Inclusive jet production cross-sections are measured in proton–proton 
collisions at a centre of-mass energy of √s = 8 TeV recorded by the ATLAS 
experiment at the Large Hadron Collider at CERN. The total integrated 
luminosity of the analysed data set amounts to 20.2 fb−1.
Double-differential cross-sections are measured for jets defined by the anti-kt
jet clustering algorithm with radius parameter of R = 0.6 and are presented as 
a function of the jet transverse momentum, in the range between 70GeV and 
2.5TeV and in six bins of the absolute jet rapidity, between 0 and 3.0.

The decorrelation model in Tab.6 of arXiv:1706.03192 is implemented.

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


#include "ATLAS_1JET_8TEV_R06_DEC.h"
#include <cmath>

//Functions to split uncertainties
double L(double x, double xmin, double xmax)
{
  double LL;

  if(x<xmin)
    LL = 0.;
  else if(x>xmax)
    LL = 1.;
  else
    LL = (x-xmin)/(xmax-xmin);
  
  return LL;
}

void splitErrors(double pT, double y, double uncertainty, int opt, double uncertainties[3])
{
  if(opt==14)
    {
      double L00 = L(log(pT),log(0.1),log(2.5));
      double L01 = sqrt(1.-pow(L(y,0.,1.5),2.));
      uncertainties[0] = L00*L01*uncertainty;

      double L10 = L(log(pT),log(0.1),log(2.5));
      double L11 = L(y,1.,3.);
      uncertainties[1] = L10*L11*uncertainty;
      
      uncertainties[2] = sqrt(pow(uncertainty,2.)-pow(uncertainties[0],2.)-pow(uncertainties[1],2.));
    }
  
  else if(opt==16)
    {
      double L00 = sqrt(1.-pow(L(log(pT),log(0.1),log(2.5)),2.));
      double L01 = sqrt(1.-pow(L(y,0.,1.5),2));
      uncertainties[0] =  L00*L01*uncertainty;

      double L10 = sqrt(1.-pow(L(log(pT),log(0.1),log(2.5)),2.));
      double L11 = L(y,1.5,3.);
      uncertainties[1] = L10*L11*uncertainty;
      
      uncertainties[2] = sqrt(pow(uncertainty,2.)-pow(uncertainties[0],2.)-pow(uncertainties[1],2.));
    }
  
  else if(opt==17)
    {
      double L00 = sqrt(1.-pow(L(log(pT),log(0.1),log(2.5)),2.));
      double L01 = sqrt(1.-pow(L(y,0.,1.5),2));
      uncertainties[0] = L00*L01*uncertainty;

      double L10 = sqrt(1.-pow(L(log(pT),log(0.1),log(2.5)),2.));
      double L11 = L(y,1.,3.);      
      uncertainties[1] = L10*L11*uncertainty;

      uncertainties[2] = sqrt(pow(uncertainty,2.)-pow(uncertainties[0],2.)-pow(uncertainties[1],2.));
    }
  else
    {
      cout << "Splitting option not available " << opt << endl;
    }
}

void ATLAS_1JET_8TEV_R06_DECFilter::ReadData()
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
    DataFile << dataPath() << "rawdata/ATLAS_1JET_8TEV_R06/" << data_file;
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
      //note: sys1 and sys2 can be: both positive, both negative,
      //negative and positive, positive and negative
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

      /*
      //Decorrelate uncertainties
      fSys[i][528].name = "UNCORR";
      fSys[i][529].name = "UNCORR";
      fSys[i][580].name = "UNCORR";
      fSys[i][581].name = "UNCORR";
      fSys[i][610].name = "UNCORR";
      fSys[i][611].name = "UNCORR";
      */

      //Split uncertainties 659++ (18 additional uncertainties, most of which are zero
      double uncertainties[3];

      for(int k=0; k<fNSys; k++)
	{
	  
	  if(k==528) //JES flavour separation
	    {
	      splitErrors(sqrt(fKin2[i]), fKin1[i], fSys[i][k].add, 14, uncertainties);

	      fSys[i][528].add  = 0.;
	      fSys[i][528].mult = 0.;
	      
	      fSys[i][659].type = MULT;   
	      fSys[i][659].name = "CORR";   
	      fSys[i][659].add  = uncertainties[0];
	      fSys[i][659].mult = fSys[i][659].add*100/fData[i];

	      fSys[i][660].type = MULT;   
	      fSys[i][660].name = "CORR";   
	      fSys[i][660].add  = uncertainties[1];
	      fSys[i][660].mult = fSys[i][660].add*100/fData[i];

	      fSys[i][661].type = MULT;   
	      fSys[i][661].name = "UNCORR";   
	      fSys[i][661].add  = uncertainties[2];
	      fSys[i][661].mult = fSys[i][661].add*100/fData[i];
	    }

	  if(k==529) //JES flavour separation
	    {
	      splitErrors(sqrt(fKin2[i]), fKin1[i], fSys[i][k].add, 14, uncertainties);

	      fSys[i][529].add  = 0.;
	      fSys[i][529].mult = 0.;
	      
	      fSys[i][662].type = MULT;   
	      fSys[i][662].name = "CORR";   
	      fSys[i][662].add  = uncertainties[0];
	      fSys[i][662].mult = fSys[i][662].add*100/fData[i];

	      fSys[i][663].type = MULT;   
	      fSys[i][663].name = "CORR";   
	      fSys[i][663].add  = uncertainties[1];
	      fSys[i][663].mult = fSys[i][663].add*100/fData[i];

	      fSys[i][664].type = MULT;   
	      fSys[i][664].name = "UNCORR";   
	      fSys[i][664].add  = uncertainties[2];
	      fSys[i][664].mult = fSys[i][664].add*100/fData[i];
	      
	    }

	  if(k==580) //JES MJB fragmentation
	    {
	      splitErrors(sqrt(fKin2[i]), fKin1[i], fSys[i][k].add, 17, uncertainties);

	      fSys[i][580].add  = 0.;
	      fSys[i][580].mult = 0.;
	      
	      fSys[i][665].type = MULT;   
	      fSys[i][665].name = "CORR";   
	      fSys[i][665].add  = uncertainties[0];
	      fSys[i][665].mult = fSys[i][665].add*100/fData[i];

	      fSys[i][666].type = MULT;   
	      fSys[i][666].name = "CORR";   
	      fSys[i][666].add  = uncertainties[1];
	      fSys[i][666].mult = fSys[i][666].add*100/fData[i];

	      fSys[i][667].type = MULT;   
	      fSys[i][667].name = "UNCORR";   
	      fSys[i][667].add  = uncertainties[2];
	      fSys[i][667].mult = fSys[i][667].add*100/fData[i];	      
	    }

	  if(k==581) //JES MJB fragmentation
	    {
	      splitErrors(sqrt(fKin2[i]), fKin1[i], fSys[i][k].add, 17, uncertainties);

	      fSys[i][581].add  = 0.;
	      fSys[i][581].mult = 0.;
	      
	      fSys[i][668].type = MULT;   
	      fSys[i][668].name = "CORR";   
	      fSys[i][668].add  = uncertainties[0];
	      fSys[i][668].mult = fSys[i][668].add*100/fData[i];

	      fSys[i][669].type = MULT;   
	      fSys[i][669].name = "CORR";   
	      fSys[i][669].add  = uncertainties[1];
	      fSys[i][669].mult = fSys[i][669].add*100/fData[i];

	      fSys[i][670].type = MULT;   
	      fSys[i][670].name = "UNCORR";   
	      fSys[i][670].add  = uncertainties[2];
	      fSys[i][670].mult = fSys[i][670].add*100/fData[i];
	    }

	  if(k==610) //JES pile-up
	    {
	      splitErrors(sqrt(fKin2[i]), fKin1[i], fSys[i][k].add, 16, uncertainties);

	      fSys[i][610].add  = 0.;
	      fSys[i][610].mult = 0.;
	      
	      fSys[i][671].type = MULT;   
	      fSys[i][671].name = "CORR";   
	      fSys[i][671].add  = uncertainties[0];
	      fSys[i][671].mult = fSys[i][671].add*100/fData[i];

	      fSys[i][672].type = MULT;   
	      fSys[i][672].name = "CORR";   
	      fSys[i][672].add  = uncertainties[1];
	      fSys[i][672].mult = fSys[i][672].add*100/fData[i];

	      fSys[i][673].type = MULT;   
	      fSys[i][673].name = "UNCORR";   
	      fSys[i][673].add  = uncertainties[2];
	      fSys[i][673].mult = fSys[i][673].add*100/fData[i];
	    }

	  if(k==611) //JES pile-up
	    {
	      splitErrors(sqrt(fKin2[i]), fKin1[i], fSys[i][k].add, 16, uncertainties);

	      fSys[i][611].add  = 0.;
	      fSys[i][611].mult = 0.;
	      
	      fSys[i][674].type = MULT;   
	      fSys[i][674].name = "CORR";   
	      fSys[i][674].add  = uncertainties[0];
	      fSys[i][674].mult = fSys[i][674].add*100/fData[i];

	      fSys[i][675].type = MULT;   
	      fSys[i][675].name = "CORR";   
	      fSys[i][675].add  = uncertainties[1];
	      fSys[i][675].mult = fSys[i][675].add*100/fData[i];

	      fSys[i][676].type = MULT;   
	      fSys[i][676].name = "UNCORR";   
	      fSys[i][676].add  = uncertainties[2];
	      fSys[i][676].mult = fSys[i][676].add*100/fData[i];
	    } 
	  
	}
    }
    
    rS.close();
    n += ndata[bin];

  }
  
}

