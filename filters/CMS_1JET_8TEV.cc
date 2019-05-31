/*

Points implemented are those reported at 
https://www.hepdata.net/record/ins1487277 
tables 1 to 6.

Full breakdown of the systematics from xfitter files. 
It contains 26 syst to be simmetrized, given as 52 column entries, 
plus the lumi sys and 3 columns related to stats.

bin 1:   0 < |y| < 0.5
========================
points    46 
real sys  26 + lumi 
art sys   37
systot    64

bin 2:   0.5 < |y| < 1.0
========================
points    46 
real sys  26 + lumi 
art sys   37
systot    64

bin 3:   1.5 < |y| < 2.0
========================
points    45 
real sys  26 + lumi 
art sys   36
systot    63

bin 4:   2.0 < |y| < 2.5
========================
points    41 
real sys  26 + lumi 
art sys   32
systot    59

bin 5:   2.5 < |y| < 3.0
========================
points    34 
real sys  26 + lumi 
art sys   25
systot    52

bin 6:   3.0 < |y| < 3.5
========================
points    27 
real sys  26 + lumi 
art sys   18
systot    45


tot points         239
tot sys per point   64  (those missing will be set to ~0)


*********************************************/


#include "CMS.h"

void CMS_1JET_8TEVFilter::ReadData()
{

  //opening files
  fstream rS, rCorr;

  //bins specification
  int nbins = 6;
  double  y[]  = {0.25, 0.75, 1.25, 1.75, 2.25, 2.75 };
  int ndata[]  = {46, 46, 45, 41, 34, 27};    
  int artsys[] = {37, 37, 36, 32, 25, 18};  //number of artificial systematics for each bin
  double stat[fNData];   

  int n = 0;              //count total number of datapoints
  const double fac = 1e7; //conversion factor from pb to 10 mub


  for(int bin=0; bin < nbins; bin++ )
  {

    string data_file = "/CMS_8TeV_jets_Ybin" + to_string(bin+1) + ".dat";
    string cov_file =  "/CMS_8TeV_jets_Ybin" + to_string(bin+1) + "___CMS_8TeV_jets_Ybin" + to_string(bin+1) + ".dat";

    //data files
    stringstream DataFile("");
    DataFile << dataPath() << "rawdata/" << fSetName << data_file;
    rS.open(DataFile.str().c_str(), ios::in);
    if (rS.fail()) {
      cerr << "Error opening data file " << DataFile.str() << endl;
      exit(-1);
    }

    //statistical correlation matrix
    stringstream DataFileCorr("");
    DataFileCorr << dataPath() << "rawdata/" << fSetName << cov_file;
    rCorr.open(DataFileCorr.str().c_str(), ios::in);

    if (rCorr.fail()) {
      cerr << "Error opening data file " << DataFileCorr.str() << endl;
      exit(-1);
    }
  
    string line;
    double ptmin, ptmax, dum;

    //skip the first lines
    getline(rS,line);
    for (int i = 0; i < 41; i++)  
      getline(rS,line);

    for (int i = n; i < n + ndata[bin]; i++)
    {
      getline(rS,line);                     
      istringstream lstream(line);          

      lstream >> dum >> dum >> dum;
      lstream >> ptmin >> ptmax;
      
      fKin1[i] = y[bin];                       // y, central value of the bin
      fKin2[i] = pow((ptmin + ptmax)*0.5,2.);  // pt2, central value of the bin
      fKin3[i] = 8000;                         // sqrt(s)

      lstream >> fData[i];                     // cross section [fb/GeV]

      //Read uncertainties
      double sys1, sys2;

      //1) Relative theoretical uncertainties due to nonperturbative corrections in the prediction (they must be ignored)
      lstream >> dum >> dum >> dum;

      //2) Absolute statistical uncertainty
      lstream >> fStat[i];

      //3) Relative systematic uncertainties due to the unfolding procedure
      //note: sys1 always positive; sys2 always negative
      lstream >> sys1 >> sys2;

      fSys[i][0].type = MULT;   
      fSys[i][0].name = "CORR";
      fSys[i][0].mult = sys1/sqrt(2.);
      fSys[i][0].add  = fSys[i][0].mult*fData[i]/100;

      fSys[i][1].type = MULT;   
      fSys[i][1].name = "CORR";
      fSys[i][1].mult = sys2/sqrt(2.);
      fSys[i][1].add  = fSys[i][1].mult*fData[i]/100;

      //4) Relative sytematic JES uncertainties, obtained from 24 independent sources of uncertainty
      //note: sys1 and sys2 can be: both positive, both negative, negative an dpositive, positive and negative
      for (int k=0; k<24; k++)
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

	  fSys[i][2+2*k].type = MULT;   
	  fSys[i][2+2*k].name = "CORR";
	  fSys[i][2+2*k].mult = sys1/sqrt(2.);
	  fSys[i][2+2*k].add  = fSys[i][2+2*k].mult*fData[i]/100;

	  fSys[i][2+2*k+1].type = MULT;   
	  fSys[i][2+2*k+1].name = "CORR";
	  fSys[i][2+2*k+1].mult = sys2/sqrt(2.);
	  fSys[i][2+2*k+1].add  = fSys[i][2+2*k+1].mult*fData[i]/100;
	}

      //5) Luminosity uncertainty
      double sys;
      lstream >> sys;
      fSys[i][50].type = MULT;   
      fSys[i][50].name = "CMSLUMI12";
      fSys[i][50].mult = sys;
      fSys[i][50].add  = fSys[i][50].mult*fData[i]/100;

      //6) Relative uncorrelated systematic uncertainty
      lstream >> dum >> sys;
      fSys[i][51].type = MULT;   
      fSys[i][51].name = "UNCORR";
      fSys[i][51].mult = sys;
      fSys[i][51].add  = fSys[i][51].mult*fData[i]/100;

      /*
      //Checking systematics
      double systotl = 0.;
      double systotr = 0.;

      for (int isys=0; isys<25; isys++)
	{
	  systotl += 2.*pow(fSys[i][2*isys].mult,2);
	  systotr += 2.*pow(fSys[i][2*isys+1].mult,2);
	}

      cout << i << "   " << fStat[i] << "  " 
	   << sqrt(systotl+pow(fSys[i][50].mult,2)+pow(fSys[i][51].mult,2))*fData[i]/100. << "   " 
	   << sqrt(systotr+pow(fSys[i][50].mult,2)+pow(fSys[i][51].mult,2))*fData[i]/100. << endl;
      */

    }
       
    /*
    //Defining covariance matrix for the specific bin
    double** covmat = new double*[artsys[bin]];
    for (int i = 0; i < artsys[bin]; i++) 
      covmat[i] = new double[artsys[bin]];

    //Reading Covariance Matrix
    for (int i = 0; i < 16; i++)
      getline(rCorr,line);
     
    for (int i = 0; i < artsys[bin]; i++){    
      for (int j = 0; j < artsys[bin]; j++) { 
        getline(rCorr,line);                       
        rCorr >> dum >> dum >> dum >> dum;
        rCorr >> covmat[i][j];
        covmat[j][i] *= stat[n+j+9]*stat[n+i+9];    
      }
      getline(rCorr,line);                       
    }

    //Generate artificial systematics
    double** syscor = new double*[artsys[bin]];
    for(int i = 0; i < artsys[bin]; i++)
      syscor[i] = new double[artsys[bin]];

    if(!genArtSys(artsys[bin],covmat,syscor))
    {
       cerr << " in " << fSetName << " : cannot generate artificial systematics" << endl;
       exit(-1);
    }

    //the first 9 points of each bin do not have any artificial systematic
    for (int i = n; i < n + 9; i++)
    {
      for (int l = realsys; l < fNSys; l++)  
      {
        fSys[i][l].add  = 1e-4;
        fSys[i][l].mult = fSys[i][l].add/fData[i]*1e2;
        fSys[i][l].type = MULT;  
        fSys[i][l].name = "CORR";
      }
    }

    // Copy the artificial systematics in the fSys matrix
    for (int i = n+9; i < n + ndata[bin]; i++)
    {
      for (int l = realsys; l < realsys + artsys[bin]; l++)  
      {
        fSys[i][l].add  = syscor[i-n-9][l-realsys]/fac;
        fSys[i][l].mult = fSys[i][l].add/fData[i]*1e2;
        fSys[i][l].type = MULT;  
        fSys[i][l].name = "CORR";
      }

      // Put the remaining systematics to 1e-4
      for (int l = realsys + artsys[bin]; l < fNSys; l++)
      {
        fSys[i][l].add  = 1e-4;
        fSys[i][l].mult = fSys[i][l].add/fData[i]*1e2;
        fSys[i][l].type = MULT;  
        fSys[i][l].name = "CORR";
      }
    }     

    
    for(int i = 0; i < artsys[bin]; i++) 
     delete[] syscor[i];
    delete[] syscor;
  
    for(int i = 0; i < artsys[bin]; i++) 
     delete[] covmat[i];
    delete[] covmat;
    */

    rS.close();
    rCorr.close();    
    n += ndata[bin];
  }    
}

