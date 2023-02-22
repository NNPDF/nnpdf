/*
Name_exp  : CMS_1JET_8TEV
Reference : Measurement and QCD analysis of double-differential inclusive jet 
            cross-sections in pp collisions at s√= 8 TeV and ratios to 2.76 
            and 7 TeV
ArXiv     : arXiv:1609.05331
Published : JHEP 1703 (2017) 156
Hepdata   : https://www.hepdata.net/record/ins1487277

Measurement of the double-differential inclusive jet cross section as a 
function of the jet transverse momentum pT and the absolute jet rapidity |y|.
The measurement is based on the data sample collected from LHC proton-proton 
collisions from the CMS detector at s√= 8 TeV. It corresponds to an integrated 
luminosity of 19.7 fb−1. Jets are reconstructed using the anti-kT clustering 
algorithm with a size parameter of 0.7 in a phase space region covering jet pT 
from 74 GeV up to 2.5 TeV and jet absolute rapidity up to |y|= 3.0. 
The low-pT jet range between 21 and 74 GeV up to |y|= 4.7 is only partially 
implemented for consistency with the applgrids.

The information on experimental ucnertainties is retrieved from the intersection
between the hepdata entry (https://www.hepdata.net/record/ins1487277) and the 
XFitter analysis (http://www.hepforge.org/archive/xfitter/1609.05331.tar.gz).
Correlations between statistical uncertainties are taken into account. 
Statistical uncertainties are correlated only between different pT bins of 
the same rapidity range due to unfolding. Different rapidity
ranges are considered as uncorrelated among themselves.
NData artificial systematic uncertainties are generated
to take into account such correlations. Null correlations are
associated to the points below the nominal pT cut (pT<74 GeV) for which the
information on correlations is not available.
On top of these uncertainties, there are 52 sources of systematiic uncertainty:
- unfolding  uncertainties (1x2) include the uncertainty induced by the 
  JER parameterization in the unfolding procedure;
- JES uncertainties (24X2) include the sum of 24 independent sources of 
  uncertainty;
- uncorr uncertainty (1) include a 1% overall uncorrelated systematic;
- lumi uncertainty (1) include an overall 2.6% luminosity systematic. 
Asymmetric uncertainties are provided for unfolding and JES. They correspond to
upwards and downwards shifts of the modeling parameters, which can result in 
+,-; -,+; +,+; -,- pairs of uncertainties. Both elements of the pair are 
implemented as separate systematic uncertainties, with their sign, after a 
1/sqrt{2} rescaling. If the two variations in the pair have the same sign, 
only the largest (in absolute value) is retained, while the other is set to 
zero. This procedure is consistent with the analysis presented in the paper and
allows one to reproduce the total left (right) systematic uncertainties
quoted on hepdata, once all the resulting left (right) sources of nuisance are
added in quadrature. 
Note that the Xfitter files contain additional nonperturbative corrections
(a rescaling factor with its left and right error). These are implemented as an
extra systematic ucnertainty determined as the difference between the measured
data and the data rescaled by the provided factor (including its uncertainty).

bin 1:   0 < |y| < 0.5
========================
points    46 
real sys  50 + lumi + uncorr 
art sys   239
systot    291

bin 2:   0.5 < |y| < 1.0
========================
points    46 
real sys  50 + lumi + uncorr
art sys   37
systot    291

bin 3:   1.5 < |y| < 2.0
========================
points    45 
real sys  50 + lumi + uncorr 
art sys   239
systot    291

bin 4:   2.0 < |y| < 2.5
========================
points    41 
real sys  50 + lumi + uncorr 
art sys   239
systot    291

bin 5:   2.5 < |y| < 3.0
========================
points    34 
real sys  50 + lumi + uncorr 
art sys   239
systot    291

bin 6:   3.0 < |y| < 3.5
========================
points    27 
real sys  50 + lumi + uncorr
art sys   239
systot    291

tot points         239
tot sys per point  291

Implemented by TG May 2019. Revised by ERN May 2019.
*/


#include "CMS_1JET_8TEV.h"

void CMS_1JET_8TEVFilter::ReadData()
{

  //opening files
  fstream rS, rCorr;

  //bins specification
  int nbins = 6;
  std::vector<double> y = {0.25, 0.75, 1.25, 1.75, 2.25, 2.75};
  std::vector<int> ndata  = {46, 46, 45, 41, 34, 27};    
  std::vector<int> artsys = {37, 37, 36, 32, 25, 18};  //number of artificial systematics for each bin
  std::vector<double> stat(fNData);   

  int n = 0;                                //count total number of datapoints
  //const double fac = 1e7;                 //conversion factor from pb to 10 mub
  const int realsys=54;                     //number of real systematics (including NP theoretical uncertainties)

  for(int bin=0; bin < nbins; bin++ )
  {

    string data_file = "/CMS_8TeV_jets_Ybin" + to_string(bin+1) + ".dat";
    stringstream DataFile("");
    DataFile << dataPath() << "rawdata/" << fSetName << data_file;
    rS.open(DataFile.str().c_str(), ios::in);
    if (rS.fail()) {
      cerr << "Error opening data file " << DataFile.str() << endl;
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

      //1) Relative theoretical uncertainties due to nonperturbative corrections in the prediction
      double np_corr, np_corr_erp, np_corr_erm;
      lstream >> np_corr >> np_corr_erp >> np_corr_erm;

      fSys[i][0].mult = (np_corr*(1. + np_corr_erp/100.) - 1)/sqrt(2.)*100.;
      fSys[i][0].add  = fSys[i][0].mult*fData[i]/100;
      fSys[i][0].type = MULT;   
      fSys[i][0].name = "SKIP";
      fSys[i][1].mult = (np_corr*(1. + np_corr_erm/100.) - 1)/sqrt(2.)*100.;	
      fSys[i][1].add  = fSys[i][1].mult*fData[i]/100;
      fSys[i][1].type = MULT;   
      fSys[i][1].name = "SKIP";

      //2) Absolute statistical uncertainty
      lstream >> stat[i];
      fStat[i] = 0.;

      //3) Relative systematic uncertainties due to the unfolding procedure
      //note: sys1 always positive; sys2 always negative
      lstream >> sys1 >> sys2;

      fSys[i][2].type = MULT;   
      fSys[i][2].name = "CORR";
      fSys[i][2].mult = sys1/sqrt(2.);
      fSys[i][2].add  = fSys[i][0].mult*fData[i]/100;

      fSys[i][3].type = MULT;   
      fSys[i][3].name = "CORR";
      fSys[i][3].mult = sys2/sqrt(2.);
      fSys[i][3].add  = fSys[i][1].mult*fData[i]/100;

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

	  fSys[i][4+2*k].type = MULT;   
	  fSys[i][4+2*k].name = "CORR";
	  fSys[i][4+2*k].mult = sys1/sqrt(2.);
	  fSys[i][4+2*k].add  = fSys[i][2+2*k].mult*fData[i]/100;

	  fSys[i][4+2*k+1].type = MULT;   
	  fSys[i][4+2*k+1].name = "CORR";
	  fSys[i][4+2*k+1].mult = sys2/sqrt(2.);
	  fSys[i][4+2*k+1].add  = fSys[i][2+2*k+1].mult*fData[i]/100;
	}

      //5) Luminosity uncertainty
      double sys;
      lstream >> sys;
      fSys[i][52].type = MULT;   
      fSys[i][52].name = "CMSLUMI12";
      fSys[i][52].mult = sys;
      fSys[i][52].add  = fSys[i][50].mult*fData[i]/100;

      //6) Relative uncorrelated systematic uncertainty
      lstream >> dum >> sys;
      fSys[i][53].type = MULT;   
      fSys[i][53].name = "UNCORR";
      fSys[i][53].mult = sys;
      fSys[i][53].add  = fSys[i][51].mult*fData[i]/100;

      /*
      //Checking systematics
      //By uncommenting the following lines, one is able to recover the values
      //of the total asymmetric (upwards and downwards) systematic ucnertainties 
      //quoted on hepdata
      
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
          
    rS.close();
    n += ndata[bin];

  }

  //Defining covariance matrix for statistical uncertainties
  double** covmat = new double*[fNData];
  for (int i = 0; i < fNData; i++) 
    covmat[i] = new double[fNData];
      
  //Initialise covariance matrix
  for (int i = 0 ; i < fNData; i++)
    {
      for (int j = 0; j < fNData; j++)
	{
	  if(i==j) covmat[i][j] = stat[i]*stat[j];
	  else covmat[i][j] = 0.;
	}
    }
  
  //Reading correlation coefficients
  for(int bin=0; bin < nbins; bin++ )
    {
      string line;
      string cov_file =  "/CMS_8TeV_jets_Ybin" + to_string(bin+1) + "___CMS_8TeV_jets_Ybin" + to_string(bin+1) + ".dat";
      stringstream DataFileCorr("");
      DataFileCorr << dataPath() << "rawdata/" << fSetName << cov_file;
      rCorr.open(DataFileCorr.str().c_str(), ios::in);
      
    if (rCorr.fail()) {
      cerr << "Error opening data file " << DataFileCorr.str() << endl;
      exit(-1);
    }  

    for (int i = 0; i < 16; i++)
      getline(rCorr,line);
    
    for (int i = 0; i < artsys[bin]; i++)
      {    
	for (int j = 0; j < artsys[bin]; j++) 
	  { 
	    double dum, rho;
	    getline(rCorr,line);                       
	    rCorr >> dum >> dum >> dum >> dum;
	    rCorr >> rho;

	    if(bin==0)      
	      covmat[9+i][9+j] = stat[9+i] * stat[9+j] * rho;
	    else if(bin==1) 
	      covmat[ndata[0]+9+i][ndata[0]+9+j] = stat[ndata[0]+9+i] * stat[ndata[0]+9+j] * rho;
	    else if(bin==2) 
	      covmat[ndata[0]+ndata[1]+9+i][ndata[0]+ndata[1]+9+j] = stat[ndata[0]+ndata[1]+9+i] * stat[ndata[0]+ndata[1]+9+j] * rho;
	    else if(bin==3) 
	      covmat[ndata[0]+ndata[1]+ndata[2]+9+i][ndata[0]+ndata[1]+ndata[2]+9+j] = stat[ndata[0]+ndata[1]+ndata[2]+9+i] * stat[ndata[0]+ndata[1]+ndata[2]+9+j] * rho;
	    else if(bin==4)
	      covmat[ndata[0]+ndata[1]+ndata[2]+ndata[3]+9+i][ndata[0]+ndata[1]+ndata[2]+ndata[3]+9+j] = stat[ndata[0]+ndata[1]+ndata[2]+ndata[3]+9+i] * stat[ndata[0]+ndata[1]+ndata[2]+ndata[3]+9+j] * rho;
	    else if(bin==5)
	      covmat[ndata[0]+ndata[1]+ndata[2]+ndata[3]+ndata[4]+9+i][ndata[0]+ndata[1]+ndata[2]+ndata[3]+ndata[4]+9+j] = stat[ndata[0]+ndata[1]+ndata[2]+ndata[3]+ndata[4]+9+i] * stat[ndata[0]+ndata[1]+ndata[2]+ndata[3]+ndata[4]+9+j] * rho;
	  }
	getline(rCorr,line);                       
      }

    rCorr.close(); 
    }   
  
  //Symmetrise covariance matrix
  for (int i = 0; i < fNData; i++)
    {
      for (int j = 0; j < fNData; j++)
	{
	  if(i!=j)
	    covmat[j][i]=covmat[i][j];
	}
    }
  
  //Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];
  
  if(!genArtSys(fNData,covmat,syscor))
    {
      cerr << " in " << fSetName << " : cannot generate artificial systematics" << endl;
      exit(-1);
    }
  
  //Assign artificial systematics to data points consistently
  for (int i = 0; i < fNData; i++)
    {
      for (int l = realsys; l < fNSys; l++)
	{
	  fSys[i][l].add = syscor[i][l-realsys];
	  fSys[i][l].mult = fSys[i][l].add/fData[i]*1e2;
	  fSys[i][l].type = ADD;
	  fSys[i][l].name = "CORR";
	}
    }

  /*
  //Check statistical uncertainties
  for (int i=0; i<fNData; i++)
  {
  double totstat=0.;
  for (int l=realsys; l<fNSys; l++)
  {
  totstat += pow(fSys[i][l].add,2);
  }
  cout << i << "   " << sqrt(totstat) << "   " << stat[i] << endl;
  }
  */
  
  for(int i = 0; i < fNData; i++) 
    {
      delete[] syscor[i];
      delete[] covmat[i];
    }
  delete[] syscor;
  delete[] covmat;
  
}

