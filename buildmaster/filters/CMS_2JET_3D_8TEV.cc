/*
Name_exp  : CMS_2JET_3D_8TEV
Reference : Measurement of the triple-differential dijet cross section 
            in proton-proton collisions at s√=8TeV and constraints 
            on parton distribution functions 
ArXiv     : arXiv:1705.02628
Published : Eur.Phys.J. C77 (2017) no.11, 746
Hepdata   : https://www.hepdata.net/record/ins1598460

Measurement of the triple-differential dijet cross section as a function of:
- average transverse momentum ptavg = (pt1 + pt2)/2
- half the rapidity separation ys = |y1 - y2|/2
- boost of dijet system yb = |y1 + y2|/2
where 1 and 2 are the two leading jets in the event.

The measurement is based on the data sample collected from LHC proton-proton 
collisions from the CMS detector at s√= 8 TeV. It corresponds to an integrated 
luminosity of 19.7 fb−1. 

Jets are reconstructed using the anti-kT clustering algorithm with radius 0.7.
Only events with at least two jets up to an absolute rapidity of |y| = 5.0 
are selected and the two jets leading in pT are required to have transverse 
momenta greater than 50 GeV and |y| < 3.0.

STATISTICAL UNCERTAINTIES: 
Correlations between statistical uncertainties are taken into account. 
Statistical uncertainties are correlated only between different pT bins of 
the same rapidity range due to unfolding. Different rapidity ranges 
are considered as uncorrelated among themselves.
NData artificial systematic uncertainties are generated
to take into account such correlations. 

SYSTEMATIC UNCERTAINTIES:
On top of these uncertainties, there are 28 sources of systematic uncertainty.
According to HepData:
- "uncorr" uncertainty (1) include a 1% overall uncorrelated systematic
- "lumi" uncertainty (1) include an overall 2.6% luminosity systematic. 
(identified by the flag CMSLUMI12, which is the same for e.g. CMS_1JET_8TEV,
since the integrated luminosity is the same for the two dataset 19.7 fb−1)
- "jererr" and "nongaussiantails" (2) are related to JER.
- JES uncertainties (24);
All these uncertainties are provided as already symmetrised.

NP UNCERTAINTIES
Nonperturbative (theoretical) corrections are implemented in the form of a
set of two additional correlated uncertainties (left and right) determined as
the difference between the central data point and the datapoint rescaled by the
correction factor (+- its error) provided by the experimentalists.

bin 1:   
0.0 < ys < 1.0
0.0 < yb < 1.0
========================
points    31 
real sys  26 + lumi + uncorr 
art sys   122
systot    150

bin 2:
1.0 < ys < 2.0
0.0 < yb < 1.0
========================
points    26 
real sys  26 + lumi + uncorr 
art sys   122
systot    150

bin 3:
2.0 < ys < 3.0
0.0 < yb < 1.0
========================
points    14 
real sys  26 + lumi + uncorr 
art sys   122
systot    150

bin 4:   
0.0 < ys < 1.0
1.0 < yb < 2.0
========================
points    23 
real sys  26 + lumi + uncorr 
art sys   122
systot    150

bin 5:
1.0 < ys < 2.0
1.0 < yb < 2.0
========================
points    17 
real sys  26 + lumi + uncorr 
art sys   122
systot    150

bin 6:
0.0 < ys < 1.0
2.0 < yb < 3.0
========================
points    11 
real sys  26 + lumi + uncorr 
art sys   122
systot    150

tot points         122            
tot sys per point  150

Implemented by GS July 2019.
*/


#include "CMS_2JET_3D_8TEV.h"


int getIndex(const int & idx, const int & bin, const std::vector<int> & ndata)
{  
  int res = 0;
  
  if(bin >= 0) res += idx;
  if(bin >= 1) res += ndata[0];
  if(bin >= 2) res += ndata[1];
  if(bin >= 3) res += ndata[2];
  if(bin >= 4) res += ndata[3];
  if(bin >= 5) res += ndata[4];
  
  return res;
}


void CMS_2JET_3D_8TEVFilter::ReadData()
{

  //opening files
  fstream rS, rNP, rCorr;

  //bins specification
  int nbins = 6;
  std::vector<int> ndata  = {31, 26, 14, 23, 17, 11};
  //number of artificial systematics for each bin
  //(equal to the number of points in the bin)
  std::vector<int> artsys = {31, 26, 14, 23, 17, 11};  
  std::vector<double> stat(fNData);   

  std::vector<double> ys = {0.5, 1.5, 2.5, 0.5, 1.5, 0.5};
  std::vector<double> yb = {0.5, 0.5, 0.5, 1.5, 1.5, 2.5};
  
  int n = 0;                                //count total number of datapoints
  //const double fac = 1e7;                   //conversion factor from pb to 10 mub
  const int realsys=30;                     //number of real systematics

  for(int bin=0; bin < nbins; bin++ ) {
    
    //Data file
    string data_file = "/CMS_2JET_3D_Ybin" + to_string(bin+1) + ".dat";
    stringstream DataFile("");
    DataFile << dataPath() << "rawdata/" << fSetName << data_file;
    rS.open(DataFile.str().c_str(), ios::in);
    if (rS.fail()) {
      cerr << "Error opening data file " << DataFile.str() << endl;
      exit(-1);
    }

    //NP corrections file
    string NP_file = "/ew_np_corrections/Table" + to_string(bin+19) + ".csv";
    stringstream NPFile("");
    NPFile << dataPath() << "rawdata/" << fSetName << NP_file;
    rNP.open(NPFile.str().c_str(), ios::in);
    if (rNP.fail()) {
      cerr << "Error opening data file " << NPFile.str() << endl;
      exit(-1);
    }

    string line, lineNP;
    double ptavg, ptavg_low, ptavg_high;
    
    double dum; //dummy variable
    
    //skip the first lines
    //In HepData .csv these lines are 12
    for (int i = 0; i < 12; i++) {
      getline(rS,line);
      getline(rNP,lineNP);
    }
    
    for (int i = n; i < n + ndata[bin]; i++) {

      getline(rS,line);
      istringstream lstream(line);
      lstream >> ptavg >> ptavg_low >> ptavg_high;

      //---------------------------------------------------------
      fKin1[i] = ys[bin];         // ys, central value of the bin
      fKin2[i] = pow(ptavg,2.);   // ptavg**2, where ptavg is central value of the bin
      fKin3[i] = yb[bin];         // yb, central value of the bin
      // fKin3[i] = 8000;           // sqrt(s) in GeV
      //---------------------------------------------------------

      lstream >> fData[i];       // cross section [pb/GeV]

      // uncorrelated systematic uncertainty
      double sys;
      lstream >> sys >> dum;
      fSys[i][0].type = MULT;   
      fSys[i][0].name = "UNCORR";
      fSys[i][0].mult = sys;
      fSys[i][0].add  = fSys[i][0].mult*fData[i]/100;

      // relative statistical uncertainty
      double multstat;
      lstream >> multstat >> dum;
      // saving the absolute value ...
      stat[i] = multstat*fData[i]/100;
      // ... and setting to zero in CommonData
      fStat[i] = 0.;

      // relative systematic uncertainties due to the unfolding procedure
      // + luminosity uncertainty
      double jer1, lum, jer2;

      lstream >> jer1 >> dum >> lum >> dum >> jer2 >> dum;

      fSys[i][1].type = MULT;   
      fSys[i][1].name = "CORR";
      fSys[i][1].mult = jer1;
      fSys[i][1].add  = fSys[i][1].mult*fData[i]/100;

      fSys[i][2].type = MULT;   
      fSys[i][2].name = "CMSLUMI12";
      fSys[i][2].mult = lum;
      fSys[i][2].add  = fSys[i][2].mult*fData[i]/100;
    
      fSys[i][3].type = MULT;   
      fSys[i][3].name = "CORR";
      fSys[i][3].mult = jer2;
      fSys[i][3].add  = fSys[i][3].mult*fData[i]/100;

      // relative sytematic JES uncertainties, obtained from 24 independent sources of uncertainty
      for (int k=0; k<24; k++) {
	
	double sys;
	lstream >> sys >> dum;
	
	fSys[i][4+k].type = MULT;   
	fSys[i][4+k].name = "CORR";
	fSys[i][4+k].mult = sys;
	fSys[i][4+k].add  = fSys[i][4+k].mult*fData[i]/100;
      }

    }

    for (int i = n; i < n + ndata[bin]; i++) {

      getline(rNP,lineNP);
      istringstream lstream(lineNP);
      double np_corr, np_corr_erp, np_corr_erm;
      char comma, perc;

      lstream >> ptavg        >> comma 
	      >> ptavg_low    >> comma
	      >> ptavg_high   >> comma
	      >> np_corr      >> comma
	      >> np_corr_erp >> perc >> comma
	      >> np_corr_erm >> perc;

      fSys[i][28].mult = (np_corr*(1. + np_corr_erp/100.) - 1)/sqrt(2.)*100.;
      fSys[i][28].add  = fSys[i][28].mult*fData[i]/100;
      fSys[i][28].type = MULT;   
      fSys[i][28].name = "SKIP";
      fSys[i][29].mult = (np_corr*(1. + np_corr_erm/100.) - 1)/sqrt(2.)*100.;	
      fSys[i][29].add  = fSys[i][29].mult*fData[i]/100;
      fSys[i][29].type = MULT;   
      fSys[i][29].name = "SKIP";

    }
      
    rS.close();
    rNP.close();
    n += ndata[bin];
  }
  
  //Defining covariance matrix for statistical uncertainties
  double** covmat = new double*[fNData];
  for (int i = 0; i < fNData; i++) 
    covmat[i] = new double[fNData];

  //Initialise covariance matrix
  for (int i = 0 ; i < fNData; i++) {
    for (int j = 0; j < fNData; j++) {
      if(i==j) covmat[i][j] = stat[i]*stat[j];
      else covmat[i][j] = 0.;
    }
  }
    
  //Reading correlation coefficients
  for(int bin=0; bin < nbins; bin++ ) {

    string line;
    string cov_file =  "/CMS_2JET_3D_Ybin" + to_string(bin+1) + "_corr.dat";
    stringstream DataFileCorr("");
    DataFileCorr << dataPath() << "rawdata/" << fSetName << cov_file;
    rCorr.open(DataFileCorr.str().c_str(), ios::in);
    
    if (rCorr.fail()) {
      cerr << "Error opening data file " << DataFileCorr.str() << endl;
      exit(-1);
    }  

    for (int i = 0; i < 12; i++) {
      getline(rCorr,line);
      // cout << line << endl;
    }

    for (int i = 0; i < artsys[bin]; i++) {  
      for (int j = 0; j < artsys[bin]; j++) {
	double dum, rho;
	
	getline(rCorr,line);                       
	istringstream lstream(line);

	lstream >> dum >> dum >> dum >> dum >> dum >> dum;
	lstream >> rho;

	covmat[getIndex(i, bin, ndata)][getIndex(j, bin, ndata)] = stat[getIndex(i, bin, ndata)] * stat[getIndex(j, bin, ndata)] * rho;
	
      }
    }

    rCorr.close();
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
  for (int i = 0; i < fNData; i++) {
    for (int l = realsys; l < fNSys; l++) {
      fSys[i][l].add = syscor[i][l-realsys];
      fSys[i][l].mult = fSys[i][l].add/fData[i]*1e2;
      fSys[i][l].type = ADD;
      fSys[i][l].name = "CORR";
    }
  }

  /*
  //Check statistical uncertainties
  for (int i=0; i<fNData; i++) {
    double totstat=0.;
    for (int l=realsys; l<fNSys; l++) {
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
