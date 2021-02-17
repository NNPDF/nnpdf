/*
Name_exp  : ATLAS_2JET_7TEV_R04
Reference : Measurement of dijet cross-sections in pp collisions at 7 TeV
            centre-of-mass energy using the ATLAS detector
ArXiv     : arXiv:1312.3524
Published : JHEP 1405 (2014) 059
Hepdata   : https://www.hepdata.net/record/ins1268975

Double-differential dijet cross-sections measured in pp collisions from the 
ATLAS experiment at 7 TeV, presented as functions of the dijet mass and half 
the rapidity separation of the two highest-pT jets. These measurements are 
obtained using data corresponding to an integrated luminosity of 4.5 fb-1, 
recorded by ATLAS detector in 2011. Jets are reconstructed using the anti-kT 
clustering algorithm for values of the jet radius parameter of 0.4, in a phase 
space region covering dijet mass up to 5 TeV and absolute rapidity up to 
|y|= 3.0.

The information on the experimental uncertainty is retrieved from hepdata.
Correlations between statistical uncertainties are taken into account.
Systematics uncertainties are fully corrrelated in dijet mass and rapidity,
statistical uncertainties are correlated only between different dijet mass bins 
of the same rapidity range. NData artificial systematic uncertainties are 
generated to take into account such correlations, starting from the covariance 
matrix obtained from hepdata.

For each rapidity range, three scenarios for the correlations between 
jet energy scale uncertainty components are given, denoted as nominal, 
stronger and weaker. They correspond to 132, 114 and 136 systematics 
respectively (+ lumi + qual). 
These scenarios are implemented as a unique set of systematics. 
The systematics names sys.name have to be changed consistently depending on the 
particular scenario one wants to consider. By default the nomial scenario is 
implemented, i.e. sys.name is set to SKIP for the stronger and weaker scenarios.

Total number of points        90
Total number of sys per point 472 (382 real + 90 art)

bin 1:   0 < |y| < 0.5
========================
points    21 
real sys  132(nominal) + 114(stronger) + 136(weaker) + lumi + qual
art sys   90

bin 2:   0.5 < |y| < 1.0
========================
points    21
real sys  132(nominal) + 114(stronger) + 136(weaker) + lumi + qual
art sys   90

bin 3:   1.0 < |y| < 1.5
========================
points    19
real sys  132(nominal) + 114(stronger) + 136(weaker) + lumi + qual
art sys   90 

bin 4:   1.5 < |y| < 2.0
========================
points    17
real sys  132(nominal) + 114(stronger) + 136(weaker) + lumi + qual
art sys   90 

bin 5:   2.0 < |y| < 2.5
========================
points    8
real sys  132(nominal) + 114(stronger) + 136(weaker) + lumi + qual
art sys   90 

bin 6:   2.5 < |y| < 3.0
========================
points    4 
real sys  132(nominal) + 114(stronger) + 136(weaker) + lumi + qual
art sys   90

Implemented by TG June 2019. 
*/


#include "ATLAS_2JET_7TEV_R04.h"

void ATLAS_2JET_7TEV_R04Filter::ReadData()
{

  fstream rS, rSs, rSw, rCorr;

  //bins specification
  int nbins = 6;
  int realsys = 384;
  std::vector<double> y = {0.25, 0.75, 1.25, 1.75, 2.25, 2.75 };
  std::vector<int> ndata = {21, 21, 19, 17, 8, 4};    

  int n = 0;                                //count total number of datapoints
  for(int bin=0; bin < nbins; bin++ )
  {

    //open data files corresponding to different scenarios
    string data_file_nominal = "/bin_" + to_string(bin+1) + "_nominal.dat";
    stringstream DataFile("");
    DataFile << dataPath() << "rawdata/" << fSetName << data_file_nominal;
    rS.open(DataFile.str().c_str(), ios::in);
    if (rS.fail()) {
      cerr << "Error opening data file " << DataFile.str() << endl;
      exit(-1);
    }

    string data_file_stronger = "/bin_" + to_string(bin+1) + "_stronger.dat";
    stringstream DataFile_s("");
    DataFile_s << dataPath() << "rawdata/" << fSetName << data_file_stronger;
    rSs.open(DataFile_s.str().c_str(), ios::in);
    if (rSs.fail()) {
      cerr << "Error opening data file " << DataFile_s.str() << endl;
      exit(-1);
    }

    string data_file_weaker = "/bin_" + to_string(bin+1) + "_weaker.dat";
    stringstream DataFile_w("");
    DataFile_w << dataPath() << "rawdata/" << fSetName << data_file_weaker;
    rSw.open(DataFile_w.str().c_str(), ios::in);
    if (rSw.fail()) {
      cerr << "Error opening data file " << DataFile_w.str() << endl;
      exit(-1);
    }
  
    string line;
    double m12min, m12max, dum;
    double sys1, sys2;
    for (int i = n; i < n + ndata[bin]; i++)
    {

      //============== nominal scenario
      getline(rS,line);                     
      istringstream lstream(line);          

      lstream >> dum >> m12min >> m12max;
      
      fKin1[i] = y[bin];                       // y, central value of the bin
      fKin2[i] = (m12min + m12max)*0.5*1000.;  // m12, central value of the bin, convert TeV to GeV
      fKin3[i] = 7000;                         // sqrt(s)

      lstream >> fData[i];                     // cross section [pb/TeV]


      //Read uncertainties

      //Statistical uncertainty (to be constructed from cov matrix)
      lstream >> dum;
      fStat[i] = 0.;                          

      //Read the 134 real systematics (nominal) 
      //note: sys1 and sys2 can be: both positive, both negative, negative and positive, positive and negative
      for (int k=0; k<66; k++)
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
	fSys[i][2*k].name = "CORR";       //set to SKIP to consider one of the other scenarios
	fSys[i][2*k].mult = sys1/sqrt(2.);
	fSys[i][2*k].add  = fSys[i][2*k].mult*fData[i]/100;

	fSys[i][2*k+1].type = MULT;   
	fSys[i][2*k+1].name = "CORR";     //set to SKIP to consider one of the other scenarios
	fSys[i][2*k+1].mult = sys2/sqrt(2.);
	fSys[i][2*k+1].add  = fSys[i][2*k+1].mult*fData[i]/100;
      }

      //Luminosity uncertainty (always simmetric, the same for all the scenarios)
      double sys;
      lstream >> sys >> dum;
      fSys[i][132].type = MULT;   
      fSys[i][132].name = "ATLASLUMI11";    
      fSys[i][132].mult = sys;
      fSys[i][132].add  = fSys[i][132].mult*fData[i]/100;

      //Qual uncertainty (always simmetric, the same for all the scenarios)
      lstream >> sys >> dum;
      fSys[i][133].type = MULT;   
      fSys[i][133].name = "CORR";    
      fSys[i][133].mult = sys1;
      fSys[i][133].add  = fSys[i][133].mult*fData[i]/100;

      //============== stronger scenario 
      getline(rSs,line);                     
      istringstream lstreams(line);
      lstreams >> dum >> dum >> dum >> dum >> dum;  // skip pt, xsec and stat values

      //Read the 114 real systematics (stronger), skip lumi and qual, read before  
      //note: sys1 and sys2 can be: both positive, both negative, negative and positive, positive and negative
      
      for (int k=0; k<57; k++)
      {
        lstreams >> sys1 >> sys2;
	  	  
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

        fSys[i][134+2*k].type = MULT;   
	fSys[i][134+2*k].name = "SKIP";   //set to CORR to consider this scenarios
	fSys[i][134+2*k].mult = sys1/sqrt(2.);
	fSys[i][134+2*k].add  = fSys[i][134+2*k].mult*fData[i]/100;

	fSys[i][134+2*k+1].type = MULT;   
	fSys[i][134+2*k+1].name = "SKIP"; //set to CORR to consider this scenarios
	fSys[i][134+2*k+1].mult = sys2/sqrt(2.);
	fSys[i][134+2*k+1].add  = fSys[i][134+2*k+1].mult*fData[i]/100;
      }

      //============== weaker
      getline(rSw,line);                     
      istringstream lstreamw(line);
      lstreamw >> dum >> dum >> dum >> dum >> dum;  // skip pt, xsec and stat values

      //Read the 136 real systematics (weaker), skip lumi and qual, read before  
      //note: sys1 and sys2 can be: both positive, both negative, negative and positive, positive and negative
      
      for (int k=0; k<68; k++)
      {
        lstreamw >> sys1 >> sys2;
	  	  
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

	fSys[i][248+2*k].type = MULT;   
	fSys[i][248+2*k].name = "SKIP";     //set to CORR to consider this scenarios
	fSys[i][248+2*k].mult = sys1/sqrt(2.);
	fSys[i][248+2*k].add  = fSys[i][248+2*k].mult*fData[i]/100;

	fSys[i][248+2*k+1].type = MULT;   
	fSys[i][248+2*k+1].name = "SKIP";   //set to CORR to consider this scenarios
	fSys[i][248+2*k+1].mult = sys2/sqrt(2.);
	fSys[i][248+2*k+1].add  = fSys[i][248+2*k+1].mult*fData[i]/100;
      }               
    }
        
    rS.close();
    rSs.close();
    rSw.close();
    n += ndata[bin];
  }

  //Defining covariance matrix for statistical uncertainties
  double** covmat = new double*[fNData];
  for (int i = 0; i < fNData; i++) 
    covmat[i] = new double[fNData];
    
  //Initialise covariance matrix
  for (int i = 0 ; i < fNData; i++)
    for (int j = 0; j < fNData; j++)
      covmat[i][j] = 0;    
  
  
  //Reading correlation coefficients
  for(int bin=0; bin < nbins; bin++ )
  {
    string line;
    string cov_file =  "/dijet_statcov/hepcov_R04_Eta" + to_string(bin) + ".txt";
    stringstream DataFileCorr("");
    DataFileCorr << dataPath() << "rawdata/" << fSetName << cov_file;
    rCorr.open(DataFileCorr.str().c_str(), ios::in);
      
    if (rCorr.fail()) 
    {
      cerr << "Error opening data file " << DataFileCorr.str() << endl;
      exit(-1);
    }  

    for (int i = 0; i < 12; i++)
      getline(rCorr,line);
    
    for (int i = 0; i < ndata[bin]; i++)
    {    
      double dum, cov;
      getline(rCorr,line);
      rCorr >> dum >> dum; 
      for (int j = 0; j < ndata[bin]; j++) 
      { 
	rCorr >> cov;
	if(bin==0)      
	  covmat[i][j] = cov;
	else if(bin==1) 
	  covmat[ndata[0]+i][ndata[0]+j] = cov;
	else if(bin==2) 
	  covmat[ndata[0]+ndata[1]+i][ndata[0]+ndata[1]+j] = cov;
	else if(bin==3) 
	  covmat[ndata[0]+ndata[1]+ndata[2]+i][ndata[0]+ndata[1]+ndata[2]+j] = cov;
	else if(bin==4)
	  covmat[ndata[0]+ndata[1]+ndata[2]+ndata[3]+i][ndata[0]+ndata[1]+ndata[2]+ndata[3]+j] = cov;
        else if(bin==5)
	  covmat[ndata[0]+ndata[1]+ndata[2]+ndata[3]+ndata[4]+i][ndata[0]+ndata[1]+ndata[2]+ndata[3]+ndata[4]+j] = cov;
      }
    }
    rCorr.close(); 
  }
  
  //Symmetrise covariance matrix (not necessary actually..it is already symmetric)
  for (int i = 0; i < fNData; i++)
    for (int j = 0; j < fNData; j++)
      if(i!=j) covmat[j][i]=covmat[i][j];
	   
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
  
  for(int i = 0; i < fNData; i++) 
  {
    delete[] syscor[i];
    delete[] covmat[i];
  }
  delete[] syscor;
  delete[] covmat;     
}
