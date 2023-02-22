/*
Experiment: CERN-LHC-ATLAS (ATLAS)
Preprinted as CERN-PH-EP-2015-239
Archived as: ARXIV:1511.04716
Auxiliary Material: https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/TOPQ-2015-06/

Record in: INSPIRE
Record in: CERN Document Server
Record in: HEPData (new site in development)

Statistical covariance matrix provided by M.A.Owen June 2016

Description of the measurement
Measurements of normalized differential cross-sections of top-quark 
pair production are presented as a function of the top-quark, 
ttbar system and event-level kinematic observables in proton-proton 
collisions at a centre-of-mass energy of sqrt(s)=8 TeV. 
The observables have been chosen to emphasize the ttbar production process 
and to be sensitive to effects of initial- and final-state radiation, 
to the different parton distribution functions, and to non-resonant 
processes and higher-order corrections. 
The dataset corresponds to an integrated luminosity of 20.3 fb^-1, 
recorded in 2012 with the ATLAS detector at the CERN Large Hadron Collider. 
Events are selected in the lepton+jets channel, requiring exactly 
one charged lepton and at least four jets with at least two of the jets 
tagged as originating from a b-quark. The measured spectra are corrected 
for detector effects and are compared to several Monte Carlo simulations. 
 
Description of the buildmaster implementation
Absolute and normalised cross sections for the distributions (lepton+jets 
channel) differential in the following variables are implemented:
1) top quark transverse momentum;       
2) top quark pair transverse momentum;  
3) top quark absolute rapidity;                  
4) top quark absolute pair rapidity;             
5) top quark pair invariant mass.       

Raw data and full breakdown of systematic uncertainties are from HepData:
http://hepdata.cedar.ac.uk/h8test/view/ins1404878
1) TABS 29-30 HepData; TABS 26-25 auxiliary material
2) TABS 25-26 HepData; TABS 32-31 auxiliary material
3) TABS 31-32 HepData; TABS 28-27 auxiliary material
4) TABS 27-28 HepData; TABS 34-33 auxiliary material
5) TABS 23-24 HepData; TABS 30-29 auxiliary material

Notes:
1) The data are provided for both absolute and normalised differential 
   distributions with 56 sources of systematic uncertainty.
2) Absolute rapidity distributions are turned into rapidity distributions:
   the central value is divided by a factor of 2 and the uncertainty is
   divided by a factor sqrt(2). This way both ATLAS rapidity distributions 
   are CMS-like. 
3) All systematic uncertainties are assumed to be multiplicative.
4) Custom uncertainty descriptions are assumed to allow for cross-correlations
   among the five differential distributions. 
5) There is an additional multiplicative systematic (luminosity) uncertainty
   +-2.8% on the unnormalized data set which is not listed in the rawdata files.
6) The 57 sources of systematics are:
   1   single top cross section
   2   W+jets scale factors
   3   fake lept. alternate real CR mu+jets
   4   eta intercalibration model (JES)
   5   effective stat. NP set 3 (JES)
   6   fake lept. alternate real CR e+jets
   7   pile-up offset mu (JES)
   8   fake lept. MC stat e+jets
   9   fake lept. MC stat mu+jets
   10  ETmiss soft jet scale
   11  alternate hard scattering model
   12  effective stat. NP set 2 (JES)
   13  electron energy scale
   14  punch-through (JES)
   15  pile-up offset NPV (JES)
   16  lepton reconstruction efficiency
   17  pile-up offset pT (JES)
   18  jet energy resolution
   19  light-jet tagging efficiency
   20  fake lept. alternate fake CR e+jets
   21  fake lept. alternate parametrization mu+jets
   22  jet reconstruction efficiency
   23  c-quark tagging efficiency
   24  diboson cross section
   25  electron energy resolution
   26  luminosity
   27  flavour composition (JES)
   28  effective detector NP set 2 (JES)
   29  effective detector NP set 3 (JES)
   30  jet vertex fraction
   31  lepton trigger efficiency
   32  b-tagged jet energy scale (JES)
   33  muon momentum scale
   34  single-particle  high pT (JES)
   35  ETmiss soft jet resolution
   36  effective detector NP set 1 (JES)
   37  fake lepton alternate parametrization e+jets
   38  effective stat. NP set 1 (JES)
   39  muon (ID) momentum resolution
   40  parton distribution function
   41  ISR/FSR + scale
   42  Z+jets cross section
   43  alternate parton shower model
   44  flavour response
   45  fake lept. alternate fake CR mu+jets
   46  muon momentum resolution
   47  effective model NP set 3 (JES)
   48  lepton identification efficiency
   49  effective mixed NP set 2 (JES)
   50  effective mixed NP set 1 (JES)
   51  b-quark tagging efficiency
   52  pile-up offset rho topology (JES)
   53  effective model NP set 4 (JES)
   54  Monte Carlo sample statistics
   55  effective model NP set 4 (JES) 
   56  effective model NP set 1 (JES)
   57  additional luminosity
*/

#include "ATLASTOPDIFF.h"
#include <math.h>

//A - NORMALISED distributions

//1) Distribution differential in top quark transverse momentum
void  ATLASTOPDIFF8TEVTPTNORMFilter::ReadData()
{
  // Opening files
  fstream f1, f2;

  //Central values, statistical and systematic uncertainties  
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/" << fSetName << "/" << fSetName << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Statistical covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/" << fSetName << "/top_pt_covariance_rel.dat";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read statistical covariance matrix
  string line;
  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
  {
    covmat[i] = new double[fNData];
    getline(f2,line);
    istringstream lstream(line);
    for(int j=0; j<fNData; j++)
    {
      lstream >> covmat[i][j];
    }
  }

  //Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
    {
      cerr << " in " << fSetName << endl;
      exit(-1);
    }
  
  //Starting filter
  for(int i=0; i<5; i++)
    {
      getline(f1,line);
    }
  
  for(int i=0; i<fNData;i++)
    {
      double pt_top, ddum;
      getline(f1,line);
      istringstream lstream(line);
      lstream >> pt_top >> ddum >> ddum; 

      fKin1[i] = pt_top;       //P_T^(top)
      fKin2[i] = Mt*Mt;       
      fKin3[i] = 8000;         //sqrt(s)

      lstream >> fData[i];     //differential distribution
      lstream >> fStat[i];     //its statistical uncertainty
      fStat[i] = 0.;
      lstream >> ddum >> ddum >> ddum >> ddum >> ddum >> ddum;

      double shift = 0.;
      
      //Artificial systematics
      for(int i=0; i<fNData; i++)
	{
	  for(int j=0; j<fNData; j++)
	    {
	      fSys[i][j].add  = syscor[i][j];
	      fSys[i][j].mult = fSys[i][j].add*100/fData[i];
	      fSys[i][j].type = MULT;
	      fSys[i][j].name = "CORR";
	    }
	}  
      
      //Real systematics
      for(int j=fNData; j<fNSys; j++)
	{
	  double sys1, sys2, right, left;
	  double stmp, dtmp;
	  string sysdescr;
	  
	  ostringstream id;
	  id << j;
	  sysdescr = "CORR";
	  
	  lstream >> sys1 >> sys2;
	  if(sys1<0) {right=sys2; left=sys1;}
	  else {right=sys1; left=sys2;}

	  //convert to relative percentage values
	  right = right/fData[i]*100;  
	  left = left/fData[i]*100;
	  symmetriseErrors(right,left,&stmp,&dtmp);
	  
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = sysdescr;
	  fSys[i][j].mult = stmp;
	  fSys[i][j].add  = fSys[i][j].mult*fData[i]/100;

	  shift += dtmp;
	}

      fData[i]*=(1.0 + shift*0.01); //Shift from asymmetric errors

    }  

  for(int i = 0; i < fNData; i++) 
    {
      delete[] syscor[i];
      delete[] covmat[i];
    }
  delete[] syscor;
  delete[] covmat;
  
  f1.close();
  f2.close();
}

//==================================================================

//2) Distribution differential in top quark pair transverse momentum
void  ATLASTOPDIFF8TEVTTPTNORMFilter::ReadData()
{
  // Opening files
  fstream f1;

  //Central values, statistical and systematic uncertainties  
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/" << fSetName << "/" << fSetName << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Starting filter
  string line;
  for(int i=0; i<5; i++)
    {
      getline(f1,line);
    }
  
  for(int i=0; i<fNData;i++)
    {
      double pt_top, ddum;
      getline(f1,line);
      istringstream lstream(line);
      lstream >> pt_top >> ddum >> ddum; 

      fKin1[i] = pt_top;       //P_T^(top)
      fKin2[i] = Mt*Mt;       
      fKin3[i] = 8000;         //sqrt(s)

      lstream >> fData[i];     //differential distribution
      lstream >> fStat[i];     //its statistical uncertainty
      lstream >> ddum >> ddum >> ddum >> ddum >> ddum >> ddum;
      
      double shift = 0.;

      for(int j=0; j<fNSys; j++)
	{
	  double sys1, sys2, right, left;
	  double stmp, dtmp;
	  string sysdescr;
	  
	  ostringstream id;
	  id << j;
	  sysdescr = "CORR";
	  
	  lstream >> sys1 >> sys2;
	  if(sys1<0) {right=sys2; left=sys1;}
	  else {right=sys1; left=sys2;}

	  //convert to relative percentage values
	  right = right/fData[i]*100;  
	  left = left/fData[i]*100;
	  symmetriseErrors(right,left,&stmp,&dtmp);
	  
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = sysdescr;
	  fSys[i][j].mult = stmp;
	  fSys[i][j].add  = fSys[i][j].mult*fData[i]/100;

	  shift += dtmp;

	}

       fData[i]*=(1.0 + shift*0.01); //Shift from asymmetric errors
      
    }  
  
  f1.close();
}

//==============================================================

//3) Distribution differential in top quark rapidity
void  ATLASTOPDIFF8TEVTRAPNORMFilter::ReadData()
{
  // Opening files
  fstream f1, f2;

  //Central values, statistical and systematic uncertainties  
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/" << fSetName << "/" << fSetName << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Statistical covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/" << fSetName << "/top_absy_covariance_rel.dat";
  f2.open(covfile.str().c_str(), ios::in);
  
  if (f2.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read statistical covariance matrix
  string line;
  double** covmat = new double*[fNData/2];
  for(int i=0; i<fNData/2; i++)
    {
      covmat[i] = new double[fNData/2];
      getline(f2,line);
      istringstream lstream(line);
      for(int j=0; j<fNData/2; j++)
	{
	  lstream >> covmat[i][j];
	}
    }

  //Generate artificial systematics
  double** syscor = new double*[fNData/2];
  for(int i = 0; i < fNData/2; i++)
    syscor[i] = new double[fNData/2];
  
  if(!genArtSys(fNData/2,covmat,syscor))
    {
      cerr << " in " << fSetName << endl;
      exit(-1);
    }  

  //Starting filter
  for(int i=0; i<5; i++)
    {
      getline(f1,line);
    }
  
  for(int i=0; i<fNData/2;i++)
    {
      double pt_top, ddum, datum, stat;
      getline(f1,line);
      istringstream lstream(line);
      lstream >> pt_top >> ddum >> ddum; 

      fKin1[i+fNData/2]   = pt_top;          //yt
      fKin1[fNData/2-1-i] = -1.0*pt_top;
      fKin2[i+fNData/2]   = Mt*Mt;       
      fKin2[fNData/2-1-i] = Mt*Mt;
      fKin3[i+fNData/2]   = 8000;            //sqrt(s)
      fKin3[fNData/2-1-i] = 8000;        

      lstream >> datum;
      fData[i+fNData/2]   = datum/2.0;       //differential distribution 
      fData[fNData/2-1-i] = datum/2.0;

      lstream >> stat;
      fStat[i+fNData/2]   = stat/sqrt(2.0);  //statistical uncertainty
      fStat[fNData/2-1-i] = stat/sqrt(2.0);

      //fStat[i+fNData/2] = 0.;
      //fStat[fNData/2-1-i] = 0.;

      lstream >> ddum >> ddum >> ddum >> ddum >> ddum >> ddum;
      
      double shift = 0.;

      /*
      //Artificial systematics
      for(int i=0; i<fNData/2; i++)
	{
	  for(int j=0; j<fNData/2; j++)
	    { 
	      //const double eps=0.0000001;

	      fSys[i+fNData/2][fNData/2-1-j].add  = syscor[i][j]/sqrt(2.0)/2.;
	      fSys[i+fNData/2][fNData/2-1-j].mult = fSys[i+fNData/2][fNData/2-1-j].add*100/fData[i+fNData/2];
	      fSys[i+fNData/2][fNData/2-1-j].type = MULT;
	      fSys[i+fNData/2][fNData/2-1-j].name = "CORR";

	      fSys[i+fNData/2][j+fNData/2].add  = syscor[i][j]/sqrt(2.0)/2.;
	      fSys[i+fNData/2][j+fNData/2].mult = fSys[i+fNData/2][j+fNData/2].add*100/fData[i+fNData/2];
	      fSys[i+fNData/2][j+fNData/2].type = MULT;
	      fSys[i+fNData/2][j+fNData/2].name = "CORR";
	      
	      fSys[fNData/2-1-i][fNData/2-1-j].add  = syscor[i][j]/sqrt(2.0)/2.;
	      fSys[fNData/2-1-i][fNData/2-1-j].mult = fSys[fNData/2-1-i][fNData/2-1-j].add*100/fData[fNData/2-1-i];
	      fSys[fNData/2-1-i][fNData/2-1-j].type = MULT;
	      fSys[fNData/2-1-i][fNData/2-1-j].name = "CORR";

	      fSys[fNData/2-1-i][j+fNData/2].add  = syscor[i][j]/sqrt(2.0)/2.;
	      fSys[fNData/2-1-i][j+fNData/2].mult = fSys[fNData/2-1-i][j+fNData/2].add*100/fData[fNData/2-1-i];
	      fSys[fNData/2-1-i][j+fNData/2].type = MULT;
	      fSys[fNData/2-1-i][j+fNData/2].name = "CORR";
	    }
	}  
      */

      //Real systematics
      for(int j=0; j<fNSys; j++)
	{
	  double sys1, sys2, right, left;
	  double stmp, dtmp;
	  string sysdescr;

	  ostringstream id;
	  id << j;
	  sysdescr = "CORR";
	  
	  lstream >> sys1 >> sys2; 
	  if(sys1<0) {right=sys2; left=sys1;}
	  else {right=sys1; left=sys2;}

	  //convert to relative percentage values
	  right = right/datum*100;  
	  left = left/datum*100;
	  symmetriseErrors(right,left,&stmp,&dtmp);
	  
	  fSys[i+fNData/2][j].type = MULT;
	  fSys[i+fNData/2][j].name = sysdescr;
	  fSys[i+fNData/2][j].mult = stmp;
	  fSys[i+fNData/2][j].add  = fSys[i+fNData/2][j].mult*fData[i+fNData/2]/100;

	  fSys[fNData/2-1-i][j].type = MULT;
	  fSys[fNData/2-1-i][j].name = sysdescr;
	  fSys[fNData/2-1-i][j].mult = stmp;
	  fSys[fNData/2-1-i][j].add  = fSys[fNData/2-1-i][j].mult*fData[fNData/2-1-i]/100;

	  shift += dtmp;

	}
          
       fData[i+fNData/2]*=(1.0 + shift*0.01);   //Shift from asymmetric errors
       fData[fNData/2-1-i]*=(1.0 + shift*0.01); //Shift from asymmetric errors

    }  

  for(int i = 0; i < fNData/2; i++) 
    {
      delete[] syscor[i];
      delete[] covmat[i];
    }
  delete[] syscor;
  delete[] covmat;
  
  f1.close();
  f2.close();
}


//=================================================================

//4) Distribution differential in top quark pair rapidity
void  ATLASTOPDIFF8TEVTTRAPNORMFilter::ReadData()
{
  // Opening files
  fstream f1, f2;

  //Central values, statistical and systematic uncertainties  
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/" << fSetName << "/" << fSetName << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Statistical covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/" << fSetName << "/ttbar_absy_covariance_rel.dat";
  f2.open(covfile.str().c_str(), ios::in);
  
  if (f2.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Read statistical covariance matrix
  string line;
  double** covmat = new double*[fNData/2];
  for(int i=0; i<fNData/2; i++)
    {
      covmat[i] = new double[fNData/2];
      getline(f2,line);
      istringstream lstream(line);
      for(int j=0; j<fNData/2; j++)
	{
	  lstream >> covmat[i][j];
	}
    }
  
  //Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData/2; i++)
    syscor[i] = new double[fNData/2];
  
  if(!genArtSys(fNData/2,covmat,syscor))
    {
      cerr << " in " << fSetName << endl;
      exit(-1);
    }  
  
  //Starting filter
  for(int i=0; i<5; i++)
    {
      getline(f1,line);
    }
  
  for(int i=0; i<fNData/2;i++)
    {
      double pt_top, ddum, datum, stat;
      getline(f1,line);
      istringstream lstream(line);
      lstream >> pt_top >> ddum >> ddum; 

      fKin1[i+fNData/2]   = pt_top;          //P_T^(top)
      fKin1[fNData/2-1-i] = -1.0*pt_top;
      fKin2[i+fNData/2]   = Mt*Mt;       
      fKin2[fNData/2-1-i] = Mt*Mt;
      fKin3[i+fNData/2]   = 8000;            //sqrt(s)
      fKin3[fNData/2-1-i] = 8000;        

      lstream >> datum;
      fData[i+fNData/2]   = datum/2.0;       //differential distribution 
      fData[fNData/2-1-i] = datum/2.0;

      lstream >> stat;
      fStat[i+fNData/2]   = stat/sqrt(2.0);  //statistical uncertainty
      fStat[fNData/2-1-i] = stat/sqrt(2.0);
      //fStat[i+fNData/2] = 0.;
      //fStat[fNData/2-1-i] = 0.;

      lstream >> ddum >> ddum >> ddum >> ddum >> ddum >> ddum;
      
      double shift = 0.;

      /*
     //Artificial systematics
      for(int i=0; i<fNData/2; i++)
	{
	  for(int j=0; j<fNData/2; j++)
	    {
	      //const double eps=0.0000001;

	      fSys[i+fNData/2][fNData/2-1-j].add  = syscor[i][j]/sqrt(2.0)/2.;
	      fSys[i+fNData/2][fNData/2-1-j].mult = fSys[i+fNData/2][fNData/2-1-j].add*100/fData[i+fNData/2];
	      fSys[i+fNData/2][fNData/2-1-j].type = MULT;
	      fSys[i+fNData/2][fNData/2-1-j].name = "CORR";

	      fSys[i+fNData/2][j+fNData/2].add  = syscor[i][j]/sqrt(2.0)/2.;
	      fSys[i+fNData/2][j+fNData/2].mult = fSys[i+fNData/2][j+fNData/2].add*100/fData[i+fNData/2];
	      fSys[i+fNData/2][j+fNData/2].type = MULT;
	      fSys[i+fNData/2][j+fNData/2].name = "CORR";
	      
	      fSys[fNData/2-1-i][fNData/2-1-j].add  = syscor[i][j]/sqrt(2.0)/2.;
	      fSys[fNData/2-1-i][fNData/2-1-j].mult = fSys[fNData/2-1-i][fNData/2-1-j].add*100/fData[fNData/2-1-i];
	      fSys[fNData/2-1-i][fNData/2-1-j].type = MULT;
	      fSys[fNData/2-1-i][fNData/2-1-j].name = "CORR";

	      fSys[fNData/2-1-i][j+fNData/2].add  = syscor[i][j]/sqrt(2.0)/2.;
	      fSys[fNData/2-1-i][j+fNData/2].mult = fSys[fNData/2-1-i][j+fNData/2].add*100/fData[fNData/2-1-i];
	      fSys[fNData/2-1-i][j+fNData/2].type = MULT;
	      fSys[fNData/2-1-i][j+fNData/2].name = "CORR";
	    }
	}  
      */

      //Real systematics
      for(int j=0; j<fNSys; j++)
	{
	  double sys1, sys2, right, left;
	  double stmp, dtmp;
	  string sysdescr;
	  
	  ostringstream id;
	  id << j;
	  sysdescr = "CORR";
	  
	  lstream >> sys1 >> sys2; 
	  if(sys1<0) {right=sys2; left=sys1;}
	  else {right=sys1; left=sys2;}

	  //convert to relative percentage values
	  right = right/datum*100;  
	  left = left/datum*100;
	  symmetriseErrors(right,left,&stmp,&dtmp);
	  
	  fSys[i+fNData/2][j].type = MULT;
	  fSys[i+fNData/2][j].name = sysdescr;
	  fSys[i+fNData/2][j].mult = stmp;
	  fSys[i+fNData/2][j].add  = fSys[i+fNData/2][j].mult*fData[i+fNData/2]/100;

	  fSys[fNData/2-1-i][j].type = MULT;
	  fSys[fNData/2-1-i][j].name = sysdescr;
	  fSys[fNData/2-1-i][j].mult = stmp;
	  fSys[fNData/2-1-i][j].add  = fSys[fNData/2-1-i][j].mult*fData[fNData/2-1-i]/100;

	  shift += dtmp;

	}
          
       fData[i+fNData/2]*=(1.0 + shift*0.01);   //Shift from asymmetric errors
       fData[fNData/2-1-i]*=(1.0 + shift*0.01); //Shift from asymmetric errors

    }  

  for(int i = 0; i < fNData/2; i++) 
    {
      delete[] syscor[i];
      delete[] covmat[i];
    }
  delete[] syscor;
  delete[] covmat;
  
  f1.close();
  f2.close();
}

//=================================================================

//5) Distribution differential in top quark pair invariant mass
void  ATLASTOPDIFF8TEVTTMNORMFilter::ReadData()
{
  // Opening files
  fstream f1, f2;

  //Central values, statistical and systematic uncertainties  
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/" << fSetName << "/" << fSetName << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

//Statistical covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/" << fSetName << "/ttbar_mass_covariance_rel.dat";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read statistical covariance matrix
  string line;
  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
  {
    covmat[i] = new double[fNData];
    getline(f2,line);
    istringstream lstream(line);
    for(int j=0; j<fNData; j++)
    {
      lstream >> covmat[i][j];
    }
  }

  //Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
    {
      cerr << " in " << fSetName << endl;
      exit(-1);
    }

  //Starting filter
  for(int i=0; i<5; i++)
    {
      getline(f1,line);
    }
  
  for(int i=0; i<fNData;i++)
    {
      double pt_top, ddum;
      getline(f1,line);
      istringstream lstream(line);
      lstream >> pt_top >> ddum >> ddum; 

      fKin1[i] = pt_top;       //P_T^(top)
      fKin2[i] = Mt*Mt;       
      fKin3[i] = 8000;         //sqrt(s)

      lstream >> fData[i];     //differential distribution
      lstream >> fStat[i];     //its statistical uncertainty
      fStat[i]=0.;
      lstream >> ddum >> ddum >> ddum >> ddum >> ddum >> ddum;
      
      double shift = 0.;

      //Artificial systematics
      for(int i=0; i<fNData; i++)
	{
	  for(int j=0; j<fNData; j++)
	    {
	      fSys[i][j].add  = syscor[i][j];
	      fSys[i][j].mult = fSys[i][j].add*100/fData[i];
	      fSys[i][j].type = MULT;
	      fSys[i][j].name = "CORR";
	    }
	}  

      //Real systematics
      for(int j=fNData; j<fNSys; j++)
	{
	  double sys1, sys2, right, left;
	  double stmp, dtmp;
	  string sysdescr;
	  
	  ostringstream id;
	  id << j;
	  sysdescr = "CORR";
	  
	  lstream >> sys1 >> sys2;
	  if(sys1<0) {right=sys2; left=sys1;}
	  else {right=sys1; left=sys2;}
	  //convert to relative percentage values
	  right = right/fData[i]*100;  
	  left = left/fData[i]*100;
	  symmetriseErrors(right,left,&stmp,&dtmp);
	  
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = sysdescr;
	  fSys[i][j].mult = stmp;
	  fSys[i][j].add  = fSys[i][j].mult*fData[i]/100;

	  shift += dtmp;

	}
          
      fData[i]*=(1.0 + shift*0.01); //Shift from asymmetric errors

    }  

  for(int i = 0; i < fNData; i++) 
    {
      delete[] syscor[i];
      delete[] covmat[i];
    }
  delete[] syscor;
  delete[] covmat;
  
  f1.close();
  f2.close();

}

/*
========================================================================
*/

//B - UNNORMALISED distributions

//1) Distribution differential in top quark transverse momentum
void  ATLASTOPDIFF8TEVTPTFilter::ReadData()
{
  // Opening files
  fstream f1, f2;

  //Central values, statistical and systematic uncertainties  
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/" << fSetName << "/" << fSetName << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Statistical covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/" << fSetName << "/top_pt_covariance_abs.dat";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read statistical covariance matrix
  string line;
  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
  {
    covmat[i] = new double[fNData];
    getline(f2,line);
    istringstream lstream(line);
    for(int j=0; j<fNData; j++)
    {
      lstream >> covmat[i][j];
    }
  }

  //Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
    {
      cerr << " in " << fSetName << endl;
      exit(-1);
    }

  //Starting filter
  for(int i=0; i<5; i++)
    {
      getline(f1,line);
    }
  
  for(int i=0; i<fNData;i++)
    {
      double pt_top, ddum;
      getline(f1,line);
      istringstream lstream(line);
      lstream >> pt_top >> ddum >> ddum; 
      
      fKin1[i] = pt_top;       //P_T^(top)
      fKin2[i] = Mt*Mt;       
      fKin3[i] = 8000;         //sqrt(s)
      
      lstream >> fData[i];     //differential distribution
      lstream >> fStat[i];     //its statistical uncertainty
      fStat[i] = 0.;
      lstream >> ddum >> ddum >> ddum >> ddum >> ddum >> ddum;
      
      double shift = 0.;
      
      //Artificial systematics
      for(int l=0; l<fNData; l++)
	{
	  for(int j=0; j<fNData; j++)
	    {
	      fSys[l][j].add  = syscor[l][j];
	      fSys[l][j].mult = fSys[l][j].add*100/fData[l];
	      fSys[l][j].type = MULT;
	      fSys[l][j].name = "CORR";
	    }
	}        

      //Real systematics
      for(int j=fNData; j<fNSys-1; j++)
	{
	  double sys1, sys2, right, left;
	  double stmp, dtmp;
	  string sysdescr;
	  
	  ostringstream id;
	  id << j;
	  sysdescr = "CORR";
	  
	  lstream >> sys1 >> sys2;
	  if(sys1<0) {right=sys2; left=sys1;}
	  else {right=sys1; left=sys2;}
	  //convert to relative percentage values
	  right = right/fData[i]*100;  
	  left = left/fData[i]*100;
	  symmetriseErrors(right,left,&stmp,&dtmp);
	  
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = sysdescr;
	  fSys[i][j].mult = stmp;
	  fSys[i][j].add  = fSys[i][j].mult*fData[i]/100;
	  
	  shift += dtmp;
	  
	}

      //Overall luminosity uncertainty
      fSys[i][64].type = MULT;
      fSys[i][64].name = "ATLASLUMI12";
      fSys[i][64].mult=2.8;
      fSys[i][64].add=fSys[i][64].mult*fData[i]/100;
      
      fData[i]*=(1.0 + shift*0.01); //Shift from asymmetric errors
      
    }  

  for(int i = 0; i < fNData; i++) 
    {
      delete[] syscor[i];
      delete[] covmat[i];
    }
  delete[] syscor;
  delete[] covmat;
  
  f1.close();
  f2.close();
}

//==================================================================

//2) Distribution differential in top quark pair transverse momentum
void  ATLASTOPDIFF8TEVTTPTFilter::ReadData()
{
  // Opening files
  fstream f1;

  //Central values, statistical and systematic uncertainties  
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/" << fSetName << "/" << fSetName << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Starting filter
  string line;
  for(int i=0; i<5; i++)
    {
      getline(f1,line);
    }
  
  for(int i=0; i<fNData;i++)
    {
      double pt_top, ddum;
      getline(f1,line);
      istringstream lstream(line);
      lstream >> pt_top >> ddum >> ddum; 

      fKin1[i] = pt_top;       //P_T^(top)
      fKin2[i] = Mt*Mt;       
      fKin3[i] = 8000;         //sqrt(s)

      lstream >> fData[i];     //differential distribution
      lstream >> fStat[i];     //its statistical uncertainty
      lstream >> ddum >> ddum >> ddum >> ddum >> ddum >> ddum;
      
      double shift = 0.;

      for(int j=0; j<fNSys-1; j++)
	{
	  double sys1, sys2, right, left;
	  double stmp, dtmp;
	  string sysdescr;
	  
	  ostringstream id;
	  id << j;
	  sysdescr = "CORR";
	  
	  lstream >> sys1 >> sys2;
	  if(sys1<0) {right=sys2; left=sys1;}
	  else {right=sys1; left=sys2;}

	  //convert to relative percentage values
	  right = right/fData[i]*100;  
	  left  = left/fData[i]*100;
	  symmetriseErrors(right,left,&stmp,&dtmp);
	  
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = sysdescr;
	  fSys[i][j].mult = stmp;
	  fSys[i][j].add  = fSys[i][j].mult*fData[i]/100;

	  shift += dtmp;

	}

      //overall luminosity uncertainty
      fSys[i][56].type = MULT;
      fSys[i][56].name = "ATLASLUMI12";
      fSys[i][56].mult=2.8;
      fSys[i][56].add=fSys[i][56].mult*fData[i]/100;
         
      fData[i]*=(1.0 + shift*0.01); //Shift from asymmetric errors
 
    }  
  
  f1.close();
}

//==============================================================

//3) Distribution differential in top quark rapidity
void  ATLASTOPDIFF8TEVTRAPFilter::ReadData()
{
  // Opening files
  fstream f1, f2;

  //Central values, statistical and systematic uncertainties  
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/" << fSetName << "/" << fSetName << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Statistical covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/" << fSetName << "/top_absy_covariance_abs.dat";
  f2.open(covfile.str().c_str(), ios::in);
  
  if (f2.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read statistical covariance matrix
  string line;
  double** covmat = new double*[fNData/2];
  for(int i=0; i<fNData/2; i++)
    {
      covmat[i] = new double[fNData/2];
      getline(f2,line);
      istringstream lstream(line);
      for(int j=0; j<fNData/2; j++)
	{
	  lstream >> covmat[i][j];
	}
    }

  //Generate artificial systematics
  double** syscor = new double*[fNData/2];
  for(int i = 0; i < fNData/2; i++)
    syscor[i] = new double[fNData/2];
  
  if(!genArtSys(fNData/2,covmat,syscor))
    {
      cerr << " in " << fSetName << endl;
      exit(-1);
    } 

  //Starting filter
  for(int i=0; i<5; i++)
    {
      getline(f1,line);
    }
  
  for(int i=0; i<fNData/2;i++)
    {
      double pt_top, ddum, datum, stat;
      getline(f1,line);
      istringstream lstream(line);
      lstream >> pt_top >> ddum >> ddum; 

      fKin1[i+fNData/2]   = pt_top;          //P_T^(top)
      fKin1[fNData/2-1-i] = -1.0*pt_top;
      fKin2[i+fNData/2]   = Mt*Mt;       
      fKin2[fNData/2-1-i] = Mt*Mt;
      fKin3[i+fNData/2]   = 8000;            //sqrt(s)
      fKin3[fNData/2-1-i] = 8000;        

      lstream >> datum;
      fData[i+fNData/2]   = datum/2.0;       //differential distribution 
      fData[fNData/2-1-i] = datum/2.0;

      lstream >> stat;
      fStat[i+fNData/2]   = stat/sqrt(2.0);  //statistical uncertainty
      fStat[fNData/2-1-i] = stat/sqrt(2.0);
      //fStat[i+fNData/2]   = 0.;            //statistical uncertainty
      //fStat[fNData/2-1-i] = 0.;

      lstream >> ddum >> ddum >> ddum >> ddum >> ddum >> ddum;
   
      double shift = 0.;

      /*
      //Artificial systematics
      for(int l=0; l<fNData/2; l++)
	{
	  for(int j=0; j<fNData/2; j++)
	    {
	     
	      fSys[l+fNData/2][j].add  = syscor[l][j]/sqrt(2.0);
	      fSys[l+fNData/2][j].mult = fSys[l+fNData/2][j].add*100/fData[l+fNData/2];
	      fSys[l+fNData/2][j].type = MULT;
	      fSys[l+fNData/2][j].name = "CORR";

	      fSys[fNData/2-1-l][j].add  = syscor[l][j]/sqrt(2.0)/2.;
	      fSys[fNData/2-1-l][j].mult = fSys[fNData/2-1-l][j].add*100/fData[fNData/2-1-l];
	      fSys[fNData/2-1-l][j].type = MULT;
	      fSys[fNData/2-1-l][j].name = "CORR";

	      
	      const double eps=0.0000001;

	      fSys[i+fNData/2][fNData/2-1-j].add  = syscor[i][j]/sqrt(2.0)/(2.+eps);
	      fSys[i+fNData/2][fNData/2-1-j].mult = fSys[i+fNData/2][fNData/2-1-j].add*100/fData[i+fNData/2];
	      fSys[i+fNData/2][fNData/2-1-j].type = MULT;
	      fSys[i+fNData/2][fNData/2-1-j].name = "CORR";

	      fSys[i+fNData/2][j+fNData/2].add  = syscor[i][j]/sqrt(2.0)/(2.+eps);
	      fSys[i+fNData/2][j+fNData/2].mult = fSys[i+fNData/2][j+fNData/2].add*100/fData[i+fNData/2];
	      fSys[i+fNData/2][j+fNData/2].type = MULT;
	      fSys[i+fNData/2][j+fNData/2].name = "CORR";
	      
	      fSys[fNData/2-1-i][fNData/2-1-j].add  = syscor[i][j]/sqrt(2.0)/2.;
	      fSys[fNData/2-1-i][fNData/2-1-j].mult = fSys[fNData/2-1-i][fNData/2-1-j].add*100/fData[fNData/2-1-i];
	      fSys[fNData/2-1-i][fNData/2-1-j].type = MULT;
	      fSys[fNData/2-1-i][fNData/2-1-j].name = "CORR";

	      fSys[fNData/2-1-i][j+fNData/2].add  = syscor[i][j]/sqrt(2.0)/2.;
	      fSys[fNData/2-1-i][j+fNData/2].mult = fSys[fNData/2-1-i][j+fNData/2].add*100/fData[fNData/2-1-i];
	      fSys[fNData/2-1-i][j+fNData/2].type = MULT;
	      fSys[fNData/2-1-i][j+fNData/2].name = "CORR";
	      	      

	    }
	}  
      */

      //Real systematics
      for(int j=0; j<fNSys-1; j++)
	{
	  double sys1, sys2, right, left;
	  double stmp, dtmp;
	  string sysdescr;
	  
	  ostringstream id;
	  id << j;    
	  sysdescr = "CORR";
	  
	  lstream >> sys1 >> sys2; 

	  if(sys1<0) {right=sys2; left=sys1;}
	  else {right=sys1; left=sys2;}

	  //convert to relative percentage values
	  right = right/datum*100;  
	  left = left/datum*100;
	  symmetriseErrors(right,left,&stmp,&dtmp);
	  
	  fSys[i+fNData/2][j].type = MULT;
	  fSys[i+fNData/2][j].name = sysdescr;
	  fSys[i+fNData/2][j].mult = stmp;
	  fSys[i+fNData/2][j].add  = fSys[i+fNData/2][j].mult*fData[i+fNData/2]/100;

	  fSys[fNData/2-1-i][j].type = MULT;
	  fSys[fNData/2-1-i][j].name = sysdescr;
	  fSys[fNData/2-1-i][j].mult = stmp;
	  fSys[fNData/2-1-i][j].add  = fSys[fNData/2-1-i][j].mult*fData[fNData/2-1-i]/100;

	  shift += dtmp;

	}

      //overall luminosity uncertainty
      fSys[i+fNData/2][56].type = MULT;
      fSys[i+fNData/2][56].name ="ATLASLUMI12";
      fSys[i+fNData/2][56].mult=2.8;
      fSys[i+fNData/2][56].add=fSys[i+fNData/2][56].mult*fData[i+fNData/2]/100;
      fSys[fNData/2-1-i][56].type = MULT;
      fSys[fNData/2-1-i][56].name = "ATLASLUMI12";
      fSys[fNData/2-1-i][56].mult=2.8;
      fSys[fNData/2-1-i][56].add=fSys[fNData/2-1-i][56].mult*fData[fNData/2-1-i]/100;
      
      fData[i+fNData/2]*=(1.0 + shift*0.01); //Shift from asymmetric errors
      fData[fNData/2-1-i]*=(1.0 + shift*0.01); //Shift from asymmetric errors

    }  

  /*
  //Artificial systematics
  for(int i=0; i<fNData/2; i++)
    {
      for(int j=0; j<fNData/2; j++)
	{
	  	  
	  fSys[i+fNData/2][fNSys-fNData/2+j].add  = syscor[i][j]/sqrt(2.0);
	  fSys[i+fNData/2][fNSys-fNData/2+j].mult = fSys[i+fNData/2][fNSys-fNData/2+j].add*100/fData[i+fNData/2];
	  fSys[i+fNData/2][fNSys-fNData/2+j].type = MULT;
	  fSys[i+fNData/2][fNSys-fNData/2+j].name = "CORR";
	  
	  fSys[fNData/2-1-i][fNSys-fNData/2+j].add  = syscor[i][j]/sqrt(2.0);
	  fSys[fNData/2-1-i][fNSys-fNData/2+j].mult = fSys[fNData/2-1-i][fNSys-fNData/2+j].add*100/fData[fNData/2-1-i];
	  fSys[fNData/2-1-i][fNSys-fNData/2+j].type = MULT;
	  fSys[fNData/2-1-i][fNSys-fNData/2+j].name = "CORR";
	  	  
	}
    } 
*/
  for(int i = 0; i < fNData/2; i++) 
    {
      delete[] syscor[i];
      delete[] covmat[i];
    }
  delete[] syscor;
  delete[] covmat;
  
  f1.close();
  f2.close();

}


//=================================================================

//4) Distribution differential in top quark pair rapidity
void  ATLASTOPDIFF8TEVTTRAPFilter::ReadData()
{
  // Opening files
  fstream f1, f2;

  //Central values, statistical and systematic uncertainties  
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/" << fSetName << "/" << fSetName << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Statistical covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/" << fSetName << "/ttbar_absy_covariance_abs.dat";
  f2.open(covfile.str().c_str(), ios::in);
  
  if (f2.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  //Read statistical covariance matrix
  string line;
  double** covmat = new double*[fNData/2];
  for(int i=0; i<fNData/2; i++)
    {
      covmat[i] = new double[fNData/2];
      getline(f2,line);
      istringstream lstream(line);
      for(int j=0; j<fNData/2; j++)
	{
	  lstream >> covmat[i][j];
	}
    }
  
  //Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData/2; i++)
    syscor[i] = new double[fNData/2];
  
  if(!genArtSys(fNData/2,covmat,syscor))
    {
      cerr << " in " << fSetName << endl;
      exit(-1);
    }  

  //Starting filter
  for(int i=0; i<5; i++)
    {
      getline(f1,line);
    }
  
  for(int i=0; i<fNData/2;i++)
    {
      double pt_top, ddum, datum, stat;
      getline(f1,line);
      istringstream lstream(line);
      lstream >> pt_top >> ddum >> ddum; 

      fKin1[i+fNData/2]   = pt_top;          //P_T^(top)
      fKin1[fNData/2-1-i] = -1.0*pt_top;
      fKin2[i+fNData/2]   = Mt*Mt;       
      fKin2[fNData/2-1-i] = Mt*Mt;
      fKin3[i+fNData/2]   = 8000;            //sqrt(s)
      fKin3[fNData/2-1-i] = 8000;        

      lstream >> datum;
      fData[i+fNData/2]   = datum/2.0;       //differential distribution 
      fData[fNData/2-1-i] = datum/2.0;

      lstream >> stat;
      fStat[i+fNData/2]   = stat/sqrt(2.0);  //statistical uncertainty
      fStat[fNData/2-1-i] = stat/sqrt(2.0);
      //fStat[i+fNData/2] = 0.;
      //fStat[fNData/2-1-i] = 0.;

      lstream >> ddum >> ddum >> ddum >> ddum >> ddum >> ddum;
      
      double shift = 0.;

      /*
      //Artificial systematics
      for(int i=0; i<fNData/2; i++)
	{
	  for(int j=0; j<fNData/2; j++)
	    {
	      const double eps=0.0000001;

	      fSys[i+fNData/2][fNData/2-1-j].add  = syscor[i][j]/sqrt(2.0)/(2.+eps);
	      fSys[i+fNData/2][fNData/2-1-j].mult = fSys[i+fNData/2][fNData/2-1-j].add*100/fData[i+fNData/2];
	      fSys[i+fNData/2][fNData/2-1-j].type = MULT;
	      fSys[i+fNData/2][fNData/2-1-j].name = "CORR";

	      fSys[i+fNData/2][j+fNData/2].add  = syscor[i][j]/sqrt(2.0)/(2.+eps);
	      fSys[i+fNData/2][j+fNData/2].mult = fSys[i+fNData/2][j+fNData/2].add*100/fData[i+fNData/2];
	      fSys[i+fNData/2][j+fNData/2].type = MULT;
	      fSys[i+fNData/2][j+fNData/2].name = "CORR";
	      
	      fSys[fNData/2-1-i][fNData/2-1-j].add  = syscor[i][j]/sqrt(2.0)/2.;
	      fSys[fNData/2-1-i][fNData/2-1-j].mult = fSys[fNData/2-1-i][fNData/2-1-j].add*100/fData[fNData/2-1-i];
	      fSys[fNData/2-1-i][fNData/2-1-j].type = MULT;
	      fSys[fNData/2-1-i][fNData/2-1-j].name = "CORR";

	      fSys[fNData/2-1-i][j+fNData/2].add  = syscor[i][j]/sqrt(2.0)/2.;
	      fSys[fNData/2-1-i][j+fNData/2].mult = fSys[fNData/2-1-i][j+fNData/2].add*100/fData[fNData/2-1-i];
	      fSys[fNData/2-1-i][j+fNData/2].type = MULT;
	      fSys[fNData/2-1-i][j+fNData/2].name = "CORR";
	    }
	}  
      */

      //Real systematics
      for(int j=0; j<fNSys-1; j++)
	{
	  double sys1, sys2, right, left;
	  double stmp, dtmp;
	  string sysdescr;
	  
	  ostringstream id;
	  id << j;
	  sysdescr = "CORR";
	  
	  lstream >> sys1 >> sys2; 
	  if(sys1<0) {right=sys2; left=sys1;}
	  else {right=sys1; left=sys2;}

	  //convert to relative percentage values
	  right = right/datum*100;  
	  left = left/datum*100;
	  symmetriseErrors(right,left,&stmp,&dtmp);
	  
	  fSys[i+fNData/2][j].type = MULT;
	  fSys[i+fNData/2][j].name = sysdescr;
	  fSys[i+fNData/2][j].mult = stmp;
	  fSys[i+fNData/2][j].add  = fSys[i+fNData/2][j].mult*fData[i+fNData/2]/100;

	  fSys[fNData/2-1-i][j].type = MULT;
	  fSys[fNData/2-1-i][j].name = sysdescr;
	  fSys[fNData/2-1-i][j].mult = stmp;
	  fSys[fNData/2-1-i][j].add  = fSys[fNData/2-1-i][j].mult*fData[fNData/2-1-i]/100;

	  shift += dtmp;

	}
          
      //overall luminosity uncertainty
      fSys[i+fNData/2][56].type = MULT;
      fSys[i+fNData/2][56].name = "ATLASLUMI12";
      fSys[i+fNData/2][56].mult=2.8;
      fSys[i+fNData/2][56].add=fSys[i+fNData/2][56].mult*fData[i+fNData/2]/100;
      fSys[fNData/2-1-i][56].type = MULT;
      fSys[fNData/2-1-i][56].name = "ATLASLUMI12";
      fSys[fNData/2-1-i][56].mult=2.8;
      fSys[fNData/2-1-i][56].add=fSys[fNData/2-1-i][56].mult*fData[fNData/2-1-i]/100;

      fData[i+fNData/2]*=(1.0 + shift*0.01); //Shift from asymmetric errors
      fData[fNData/2-1-i]*=(1.0 + shift*0.01); //Shift from asymmetric errors
    
    }  

  for(int i = 0; i < fNData/2; i++) 
    {
      delete[] syscor[i];
      delete[] covmat[i];
    }
  delete[] syscor;
  delete[] covmat;
  
  f1.close();
  f2.close();
}

//=================================================================

//5) Distribution differential in top quark pair invariant mass
void  ATLASTOPDIFF8TEVTTMFilter::ReadData()
{
  // Opening files
  fstream f1, f2;

  //Central values, statistical and systematic uncertainties  
  stringstream datafile("");
  datafile << dataPath() 
	   << "rawdata/" << fSetName << "/" << fSetName << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Statistical covariance matrix
  stringstream covfile("");
  covfile << dataPath()
	  << "rawdata/" << fSetName << "/ttbar_mass_covariance_abs.dat";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Read statistical covariance matrix
  string line;
  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
  {
    covmat[i] = new double[fNData];
    getline(f2,line);
    istringstream lstream(line);
    for(int j=0; j<fNData; j++)
    {
      lstream >> covmat[i][j];
    }
  }

  //Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
    {
      cerr << " in " << fSetName << endl;
      exit(-1);
    }

  //Starting filter
  for(int i=0; i<5; i++)
    {
      getline(f1,line);
    }
  
  for(int i=0; i<fNData;i++)
    {
      double pt_top, ddum;
      getline(f1,line);
      istringstream lstream(line);
      lstream >> pt_top >> ddum >> ddum; 

      fKin1[i] = pt_top;       //P_T^(top)
      fKin2[i] = Mt*Mt;       
      fKin3[i] = 8000;         //sqrt(s)

      lstream >> fData[i];     //differential distribution
      lstream >> fStat[i];     //its statistical uncertainty
      fStat[i] = 0.;
      lstream >> ddum >> ddum >> ddum >> ddum >> ddum >> ddum;
      
      double shift = 0.;

      //Artificial systematics
      for(int i=0; i<fNData; i++)
	{
	  for(int j=0; j<fNData; j++)
	    {
	      fSys[i][j].add  = syscor[i][j];
	      fSys[i][j].mult = fSys[i][j].add*100/fData[i];
	      fSys[i][j].type = MULT;
	      fSys[i][j].name = "CORR";
	    }
	}  
      
      //Real systematics
      for(int j=fNData; j<fNSys-1; j++)
	{
	  double sys1, sys2, right, left;
	  double stmp, dtmp;
	  string sysdescr;
	  
	  ostringstream id;
	  id << j;
	  sysdescr = "CORR";
	  
	  lstream >> sys1 >> sys2;
	  if(sys1<0) {right=sys2; left=sys1;}
	  else {right=sys1; left=sys2;}

	  //convert to relative percentage values
	  right = right/fData[i]*100;  
	  left  = left/fData[i]*100;
	  symmetriseErrors(right,left,&stmp,&dtmp);
	  
	  fSys[i][j].type = MULT;
	  fSys[i][j].name = sysdescr;
	  fSys[i][j].mult = stmp;
	  fSys[i][j].add  = fSys[i][j].mult*fData[i]/100;

	  shift += dtmp;

	}

      //overall luminosity uncertainty
      fSys[i][63].type = MULT;
      fSys[i][63].name = "ATLASLUMI12";
      fSys[i][63].mult=2.8;
      fSys[i][63].add=fSys[i][63].mult*fData[i]/100;
         
      fData[i]*=(1.0 + shift*0.01); //Shift from asymmetric errors
 
    }  

  for(int i = 0; i < fNData; i++) 
    {
      delete[] syscor[i];
      delete[] covmat[i];
    }
  delete[] syscor;
  delete[] covmat;
  
  f1.close();
  f2.close();
}
