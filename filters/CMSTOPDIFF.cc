/**
 
This file describes the buildmaster implementation of the
differential distributions for top quark pair production at 8 TeV
from ATLAS and CMS

ATLAS => http://arxiv.org/abs/1511.04716
CMS => http://arxiv.org/abs/1505.04480

There are two different ways to include this data, with
normalized or with absolute distributions

By modifying this file, one can select the normalized or the absolute
distributions

*/

#include "CMSTOPDIFF.h"

// Raw data available from HepData
// Data is available separately for the lepton+jets and
// for the dilepton channel
// http://hepdata.cedar.ac.uk/view/ins1370682

//-------------------------------------------------------------------------------

// CMS top quark differential distributions
// 8 TeV, top quark pt distribution
void  CMSTOPDIFF8TEVTPTFilter::ReadData()
{
  // Opening files
  fstream f1;

  //
  // The raw data for (1/sigma) dsigma/dpt is taken from
  // Table 15 ( Page 41 of preprint, table A.6, 1st variable in table, pt(top).
  //
  //   A	A	A	pt_top	pt_top	pttop  	y	dy+	dy-	dy+	dy-	
  //
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
	   << fSetName << "/CMSTOPDIFF8TEVTPT.dat";
  f1.open(datafile.str().c_str(), ios::in);
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
  }
  string line;
  for(int i=0;i<10;i++) getline(f1,line);
  for (int i = 0; i < fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
            
      int idum=0;
      lstream >> idum; 
      if(idum!=(i+1)){
	cout<<"Error in CMSTOPDIFF.cc"<<endl;
	cout<<"idum = "<<idum<<endl;
	exit(-10);
      }

      double adum=0;
      for(int j=0;j<5;j++)lstream >> adum;
      double pt_top=0.0;
      lstream >> pt_top;

      fKin1[i] = pt_top; // <pt_top> as characteristic kin variable
      fKin2[i] = Mt;
      fKin3[i] = 8000;
      
      lstream >> fData[i]; // Central values
      double stat=0;
      double sys=0;
      lstream>> stat >>adum; // Statistical error
      lstream >> sys >>adum; // Systematic error (absolute value)

      fStat[i]=pow(pow(stat,2.0) + pow(sys,2.0), 0.5); // 
      // For the time being add in quadrature statistical and systematic uncertainties
      
      //cout<<i<<" "<<fKin2[i]<<" "<<fData[i]<<" "<<1e2*fStat[i]/fData[i]<<endl;
      
    }
  
  f1.close();
}

///////////////////////////////////////////////////////////

//-------------------------------------------------------------------------------

// CMS top quark differential distributions
// 8 TeV, top quark rapidity distribution
void  CMSTOPDIFF8TEVTRAPFilter::ReadData()
{
  // Opening files
  fstream f1;

  //
  // The raw data for (1/sigma) dsigma/dy(top) is taken from
  // Table 17 ( Page 42 of preprint, table A.7, 1st variable in table, y(top).
  //
  //   A	A	A	y_top	y_top	y_top  	y	dy+	dy-	dy+	dy-	
  //
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
	   << fSetName << "/CMSTOPDIFF8TEVTRAP.dat";
  f1.open(datafile.str().c_str(), ios::in);
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
  }
  string line;
  for(int i=0;i<10;i++) getline(f1,line);
  for (int i = 0; i < fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
            
      int idum=0;
      lstream >> idum; 
      if(idum!=(i+1)){
	cout<<"Error in CMSTOPDIFF.cc"<<endl;
	cout<<"idum = "<<idum<<endl;
	exit(-10);
      }

      double adum=0;
      for(int j=0;j<5;j++)lstream >> adum;
      double y_top=0.0;
      lstream >> y_top;

      fKin1[i] = y_top; // <y_top> as characteristic kin variable
      fKin2[i] = Mt;
      fKin3[i] = 8000;   

      lstream >> fData[i]; // Central values
      double stat=0;
      double sys=0;
      lstream>> stat >>adum; // Statistical error
      lstream >> sys >>adum; // Systematic error (absolute value)

      fStat[i]=pow(pow(stat,2.0) + pow(sys,2.0), 0.5); // 
      // For the time being add in quadrature statistical and systematic uncertainties
      
      // cout<<i<<" "<<fKin2[i]<<" "<<fData[i]<<" "<<1e2*fStat[i]/fData[i]<<endl;
      
    }
  
  f1.close();
}

///////////////////////////////////////////////////////////

// CMS top quark differential distributions
// 8 TeV, top quark pair invariant mass distribution
// These are differential distributions for the lepton+jets channel
// Both absolute distributions and dilepton final state distributions also available
void  CMSTOPDIFF8TEVTTMFilter::ReadData()
{
  // Opening files
  fstream f1;

  //
  // The raw data for (1/sigma) dsigma/dm(tt) is taken from
  // Table 23 ( Page 44 of preprint, table A.9, 3rd variable in table, m(ttbar)
  //
  //   A	A	A	m_tt	m_tt	m_tt  	y	dy+	dy-	dy+	dy-	
  //
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
	   << fSetName << "/CMSTOPDIFF8TEVTTM.dat";
  f1.open(datafile.str().c_str(), ios::in);
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
  }
  string line;
  for(int i=0;i<10;i++) getline(f1,line);
  for (int i = 0; i < fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
            
      int idum=0;
      lstream >> idum; 
      if(idum!=(i+1)){
	cout<<"Error in CMSTOPDIFF.cc"<<endl;
	cout<<"idum = "<<idum<<endl;
	exit(-10);
      }

      double adum=0;
      for(int j=0;j<5;j++)lstream >> adum;
      double m_tt=0.0;
      lstream >> m_tt;

      fKin1[i] = m_tt; // <m_tt> as characteristic kin variable
      fKin2[i] = Mt;
      fKin3[i] = 8000;

      lstream >> fData[i]; // Central values
      double stat=0;
      double sys=0;
      lstream>> stat >>adum; // Statistical error
      lstream >> sys >>adum; // Systematic error (absolute value)

      fStat[i]=pow(pow(stat,2.0) + pow(sys,2.0), 0.5); // 
      // For the time being add in quadrature statistical and systematic uncertainties
      
      // cout<<i<<" "<<fKin2[i]<<" "<<fData[i]<<" "<<1e2*fStat[i]/fData[i]<<endl;
      
    }
  
  f1.close();
}

///////////////////////////////////////////////////////////

// CMS top quark differential distributions
// 8 TeV, top quark pair rapidity distribution
void  CMSTOPDIFF8TEVTTRAPFilter::ReadData()
{
  // Opening files
  fstream f1;

  //
  // The raw data for (1/sigma) dsigma/dy(top) is taken from
  // Table 17 ( Page 42 of preprint, table A.7, 1st variable in table, y(top).
  //
  //   A	A	A	y_top	y_top	y_top  	y	dy+	dy-	dy+	dy-	
  //
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
	   << fSetName << "/CMSTOPDIFF8TEVTTRAP.dat";
  f1.open(datafile.str().c_str(), ios::in);
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
  }
  string line;
  for(int i=0;i<10;i++) getline(f1,line);
  for (int i = 0; i < fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
            
      int idum=0;
      lstream >> idum; 
      if(idum!=(i+1)){
	cout<<"Error in CMSTOPDIFF.cc"<<endl;
	cout<<"idum = "<<idum<<endl;
	exit(-10);
      }

      double adum=0;
      for(int j=0;j<5;j++)lstream >> adum;
      double y_tt=0.0;
      lstream >> y_tt;

      fKin1[i] = y_tt; // <y_tt> as characteristic kin variable
      fKin2[i] = Mt;
      fKin3[i] = 8000;   

      lstream >> fData[i]; // Central values
      double stat=0;
      double sys=0;
      lstream>> stat >>adum; // Statistical error
      lstream >> sys >>adum; // Systematic error (absolute value)

      fStat[i]=pow(pow(stat,2.0) + pow(sys,2.0), 0.5); // 
      // For the time being add in quadrature statistical and systematic uncertainties
      
      // cout<<i<<" "<<fKin2[i]<<" "<<fData[i]<<" "<<1e2*fStat[i]/fData[i]<<endl;
      
    }
  
  f1.close();
}

///////////////////////////////////////////////////////////

// CMS top quark differential distributions
// 8 TeV, top quark pair transverse momentum
void  CMSTOPDIFF8TEVTTPTFilter::ReadData()
{
  // Opening files
  fstream f1;

  //
  // The raw data for (1/sigma) dsigma/dy(top) is taken from
  // Table 17 ( Page 42 of preprint, table A.7, 1st variable in table, y(top).
  //
  //   A	A	A	y_top	y_top	y_top  	y	dy+	dy-	dy+	dy-	
  //
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
	   << fSetName << "/CMSTOPDIFF8TEVTTPT.dat";
  f1.open(datafile.str().c_str(), ios::in);
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
  }
  string line;
  for(int i=0;i<10;i++) getline(f1,line);
  for (int i = 0; i < fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
            
      int idum=0;
      lstream >> idum; 
      if(idum!=(i+1)){
	cout<<"Error in CMSTOPDIFF.cc"<<endl;
	cout<<"idum = "<<idum<<endl;
	exit(-10);
      }

      double adum=0;
      for(int j=0;j<5;j++)lstream >> adum;
      double pt_tt=0.0;
      lstream >> pt_tt;

      fKin1[i] = pt_tt; // <pt_tt> as characteristic kin variable
      fKin2[i] = Mt;
      fKin3[i] = 8000;    

      lstream >> fData[i]; // Central values
      double stat=0;
      double sys=0;
      lstream>> stat >>adum; // Statistical error
      lstream >> sys >>adum; // Systematic error (absolute value)

      fStat[i]=pow(pow(stat,2.0) + pow(sys,2.0), 0.5); // 
      // For the time being add in quadrature statistical and systematic uncertainties
      
      // cout<<i<<" "<<fKin2[i]<<" "<<fData[i]<<" "<<1e2*fStat[i]/fData[i]<<endl;
      
    }
  
  f1.close();
}

///////////////////////////////////////////////////////////

