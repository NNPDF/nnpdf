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

#include "ATLASTOPDIFF.h"

// Raw data available here
// https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/TOPQ-2015-06/#auxstuff

//-------------------------------------------------------------------------------

// ATLAS top quark differential distributions
// 8 TeV, top quark pt distribution
void  ATLASTOPDIFF8TEVTPTFilter::ReadData()
{
  // Opening files
  fstream f1;

  //
  // The raw data for (1/sigma) dsigma/dpt is taken from
  // https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/TOPQ-2015-06/tabaux_025.pdf
  // Format
  // pt_min  pt_max   (1/sigma)dsigma/dpt   Total_Unc(%)
  //
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
	   << fSetName << "/ATLASTOPDIFF8TEVTPT.dat";
  f1.open(datafile.str().c_str(), ios::in);
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
  }
  string line;
  for (int i = 0; i < fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
            
      int idum=0;
      lstream >> idum; 
      if(idum!=(i+1)){
	cout<<"Error in ATLASTOPDIFF.cc"<<endl;
	cout<<"idum = "<<idum<<endl;
	exit(-10);
      }

      double pt_min=0.0;
      lstream >> pt_min;
      double pt_max=0.0;
      lstream >> pt_max;

      fKin1[i] = 0.5*(pt_min + pt_max); // <pt_top> as characteristic kin variable
      fKin2[i] = Mt;
      fKin3[i] = 8000;

      lstream >> fData[i]; // Central values
      lstream >> fStat[i]; // Total statistical+systematic error (in percent)
      
    }
  
  f1.close();
}

///////////////////////////////////////////////////////////


//-------------------------------------------------------------------------------

// ATLAS top quark differential distributions
// 8 TeV, top quark rapidity distribution
// Using here normalized ditributions, absolute distributions also available
void  ATLASTOPDIFF8TEVTRAPFilter::ReadData()
{
  // Opening files
  fstream f1;

  //
  // The raw data for (1/sigma) dsigma/dpt is taken from
  // https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/TOPQ-2015-06/tabaux_027.pdf
  // Format
  // i  yt_min  yt_max   (1/sigma)dsigma/dyt   Total_Unc(%)
  //
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
	   << fSetName << "/ATLASTOPDIFF8TEVTRAP.dat";
  f1.open(datafile.str().c_str(), ios::in);
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
  }
  cout<<datafile.str()<<endl;
  string line;
  for (int i = 0; i < fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      cout<<line<<endl;
      
      int idum=0;
      lstream >> idum;
      if(idum!=(i+1)){
	cout<<"Error in ATLASTOPDIFF.cc TRAP"<<endl;
	cout<<"idum = "<<idum<<endl;
	exit(-10);
      }

      double yt_min=0.0;
      lstream >> yt_min;
      double yt_max=0.0;
      lstream >> yt_max;

      fKin1[i] = 0.5*(yt_min + yt_max); // <y_top> as characteristic kin variable
      fKin2[i] = Mt;
      fKin3[i] = 8000;

      lstream >> fData[i]; // Central values
      lstream >> fStat[i]; // Total statistical+systematic error (in percent)
      
    }
  
  f1.close();
}

///////////////////////////////////////////////////////////

//-------------------------------------------------------------------------------

// ATLAS top quark differential distributions
// 8 TeV, top quark pair pt distribution
// Normalized to the total cross-section - absolute distributions also available
void  ATLASTOPDIFF8TEVTTPTFilter::ReadData()
{
  // Opening files
  fstream f1;

  //
  // The raw data for (1/sigma) dsigma/dp_tt is taken from
  // https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/TOPQ-2015-06/tabaux_031.pdf
  // Format
  // i  p_tt_min  p_tt_max   (1/sigma)dsigma/dp_tt   Total_Unc(%)
  //
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
	   << fSetName << "/ATLASTOPDIFF8TEVTTPT.dat";
  f1.open(datafile.str().c_str(), ios::in);
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
  }
  string line;
  for (int i = 0; i < fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
            
      int idum=0;
      lstream >> idum; 
      if(idum!=(i+1)){
	cout<<"Error in ATLASTOPDIFF.cc TTPT"<<endl;
	cout<<"idum = "<<idum<<endl;
	exit(-10);
      }

      double ptt_min=0.0;
      lstream >> ptt_min;
      double ptt_max=0.0;
      lstream >> ptt_max;

      fKin1[i] = 0.5*(ptt_min + ptt_max); // <ptt_top> as characteristic kin variable
      fKin2[i] = Mt;
      fKin3[i] = 8000;

      lstream >> fData[i]; // Central values
      lstream >> fStat[i]; // Total statistical+systematic error (in percent)
      
    }
  
  f1.close();
}

///////////////////////////////////////////////////////////

//-------------------------------------------------------------------------------

// ATLAS top quark differential distributions
// 8 TeV, top quark pair  rapidity distribution
// Using normalized cross-sections, absolute distributions also available
void  ATLASTOPDIFF8TEVTTRAPFilter::ReadData()
{
  // Opening files
  fstream f1;

  //
  // The raw data for (1/sigma) dsigma/dpt is taken from
  // https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/TOPQ-2015-06/tabaux_033.pdf
  // Format
  // i  y_tt_min  y_tt_max   (1/sigma)dsigma/dy_tt   Total_Unc(%)
  //
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
	   << fSetName << "/ATLASTOPDIFF8TEVTTRAP.dat";
  f1.open(datafile.str().c_str(), ios::in);
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
  }
  string line;
  for (int i = 0; i < fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
            
      int idum=0;
      lstream >> idum; 
      if(idum!=(i+1)){
	cout<<"Error in ATLASTOPDIFF.cc TTRAP"<<endl;
	cout<<"idum = "<<idum<<endl;
	exit(-10);
      }

      double y_tt_min=0.0;
      lstream >> y_tt_min;
      double y_tt_max=0.0;
      lstream >> y_tt_max;

      fKin1[i] = 0.5*(y_tt_min + y_tt_max); // <y_tt> as characteristic kin variable
      fKin2[i] = Mt;
      fKin3[i] = 8000;

      lstream >> fData[i]; // Central values
      lstream >> fStat[i]; // Total statistical+systematic error (in percent)
      
    }
  
  f1.close();
}

///////////////////////////////////////////////////////////

//-------------------------------------------------------------------------------

// ATLAS top quark differential distributions
// 8 TeV, top quark pair invariant mass distribution
// Normalized to the total cross-section - absolute distributions also available
void  ATLASTOPDIFF8TEVTTMFilter::ReadData()
{
  // Opening files
  fstream f1;

  //
  // The raw data for (1/sigma) dsigma/dmtt is taken from
  // https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/TOPQ-2015-06/tabaux_029.pdf
  // Format
  // i   mtt_min  mtt_max   (1/sigma)dsigma/dmtt   Total_Unc(%)
  //
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
	   << fSetName << "/ATLASTOPDIFF8TEVTTM.dat";
  f1.open(datafile.str().c_str(), ios::in);
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
  }
  string line;
  cout<<fNData<<endl;
  for (int i = 0; i < fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      
      int idum=0;
      lstream >> idum; 
      if(idum!=(i+1)){
	cout<<"Error in ATLASTOPDIFF.cc TTM"<<endl;
	cout<<"idum = "<<idum<<endl;
	exit(-10);
      }

      double mtt_min=0.0;
      lstream >> mtt_min;
      double mtt_max=0.0;
      lstream >> mtt_max;

      fKin1[i] = 0.5*(mtt_min + mtt_max); // <mtt_top> as characteristic kin variable
      fKin2[i] = Mt;
      fKin3[i] = 8000;

      lstream >> fData[i]; // Central values
      lstream >> fStat[i]; // Total statistical+systematic error (in percent)
      
    }
    
  f1.close();
}

///////////////////////////////////////////////////////////

