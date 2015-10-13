/**
 *     EXPeriment      = SLAC E49, E61, E87, E89
 *     REACtion        = e proton --> e X
 *     Author          = Whitlow
 *     REFerence       = PL B282(92)475 and SLAC-357 (1990) (Ph.D)
 *     Additional info : Re-analysed data - BCDMS binning
 *     The proton structure functions from a re-analysis of the
 *     SLAC experiments E49, E61, E87, and E89 in electron proton
 *     deep inelastic scattering binned in the same x ranges as
 *     the BCDMS data.
 *
 *     mean          mean            F2              errors
 *     x           Q**2                       stat.      sys.
 *
 *     EXPeriment      = SLAC E49, E61, E87, E89
 *     REACtion        = e deuterium --> e X
 *     Plab            = 4.50-24.50 GeV
 *     Author          = Whitlow
 *     REFerence       = SLAC-357 (1990) (Ph.D)
 *     Additional info : Re-analysed data - BCDMS binning
 *     The nucleon structure functions from a re-analysis of the
 *     SLAC experiments E49, E61, E87, E89, E137 and E140 in electron
 *     deuterium deep inelastic scattering binned in the same
 *     x ranges as the BCDMS data.
 *
 *     x           Q**2             F2             errors
 *     (GeV**2)                     stat.      sys.
 *
 *     The systematic error is taken fully correlated for all the data points p and d
 *
 *     The absolute normalization is taken to be 2.1% for p and 1.7% for d
 *     The relative normalization of 1.1% is taken into account
 *
 *     NB: There are not info about y and since the data points that were taken at
 *     diffferent beam energies are merged together in the same sample
 *     there are no way to evaluate it. Thus it's set arbitrary equal to 0
 *     simply to fill somehow the entry in datawarehouse.res
 */

#include "SLAC.h"

void SLACPFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/slac_p.data";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  double relnor;
  
  relnor = 1.1*0.5;    //relative normalisation of 1.1% between targets
  
  string line;
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    lstream >> fKin1[i];  //x
    lstream >> fKin2[i];  //q2
    
    lstream >> fData[i];  //obs 
    
    //  SLAC gives errors in absolute value
    //  and we assume the sys errors to be uncorrelated in order
    //  to overestimate them and give less weight
    //  to these data points in evaluating chi2

    lstream >> fStat[i];
    
    lstream >> fSys[i][0].add;
    fSys[i][0].mult = fSys[i][0].add*100/fData[i];
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    
    fSys[i][1].mult = 2.1;  //absnorm
    fSys[i][1].add = fSys[i][0].mult*fData[i]*1e-2;
    fSys[i][1].type = MULT;
    fSys[i][1].name = "SLACNORM";
    
    fSys[i][2].mult = relnor;     //relnorm
    fSys[i][2].add = fSys[i][1].mult*fData[i]*1e-2;
    fSys[i][2].type = MULT;
    fSys[i][2].name = "SLACRELNORM";
  }
  
  f1.close();
}


/**
 *     EXPeriment      = SLAC E49, E61, E87, E89
 *     REACtion        = e proton --> e X
 *     Author          = Whitlow
 *     REFerence       = PL B282(92)475 and SLAC-357 (1990) (Ph.D)
 *     Additional info : Re-analysed data - BCDMS binning
 *     The proton structure functions from a re-analysis of the
 *     SLAC experiments E49, E61, E87, and E89 in electron proton
 *     deep inelastic scattering binned in the same x ranges as
 *     the BCDMS data.
 *
 *     mean          mean            F2              errors
 *     x           Q**2                       stat.      sys.
 *
 *     EXPeriment      = SLAC E49, E61, E87, E89
 *     REACtion        = e deuterium --> e X
 *     Plab            = 4.50-24.50 GeV
 *     Author          = Whitlow
 *     REFerence       = SLAC-357 (1990) (Ph.D)
 *     Additional info : Re-analysed data - BCDMS binning
 *     The nucleon structure functions from a re-analysis of the
 *     SLAC experiments E49, E61, E87, E89, E137 and E140 in electron
 *     deuterium deep inelastic scattering binned in the same
 *     x ranges as the BCDMS data.
 *
 *     x           Q**2             F2             errors
 *     (GeV**2)                     stat.      sys.
 *
 *     The systematic error is taken fully correlated for all the data points p and d
 *
 *     The absolute normalization is taken to be 2.1% for p and 1.7% for d
 *     The relative normalization of 1.1% is taken into account
 *
 *     NB: There are not info about y and since the data points that were taken at
 *     diffferent beam energies are merged together in the same sample
 *     there are no way to evaluate it. Thus it's set arbitrary equal to 0
 *     simply to fill somehow the entry in datawarehouse.res
 */
void SLACDFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/slac_d.data";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  double relnor;
  
  relnor = 1.1*0.5;    //relative normalisation of 1.1% between targets
  
  string line;
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    lstream >> fKin1[i];  //x
    lstream >> fKin2[i];  //q2
    
    lstream >> fData[i];  //obs 
    
    //  SLAC gives errors in absolute value
    //  and we assume the sys errors to be uncorrelated in order
    //  to overestimate them and give less weight
    //  to these data points in evaluating chi2
   
    lstream >> fStat[i];
    
    lstream >> fSys[i][0].add;
    fSys[i][0].mult = fSys[i][0].add*100/fData[i];
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    
    fSys[i][1].mult = 1.7;  //absnorm
    fSys[i][1].add = fSys[i][0].mult*fData[i]*1e-2;
    fSys[i][1].type = MULT;
    fSys[i][1].name = "SLACNORM";
    
    fSys[i][2].mult = -relnor;     //relnorm
    fSys[i][2].add = fSys[i][1].mult*fData[i]*1e-2;
    fSys[i][2].type = MULT;
    fSys[i][2].name = "SLACRELNORM";
  }
  
  f1.close();
}

