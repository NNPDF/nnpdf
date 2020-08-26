/*
WARNING:
THE USE OF THIS FILTER IS DEPRECATED IN FAVOUR OF NUTECFe.cc
*/

/**
 *
 *     nutev-dimuon.f
 *
 *     The NuteV dimuon cross sectio measurements
 *
 *     The NuTev dimuon table is the cross section
 *     of neutrino and antineutrino produced charged current dimuons
 *     from charm production such that the muon from the decay of the
 *     charmed meson has an energy greater than 5 GeV (hence "forward").
 *     The data has been corrected for detector smearing and acceptance
 *     effects, as well as for the background due to semileptonic pion
 *     and kaon decays.  Rather than a model dependent charm production
 *     cross section, we are providing a model independent neutrino DIMUON
 *     production cross section.  This has the advantage that users
 *     of the data are not subjected to our choice of fragmentation and
 *     decay models.  They are however subject to the freedom and
 *     responsibility of providing their own.  Further details regarding the
 *     table and its use are outlined in the references below.  They are
 *     included with this package as gentle encouragement that they be
 *     read.  It should be noted that the results included with this package
 *     supercede the NuTeV cross section values included in reference 1.)
 *
 *     References:
 *
 *     1.)         M. Goncharov et al, PRECISE MEASUREMENT OF DIMUON
 *                  PRODUCTION CROSS-SECTIONS IN MUON NEUTRINO
 *                  FE AND MUON ANTI-NEUTRINO FE DEEP INELASTIC
 *                  SCATTERING AT THE TEVATRON
 *                  Phys.Rev.D64:112006,2001, hep-ex/0102049
 *
 *     2.)         D. Mason, NUTEV STRANGE/ANTISTRANGE SEA MEASUREMENTS
 *                  FROM NEUTRINO CHARM PRODUCTION
 *                  DIS 2005 Proceedings pp 851-854
 *
 *
 *
 *     NuTeVtable.dat        the data
 *     NuTeVdimuprd.ps       reference for the data
 *     IS05proceedings.ps   reference for this update
 *
 *     Data file:
 *
 * CROSS SECTION VALUES REPORTED HERE ARE:
 *
 *                PI          dsigma (E,x,y)
 *     100 * ------------ *  ---------------
 *           G^2 M_p E_nu         dx dy
 *
 *     The data file is formatted:  'f3.0,f7.2,2f6.3,13f8.4'
 *     There are 90 (2 polarities * 3 nu energies * 3 y * 5 x) data points.
 *     The data columns are as follows:
 *
 *     1.    Polarity, 1=neutrino, 2=antineutrino
 *     2.    Neutrino energy of bin
 *     3.    y, inelasticity of bin
 *     4.    Bjorken x of bin
 *     5.    Cross section value
 *     6.    Cross section statistical error
 *     7.    Cross section total systematic error
 *     8-15. Individual systematics to make up value in 7. as follows:
 *         8  pi-K model numode
 *         9  pi-K model nubarmode
 *        10  toroid muon energy scale 1%
 *        11  hadronic energy scale 0.5%
 *        12  R_Longitudinal 20%
 *        13  Muon range-out energy 2.5%
 *        14  MC statistics
 *        15  nu/nubar relative normalization 0.67%
 *     16.   Example 5GeV muon energy cut acceptance correction *
 *     17.   DOF for each bin.
 *
 *
 *     Important assumption ->
 *     The included acceptance corrections assume leading order kinematics,
 *     and were calculated within our leading order model (see ref. 1.)
 *
 *     The updated NLO acceptance corrections can be determined
 *     from D. Mason's Ph. D. Thesis. Compare LO and NLO
 *     acceptances for each different value of Q2
 *
 *     How to use the experimental NuTev dimuon data:
 *
 *     Fitting a charm cross section to this table entails modeling
 *     the fragmentation and semi-muonic decay of the charmed particle
 *     to find an acceptance correction for the 5 GeV cut on the decay muon,
 *     as outlined in reference 1.). As an example an acceptance correction
 *     is included in row 16 assuming NuTeV's LO model.  The errors on
 *     each cross section point are constructed so that one may assume the
 *     bins to be independent of each other, and the correlations are taken
 *     into account via the degree of freedom  for each bin.  A good fit
 *     should give you a chi^2 around the sum of the DOF's of each bin, or
 *     around 40.  A chi^2 of 90 is not a good fit.
 *
 *     Questions regarding the use of this table may be addressed to:
 *
 *     David Mason
 *     dmason@jeckyll.uoregon.edu
 *
 *     The number for this observable (Dimuon cross sections from
 *     NuteV) is ITARGET=21 for nuetrinos and ITARGET=22 for
 *     antineutrinos
 *
 *     NuTev target is iron, so it can be considered up to
 *     a good approximation an isoscalar target. Non-isoscalarity
 *     corrections can be easily applied at the level
 *     of PDF decomposition
 *
 *     Note that this experiment must be divided into
 *     two data sets, one for neutrino DIS dimuons and the
 *     other for antineutrino DIS dimuons
 *
 *     The statistical uncertainties in the neutrino set are
 *     much smaller than those in the antineutrino set, since
 *     in the first case the sample is almost five times larger
 *     than in the second
 *
 *     For the time being CHORUS does not provide model-independent
 *     infor on neutrino dimuon cross sections but only
 *     a LO analysis, which cannot be used for us
 *     Ref: Nuclear Physics B 798 2008 1-16
 *
 *     The NuTev dimuon cross section table is indeed shown to be
 *     model independent in D. Mason's Ph.D. thesis, where it
 *     is shown that cross sections extracted from the NLO MC
 *     are essentially the same as those from the LO MC ones,
 *     which are publicly available
 *
 *     The values of Q2 and W2 (required to be known both
 *     to implement acceptance corrections and to check if kinematical
 *     cuts should be imposed also for the NuTeV dimuon cross section
 *     data)
 *     The kinematical relations are the followning
 *     Q2 = s * x * y = 2 * MN * Enu * x * y
 *     W2 = Q2 * ( 1 - x ) / x
 *
 *     To know Q2 is also important for the kinematical
 *     plance coverage plots!
 *
 *     With the standard kinematical cuts, 4 points out of the total
 *     sample of 90 dimuon cross section points will be excluded,
 *     since they correspond to too small values of Q2
 *     Update selkin.f with also Nutev dimuon data included
 *     for the kinematical cuts
 *
 *     From the NuTeV acceptance correction plots we observe that
 *     Q2 is bounded to be in the range
 *     2 GeV2 < Q2 < 250 GeV2
 *     This is certainly wrong, there is a problem with that
 *     plot? Repeat with the acceptance corrections from
 *     the NuTev table!
 *
 *
 *     The kinematic variables are chosen to be (x,Q2,y), in order
 *     to be consistent with the other DIS observables
 *     The conversion between Q2 and E is straightforward, for
 *     fixed target kinematics
 *
 *     Q2 = 2 MN Enu *  ( x * y )
 *     Enu = Q2  / ( 2 * MN * x * y )
 *
 *     Since we are going to use in any case the NLO acceptance
 *     corrections from NuTeV, we correct experimental data
 *     for these at the level of the filtering, to avoind
 *     having to do so at the level of the computation of the
 *     observable
 *
 *     Similar considerations for the branching ratio of averaged
 *     charmed mesons into muons
 *     BRC = 0.099 \pm 0.012, same as used in NuTeV dimuon analysis
 *     from  a reanalysis of 125 charm events measured by FNAL E531
 *
 *     Effects of finite charm mass, which are potentially important,
 *     can be implemented with a slow reescaling as in the
 *     case of the I-ZM-VFN for all other inclusive observables,
 *     just a modification of the x variable
 *
 *     Important to understand if systematic uncertainties are
 *     completely correlated among them, or if we choose them
 *     to be completely uncorrelated
 *     In principle we have the full covariance matrix so we can
 *     choose to keep full dependence with the correlations
 *
 *     Check the effective number of degrees of freedom quoted
 *     by the NuTev people, divided by Ndat it should give
 *     the expected reasonable value of the chi2 when all
 *     uncertainties are added in quadrature.
 *     DOFtot/fNData = 0.4867 for all 90 Nutev dimuon data points
 *         (some of them cut by kinematics)
 *
 *     So finally the Nutev dimuon data takes errors in quadrature
 *
 *     Possibility of adding nuclear corrections
 *
 *     No absolute normalization uncertaintis are mentioned
 *     in the NuTeV dimuon papers, but a 2.1% normalization is
 *     quoted in the NuTeV inclusive structure function papers,
 *     Phys.Rev.D74:012008,2006.
 *     e-Print: hep-ex/0509010
 *
 *     NLO acceptance tables have been provided
 *     by D. Mason both for the NuTeV and the
 *     CCFR data and can be found in the data files
 *
 *     nf20-1.25-0.60.nu.cor
 *     nf20-1.25-0.60.bar.cor
 *     cnf22-1.25-0.60.nu.cor
 *     cnf22-1.25-0.60.bar.cor
 *
 */


#include "NUTEV.h"

void NTVNUDMNFilter::ReadData()
{
  // Opening files
  fstream f1, f2;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/NuTeVtable.dat";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  stringstream datafile2("");
  datafile2 << dataPath() << "rawdata/"
  << fSetName << "/nf20-1.25-0.60.nu.cor";
  f2.open(datafile2.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << datafile2.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  double mn = 0.938;
  double BrC = 0.099;
  
  double acc_cor[fNData];
  string line;

  //Get NLO acceptance tables
  for (int i = 0; i < fNData; i++)
  {
    getline(f2,line);
    istringstream lstream(line);
    lstream >> acc_cor[i];
  } 
  
  double tmp, DOF, enu;
 
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> tmp;  // Neutrino polarity (not used)
    
    lstream >> enu;  //Neutrion energy
    lstream >> fKin3[i];    // Inelasticity y
    lstream >> fKin1[i];   // Bjorken x
    
    
    fKin2[i] = 2.0*mn*enu*fKin1[i]*fKin3[i];    // Q2
    
    // Dimuon neutrino double differencial cross section
    // Corrected for NLO acceptance and branching ratio
    // to yield the NLO charm production cross section
    lstream >> fData[i];
    fData[i] /= acc_cor[i]*BrC;
    
    // Statistical uncertainties  
    lstream >> fStat[i];
    fStat[i] /= acc_cor[i]*BrC;        // Be careful with the acceptance corrections here 
    
    // Systematic uncertainties (treated as uncorrelated)
    lstream >> fSys[i][0].add;
    fSys[i][0].add /= acc_cor[i]*BrC;
    fSys[i][0].mult = fSys[i][0].add*100/fData[i];
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
       
    // Normalization uncertainty
    fSys[i][1].mult = 2.1;
    fSys[i][1].add = fSys[i][1].mult*fData[i]*1e-2;
    fSys[i][1].type = MULT;
    fSys[i][1].name = "NUTEVNORM";
    
    for(int i = 0; i < 8; i++)
      lstream >> tmp;           //Individual systematics (not used)
    
    lstream >> tmp;    //LO acceptance correction (not used)
    lstream >> DOF;    //DOF (not used)
  }
  
  f1.close();
  f2.close();
}

/**
 * See filterNTVNUDMN()
 */
void NTVNBDMNFilter::ReadData()
{
  // Opening files
  fstream f1, f2;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/NuTeVtable.dat";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  stringstream datafile2("");
  datafile2 << dataPath() << "rawdata/"
  << fSetName << "/nf20-1.25-0.60.bar.cor";
  f2.open(datafile2.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << datafile2.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  // Starting filter
  double mn = 0.938;
  double BrC = 0.099;
  
  double acc_cor[fNData];
  string line;

  //Get NLO acceptance tables
  for (int i = 0; i < fNData; i++)
  {
    getline(f2,line);
    istringstream lstream(line);
    lstream >> acc_cor[i];
  } 
  
  double tmp, DOF, enu;
 
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> tmp;  // Neutrino polarity (not used)
    
    lstream >> enu;  //Neutrion energy
    lstream >> fKin3[i];    // Inelasticity y
    lstream >> fKin1[i];   // Bjorken x
    
    
    fKin2[i] = 2.0*mn*enu*fKin1[i]*fKin3[i];    // Q2
    
    // Dimuon neutrino double differencial cross section
    // Corrected for NLO acceptance and branching ratio
    // to yield the NLO charm production cross section
    lstream >> fData[i];
    fData[i] /= acc_cor[i]*BrC;
    
    // Statistical and systematic uncertainties  
    lstream >> fStat[i];
    fStat[i] /= acc_cor[i]*BrC;        // Be careful with the acceptance corrections here 
    
    // Systematic uncertainties (treated as uncorrelated)
    lstream >> fSys[i][0].add;
    fSys[i][0].add /= acc_cor[i]*BrC;
    fSys[i][0].mult = fSys[i][0].add*100/fData[i];
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
       
    // Normalization uncertainty
    fSys[i][1].mult = 2.1;
    fSys[i][1].add = fSys[i][1].mult*fData[i]*1e-2;
    fSys[i][1].type = MULT;
    fSys[i][1].name = "NUTEVNORM";
    
    for(int i = 0; i < 8; i++)
      lstream >> tmp;           //Individual systematics (not used)
    
    lstream >> tmp;    //LO acceptance correction (not used)
    lstream >> DOF;    //DOF (not used)
  }
  
  f1.close();
  f2.close();
}

