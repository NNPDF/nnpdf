/*
  NOTE ADDED by ERN January 2021

  The implementation of the CDF Z rapidity distributions has been updated to 
  reflect v4 in the arXiv. The difference w.r.t. the previous implementation
  consists in th fact that the last two bins have been collapsed together;
  central values and ucnertainties are also slightly different.
  PLEASE note that the CDFZRAP data set is DEPRECATED in favour of the
  CDFRAP_NEW data set
*/


/**
 * ABULENCIA 07 PR D75,092006
 * Measured inclusive jet differential cross section as a function of PT
 * with the jet resolution parameter D = 0.7.
 * PBAR P --> JET X SQRT(S) IN GEV 1960
 *
 * Run II inclusive jet cross sections from the CDF collaboration
 * at the Tevatron with the inclusive kt algorithm
 *
 * All jet-related data can be found at the Durham data base
 *
 * http://durpdg.dur.ac.uk/cgi-  bin/hepdata/reacsearch/TESTREAC/fsp+jet%25+or+fsp+2jet%25+or+fsp+3jet%25+or+fsp+4jet%25/NORMAL&skip=0
 *
 * This particular data set used measurements published in
 * PR D75,092006 (2007) 	 Preprinted as HEP-EX/0701051
 *
 * The data has been obtained from HEPDATA as well:
 * http://durpdg.dur.ac.uk/cgi-bin/hepdata/tabkum/TABLE/11780/999/1/1
 *
 * Format
 *
 *     NBINRAP
 *     do IBINRAP=1,NBINRAP
 *
 *       ABS(YRAP(P=3))_low , ABS(YRAP(P=3))_high
 *         - x = PT(P=3) IN GEV
 *         - y = D2(SIG)/DPT/DYRAP IN NB/GEV
 *       NDATBIN(IBINRAP)
 *       do IDAT = 1,NDATBIN(IBINRAP)
 *          x1   x2   y  dstat   +dsys  -dsys)
 *       enddo
 *     enddo
 *
 * Note that for this data set, the uncertainties are
 * given in absolute value
 * they need to be converted to percentage values
 *
 * The details of the experimental correlation matrix
 * are given in  PR D75,092006 (2007), in particular in
 * Appendix A and in Section X
 *
 * The various sources of systematinc uncertainties are ...
 *
 * 1.- Absolute jet energy scale in the calorimeter
 * (dominant systematic uncertainty)
 *
 * 2.- Uncertainty in definition of exclusive dijet sample
 *
 * 3.- Jet momentum resolution
 *
 * 4.- Modellization of the parton cascades for the
 *     unfolding procedure
 *
 * For the first fits, add systematic and statistical
 * errors in quadrature
 *
 * Contact Person: Mario Martinez-Perez from IFAE
 *
 * The systematic uncertainty from the total integrated
 * luminosity is 5.8%
 *
 * The 24 independent sources of systematic uncertaities are
 * 4 from energy scale (table VIII )
 * 3 beta_data/beta_mc * 5 rapidity ranges = 15
 * 1 resolution
 * 1 unfolding
 * 1 pt-spectra
 * 1 delta_mi
 * 1 Chad in pQCD
 *
 * Plus the luminosity normalization uncertainty which
 * has a different treatment than the rest of the
 * correlated systematic uncertainties
 *
 * Therefore the labelling of systematic uncertainties is
 *
 *     1.- JES-A
 *     2.- JES-B
 *     3.- JES-C
 *     4.- JES-D
 *     5.- BETA_A-rap1
 *     6.- BETA_B-rap1
 *     7.- BETA_C-rap1
 *     8.- BETA_A-rap2
 *     9.- BETA_B-rap2
 *     10.- BETA_C-rap2
 *     11.- BETA_A-rap3
 *     12.- BETA_B-rap3
 *     13.- BETA_C-rap3
 *     14.- BETA_A-rap4
 *     15.- BETA_B-rap4
 *     16.- BETA_C-rap4
 *     17.- BETA_A-rap5
 *     18.- BETA_B-rap5
 *     19.- BETA_C-rap5
 *     20.- RES
 *     21.- UNFOL
 *     22.- PTSPEC
 *     23.- DELTA_MI
 *     24.- CHAD
 *
 * Recall than systematic uncertainties are saved in percent,
 * and that sometimes asymmetric combinations are required
 *
 */

#include "CDF.h"

void CDFR2KTFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3, f4;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/CDF-RunII-kt.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream sys("");
  sys << dataPath() << "rawdata/"
  << fSetName << "/CDF-RunII-kt-TOT.sys";
  f2.open(sys.str().c_str(), ios::in);

  if (f2.fail()) {
    cerr << "Error opening data file " << sys.str() << endl;
    exit(-1);
  }

  stringstream jesfile("");
  jesfile << dataPath() << "rawdata/"
  << fSetName << "/CDF-RunII-kt-JES.sys";
  f3.open(jesfile.str().c_str(), ios::in);

  if (f3.fail()) {
    cerr << "Error opening data file " << jesfile.str() << endl;
    exit(-1);
  }

  stringstream chadfile("");
  chadfile << dataPath() << "rawdata/"
  << fSetName << "/CDF-RunII-kt-CHAD.sys";
  f4.open(chadfile.str().c_str(), ios::in);

  if (f4.fail()) {
    cerr << "Error opening data file " << chadfile.str() << endl;
    exit(-1);
  }

  // Starting filter
  string line;
  double etamin,eta,ptmin,ptmax;
  int nrapbin, ndatbin[5];
  double s = 1960;

  int index = 0;
  getline(f1,line);
  istringstream lstream(line);
  lstream >> nrapbin;
  for (int i = 0; i < nrapbin; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    lstream >> etamin >> eta;
    eta = (eta+etamin)*0.5;

    getline(f1,line);
    istringstream lstream2(line);
    lstream2 >> ndatbin[i];
    for (int j = 0; j < ndatbin[i]; j++)
    {
      getline(f1,line);
      istringstream lstream(line);

      fKin1[index] = eta;                   //eta
      lstream >> ptmin >> ptmax;
      fKin2[index] = (ptmax+ptmin)*(ptmax+ptmin)*0.25;    //pt2
      fKin3[index] = s;                     //sqrt(s)

      lstream >> fData[index];             //obs
      lstream >> fStat[index];             //stat

      index++;
    }
  }

   if (index != fNData)
  {
    cerr << "Mismatch in the number of data points in CDFR2KT" << endl;
    exit(-1);
  }

  // Pre-zero all rapidity-bin-specific systematic errors
  // The following section only fills the systemaics for the relevant
  // rapidity bins, so we need this to initialise all the values to zero.
  for (int i=0; i<fNData; i++)
    for (int k = 5; k <= 19; k++)
        fSys[i][k].mult = 0;
  /*
   Reading decomposition of the systematic uncertainties associated
   to the jet energy scale
   Note that this table gives the percentage contribution to
   each of the four sources of correlated uncertainty which contribute
   to the systematic uncertainty associated to Jet Energy Scale (JES)
 */
  double jes[17][4];
  double sum;
  for (int j = 0; j < 17; j++)
  {
    getline(f3,line);
    istringstream lstream(line);
    lstream >> ptmin >> ptmax;

    sum = 0;
    for (int k = 0; k < 4; k++)
    {
      lstream >> jes[j][k];
      jes[j][k]*=1e-2;
      sum+= jes[j][k]*jes[j][k];
    }

    if (fabs(sum-1e0) >= 5e-3)
    {
      cerr << "Problem in JES systematic for CDF R2 kt!" << endl;
      cerr << "SUM = " << sum << endl;
      exit(-1);
    }
  }

  // Get systematic uncertainties
  /* The format of the systematic uncertainties files is as follows
   See Table II-III of CDF paper

   PTMIN PTMAX JES+ JES- beta_a+ beta_a- beta_b+ beta_b-
   beta_c+ beta_c- Res+ Res- Unfold Ptspect dmi+ dmi-

   Note that systematic uncertainties are given in percentage
   */
  double dtmp,stmp,up,down;
  double shift[fNData];

  index = 0;
  getline(f2,line);
  for (int i = 0; i < nrapbin; i++)
  {
    getline(f2,line);
    getline(f2,line);
    for (int j = 0; j < ndatbin[i]; j++)
    {

      getline(f2,line);
      istringstream lstream(line);
      lstream >> ptmin >> ptmax;

      shift[index] = 0;

      lstream >> up >> down;         // Get symmetrized systematic errors for JES
      symmetriseErrors(up, down, &stmp, &dtmp);
      shift[index]+= dtmp;
      for (int k = 0; k < 4; k++)
        fSys[index][k+1].mult = stmp*jes[j][k];

      for (int k = 0; k < 3; k++)
      {
        lstream >> up >> down;       // Get symmetrized systematic error for beta ratios
        symmetriseErrors(up, down, &stmp, &dtmp);
        fSys[index][5+3*i+k].mult = stmp;
        shift[index]+= dtmp;
      }

      lstream >> up >> down;         // Get systematic error for resolution
      symmetriseErrors(up, down, &stmp, &dtmp);
      fSys[index][20].mult = stmp;
      shift[index]+= dtmp;

      lstream >> fSys[index][21].mult;  // Get systematic error for unfolding
      lstream >> fSys[index][22].mult;  // Get systematic error for ptspectra

      lstream >> up >> down;         // Get systematic error for delta-MI
      symmetriseErrors(up, down, &stmp, &dtmp);
      fSys[index][23].mult = stmp;
      shift[index]+= dtmp;

      index++;
    }
  }

  /*
   Get corrections due to the hadronization effects
   and their associated systematic uncertainty
   The systematic uncertainty associated to hadronization effects
   is only correlated within all the pt bins within each
   of the various rapidity ranges
   */
  double chad, syschad;

  index = 0;
  getline(f4,line);
  for (int i = 0; i < nrapbin; i++)
  {
    getline(f4,line);
    getline(f4,line);
    for (int j = 0; j < ndatbin[i]; j++)
    {
      getline(f4,line);
      istringstream lstream(line);
      lstream >> ptmin >> ptmax;

      lstream >> chad >> syschad;

      // divide data to np corrections
      fData[index] /= chad;
      fStat[index] /= chad;

      fSys[index][24].mult = syschad/chad*100;

      index++;
    }
  }

  for (int i = 0; i < fNData; i++)
  {
    // Normalization uncertainty from integrated luminosity
    fSys[i][0].mult = 5.8;
    fSys[i][0].add = fSys[i][0].mult*fData[i]*1e-2;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CDFLUMI";

    for (int l = 1; l < fNSys; l++)
    {
      fSys[i][l].add = fSys[i][l].mult*fData[i]*1e-2;
      fSys[i][l].type = MULT;
      fSys[i][l].name = "CORR";
    }

    // Shift due to the presence of asymmetric systematic errors
    fData[i] *= (1.0 + shift[i]*0.01);

  }

  f1.close();
  f2.close();
  f3.close();
  f4.close();

}



/**
 *
 *     cdf-Zrapdist.f
 *
 *     Z rapidity distribution measurements at the Tevatron
 *     by the CDF collaboration
 *
 *     Actually the measured cross section is for the reaction
 *     dsigma/dy in ppbar -> Z0/gamma* -> e+e-
 *     in the mass range 66 < M_{Z0/gamma*} < 116
 *
 *     The data has been extracted from
 *
 *     "dsigma/dy distribution of Drell-Yan dielectron pairs"
 *     by J. Han, A. Bodek, W. Sakumoto, Y. Chung (Rochester University)
 *     CDF internal note, updated May 1 2008
 *
 *     There is an update from
 *
 *     dsigma/dy Distribution of Drell-Yan Dielectron Pairs with 2.1 fb -1
 *     Authors: Jiyeon Han, Arie Bodek, Willis Sakumoto, Y.S.Chung
 *     Last Updated: May, 14th, 2009
 *
 *     but no official publication yet
 *
 *     The full correlation matrix is available, with a total of 11
 *     independent correlated systematic uncertainties.
 *
 *     1.- Luminosity
 *     2.- Background on Z (CC)
 *     3.- Background on Z (CP)
 *     4.- Background on Z (PP)
 *     5.- Central electron ID efficiency
 *     6.- Plug electron ID efficiency
 *     7.- Central material effect
 *     8.- Plug material effect
 *     9.- Z vertex finding efficiency
 *     10.- Electron tracking efficiency
 *     11.- No tracks events fraction
 *
 *     The systematic correlated uncertainty is essentially
 *     dominated by the luminosity uncertainty, but in any
 *     case the experimental correlation matrix will be
 *     constructed in the proper way
 *
 *     Run II data, so sqrt(s)=1960 GeV
 *
 *     Data in file:
 *     data/tev-Zrapdist/CDF-Zrapdist.data (from the 2008 update)
 *
 *     with the format
 *
 *     yZ  dsigma/dyZ/0.1 [pb]   stat  sys
 *
 *     Note that differently from D0, these data set does not
 *     divide for the total integrated cross section for Z production!
 *
 *     where the  errors are given in abolute value
 *     so that they must be converted to percentage
 *
 *     For this process the distribution is symmetric in rapidity, so
 *     only positive rapidities are considered
 *
 *     An additional 6% luminosity uncertainty has to be considered,
 *     as is the default by the tevatron
 *
 *
 *     Published data from D0 can be found in the following reference:
 *
 *     Published data can be found in
 *     Subjects: 	High Energy Physics - Experiment (hep-ex)
 *     Journal reference: 	Phys.Rev.D76:012003,2007
 *     DOI: 	10.1103/PhysRevD.76.012003
 *     Report number: 	Fermilab-Pub-07-040-E
 *     Cite as: 	arXiv:hep-ex/0702025v1
 *
 *     which includes the full experimental systematic covariance
 *     matrix, and which should probably be added as a
 *     new experiment in the analysis
 *
 */

/* Note added by ERN on Janaury 2021
  
   The data set has been updated to reflect the results published in
      Phys. Lett. B692 232 (2010) 
   The last two rapidity bins have been collapsed together; central values and
   uncertainties have been updated. The full breakdown of systematic 
   uncertainties has been produced with the error_propagator_g++.C provided
   from the following link 
     https://www-cdf.fnal.gov/physics/ewk/2009/dszdy/dszdy_sys.htm
   From NNPDF4.0 onwards, the use of the CDFZRAP data set is DEPRECATED in
   favour of the CDFZRAP_NEW data set.
*/

void CDFZRAPFilter::ReadData()
{
  // Opening files
  fstream f1;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/CDF-Zrapdist-2009.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Starting filter
  string line;
  const double MZ2 = pow(MZ, 2.0);
  const double s = 1960;
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);

    lstream >> fKin1[i];   //y
    fKin2[i] = MZ2;   //Mass Z squared
    fKin3[i] = s;     //sqrt(s)

    lstream >> fData[i];
    lstream >> fStat[i];

    fSys[i][0].mult = 6.0;  //luminosity
    lstream >> fSys[i][0].add;         //absolute value of 6% luminosity
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CDFLUMI";

    for (int isys = 1; isys < fNSys; isys++)
    {
      lstream >> fSys[i][isys].add;      //systematics
      fSys[i][isys].mult = fSys[i][isys].add*100/fData[i];
      fSys[i][isys].type = MULT;
      fSys[i][isys].name = "CORR";
    }
  }

  f1.close();
}

void CDFZRAP_NEWFilter::ReadData()
{
  // Opening files
  fstream f1;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/CDFZRAP/CDF-Zrapdist-2009_new.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Starting filter
  string line;
  const double MZ2 = pow(MZ, 2.0);
  const double s = 1960;
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);

    lstream >> fKin1[i];   //y
    fKin2[i] = MZ2;   //Mass Z squared
    fKin3[i] = s;     //sqrt(s)

    lstream >> fData[i];
    lstream >> fStat[i];

    //Luminosity uncertainty (6%)
    lstream >> fSys[i][0].add;
    fSys[i][0].mult = fSys[i][0].add/fData[i]*100.;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CDFLUMI";

    for (int isys = 1; isys < fNSys; isys++)
    {
      lstream >> fSys[i][isys].add;      //systematics
      fSys[i][isys].mult = fSys[i][isys].add*100/fData[i];
      fSys[i][isys].type = MULT;
      fSys[i][isys].name = "CORR";
    }
  }

  f1.close();
}

/**
 *
 *     Direct W rapidity asymmetry measurements at the Tevatron
 *
 *     Includes the same information as previous measurements of
 *     the electron rapidity asymmetry, but it is easier to
 *     compute theoretically
 *
 *     The asymmetry is defined as
 *
 *     A(yW) = ( dsigma(W+)/dyW - dsigma(W-)/dyW ) /
 *             ( dsigma(W+)/dyW + dsigma(W-)/dyW )
 *
 *     which in the leading order parton model can be written as
 *
 *     A(yW) = ( u(x1)d(x2) - d(x1)u(x2) ) /
 *             ( u(x1)d(x2) + d(x1)u(x2) )
 *
 *     where x_{1,2} = ( MW/sqrt(s) ) * exp( +- yW ) again
 *     for leading order kinematics
 *
 *     The data has been extracted from
 *
 *     "Direct Measurement of W Boson Charge Asymmetry
 *     with 1 1/fb of Run II Data"
 *     CDF note 8942
 *
 *     Data in file:
 *     data/tev-Wasymm/tev-Wasymm.data
 *
 *     with the format
 *
 *     yWmin yWmax <yW>   A(yW)  sys  sys+stat
 *
 *     where the  errors are given in abolute value
 *     so that they must be converted to percentage
 *
 *
 *     An update of this data set has been presented in
 *     arXiv:0901.2169 [hep-ex]
 *     "Direct Measurement of the W Production Charge Asymmetry in p anti-p Collisions at s**(1/2) = 1.96-TeV"
 *     FERMILAB-PUB-09-017-E
 *
 *     where the latest CDF results for the W production asymmetry
 *     are presented, including the full covariance matrix
 *     of the systematic uncertainties
 *
 *     The covariance matrix is given in the file
 *
 *     data/tev-Wasymm/tev-Wasymm-sys.data
 *
 *     with the format
 *
 *     yWmin yWmax (syst(ii),ii=1,7) stat
 *
 *     for all seven systematic uncertainties in the measurement
 *     which need to be multiplied by 0.01 to get their
 *     absolute value
 *
 *     The 7 systematic correlated uncertainties are
 *     1.- Charge MisID
 *     2.- Backgrounds
 *     3.- Energy scales and resolution
 *     4.- Recoil model
 *     5.- Electron trigger
 *     6.- Electron ID
 *     7.- PDFs
 *
 *     ------------------
 *
 *     We need to check if there are analogous results from DO
 *     and if in global PDF fits this is the only relevant data set
 *
 *     Electroweak convenor at D0: Heidi Schelmann, schellman@fnal.gov
 *
 *     For the W asymmetry, D0 has only measurements of the lepton
 *     decay distribution, and not of the W asymmetry itself like CDF,
 *     and are described in two publications:
 *     arXiv:0807.3367v1 [hep-ex]
 *     arXiv:0709.4254.
 *     One is muon and the other is electron
 *     - they are essentially completely independent.
 *
 */
void CDFWASYMFilter::ReadData()
{
  // Opening files
  fstream f1, f2;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/tev-Wasymm.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream datafile2("");
  datafile2 << dataPath() << "rawdata/"
  << fSetName << "/tev-Wasymm-sys.data";
  f2.open(datafile2.str().c_str(), ios::in);

  if (f2.fail()) {
    cerr << "Error opening data file " << datafile2.str() << endl;
    exit(-1);
  }

  // Starting filter
  string line;
  double tmp;
  const double MW2 = pow(MW, 2.0);
  const double s = 1960;
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);       //tev-Wasymm.data
    istringstream lstream(line);

    lstream >> tmp >> tmp;  //bottom and top of y bins

    lstream >> fKin1[i];    //<y>
    fKin2[i] = MW2;         //Mass W squared
    fKin3[i] = s;           //sqrt(s)

    lstream >> fData[i];

    getline(f2,line);       //tev-Wasymm-sys.data
    istringstream lstream2(line);

    lstream2 >> tmp >> tmp;  //bottom and top of y bins

    for (int isys = 0; isys < fNSys; isys++)
    {
      lstream2 >> fSys[i][isys].add;      //systematics
      fSys[i][isys].add*= 0.01;          //values in file are (x10^-2)
      fSys[i][isys].mult = fSys[i][isys].add*100/fData[i];  //convert to percentage
      fSys[i][isys].type = ADD;
      fSys[i][isys].name = "CORR";
    }

    lstream2 >> fStat[i];   //statistical uncertainty (abs)
    fStat[i]*= 0.01;
  }

  f1.close();
}

/**
 *
 *     Direct W rapidity asymmetry measurements at the Tevatron
 *
 *     Includes the same information as previous measurements of
 *     the electron rapidity asymmetry, but it is easier to
 *     compute theoretically
 *
 *     The asymmetry is defined as
 *
 *     A(yW) = ( dsigma(W+)/dyW - dsigma(W-)/dyW ) /
 *             ( dsigma(W+)/dyW + dsigma(W-)/dyW )
 *
 *     which in the leading order parton model can be written as
 *
 *     A(yW) = ( u(x1)d(x2) - d(x1)u(x2) ) /
 *             ( u(x1)d(x2) + d(x1)u(x2) )
 *
 *     where x_{1,2} = ( MW/sqrt(s) ) * exp( +- yW ) again
 *     for leading order kinematics
 *
 *     The data has been extracted from
 *
 *     "Direct Measurement of W Boson Charge Asymmetry
 *     with 1 1/fb of Run II Data"
 *     CDF note 8942
 *
 *     Data in file:
 *     data/tev-Wasymm/tev-Wasymm.data
 *
 *     with the format
 *
 *     yWmin yWmax <yW>   A(yW)  sys  sys+stat
 *
 *     where the  errors are given in abolute value
 *     so that they must be converted to percentage
 *
 *
 *     An update of this data set has been presented in
 *     arXiv:0901.2169 [hep-ex]
 *     "Direct Measurement of the W Production Charge Asymmetry in p anti-p Collisions at s**(1/2) = 1.96-TeV"
 *     FERMILAB-PUB-09-017-E
 *
 *     where the latest CDF results for the W production asymmetry
 *     are presented, including the full covariance matrix
 *     of the systematic uncertainties
 *
 *     The covariance matrix is given in the file
 *
 *     data/tev-Wasymm/tev-Wasymm-sys.data
 *
 *     with the format
 *
 *     yWmin yWmax (syst(ii),ii=1,7) stat
 *
 *     for all seven systematic uncertainties in the measurement
 *     which need to be multiplied by 0.01 to get their
 *     absolute value
 *
 *     The 7 systematic correlated uncertainties are
 *     1.- Charge MisID
 *     2.- Backgrounds
 *     3.- Energy scales and resolution
 *     4.- Recoil model
 *     5.- Electron trigger
 *     6.- Electron ID
 *     7.- PDFs
 *
 *     ------------------
 *
 *     We need to check if there are analogous results from DO
 *     and if in global PDF fits this is the only relevant data set
 *
 *     Electroweak convenor at D0: Heidi Schelmann, schellman@fnal.gov
 *
 *     For the W asymmetry, D0 has only measurements of the lepton
 *     decay distribution, and not of the W asymmetry itself like CDF,
 *     and are described in two publications:
 *     arXiv:0807.3367v1 [hep-ex]
 *     arXiv:0709.4254.
 *     One is muon and the other is electron
 *     - they are essentially completely independent.
 *
 */
