/**
 * Run II inclusive jet cross sections from the D0 collaboration
 * at the Tevatron with the cone algorithm (D0 MidPoint)
 *
 * Measured inclusive jet differential cross section as a
 * function of PT for the absolute rapidity
 * range 0.0 to 0.4 with the cone radius R = 0.7.
 * PBAR P --> JET X
 * SQRT(S) IN GEV 1960
 *
 * All jet-related data can be found at the Durham data base
 *
 * http://durpdg.dur.ac.uk/cgi-bin/hepdata/reacsearch/TESTREAC/fsp+jet%25+or+fsp+2jet%25+or+fsp+3jet%25+or+fsp+4jet%25/NORMAL&skip=0
 *
 * This particular data set used measurements published in
 * Preprinted as HEP-EX/0802.2400
 *
 * The use of Run II jet data is better than the Run I one
 * because of the better understanding of the low pt region
 * (cfr. Mario Martinez-Perez)
 *
 * Format
 *
 * NBINRAP
 * do IBINRAP=1,NBINRAP
 *
 * x1, x2:  ABS(YRAP(P=3))_low - ABS(YRAP(P=3))_high
 * - x = PT(P=3) IN GEV
 * - y = D2(SIG)/DPT/DYRAP IN PB/GEV
 * NDATBIN(IBINRAP)
 *   do IDAT = 1,NDATBIN(IBINRAP)
 *     x1   x2      x    y  dstat (%)  +dsys (%) -dsys(&)
 *   enddo
 * enddo
 *
 * On top of the quoted systematic uncertainties, there is a luminosity
 * uncertainty of 6.1% fully correlated in all the bins of y and pt
 *
 * For the first fits add in quadrature systematic and
 * statistical uncertainties
 *
 * The correlation matrix for the systematic uncertainties
 * for this data set is given in the file
 * D0-jets-cone-RunII-sys.data
 *
 * There is a total of 23 systematic uncertainties
 * Description of systematic sources.
 * Source    Description
 * duncorr   Uncorrelated uncertainty
 * dsys001   EM energy scale
 * dsys002   Photon energy scale
 * dsys003   High p_{T} extrapolation
 * dsys004   #eta-intercalibration
 * dsys005   Detector showering
 * dsys006   Luminosity
 * dsys007   #eta-intercalibration
 * dsys008   #eta-intercalibration
 * dsys009   #eta-intercalibration
 * dsys010   JES resolution bias
 * dsys011   Resolution method
 * dsys012   Non-Gaussian tails
 * dsys013   Zero-suppression
 * dsys014   Resolution
 * dsys015   #eta-intercalibration fit
 * dsys016   JES MPF bias
 * dsys017   JES MPF bias
 * dsys018   Rapidity unfolding
 * dsys019   Trigger matching
 * dsys020   Dijet response fit
 * dsys021   Dijet response fit
 * dsys022   Trigger matching
 * dsys023   CC response fit
 *
 * The full correlation matrix of this experiment is easily available
 * Check that total quoted systematic uncertainty corresponds
 * to the sum in quadrature of all indepenedent correlated
 * sources of systematic errors
 *
 * Format of file "D0-jets-cone-RunII-sys.data"
 *
 * NBINRAP
 * do IBINRAP=1,NBINRAP
 *
 * ABS(YRAP(P=3))_low   ABS(YRAP(P=3))_high
 * NDATBIN(IBINRAP)
 *  do IDAT = 1,NDATBIN(IBINRAP)
 *    pt(min)   pt(max)   duncorr(%,+-)   dsys001(%,+-)   dsys002(%,+-)   dsys003(%,+-)   dsys004(%,+-)   dsys005(%,+-)
 *  enddo
 *  do IDAT = 1,NDATBIN(IBINRAP)
 *    pt(min)   pt(max)   dsys006(%,+-)   dsys007(%,+-)   dsys008(%,+-)   dsys009(%,+-)   dsys010(%,+-)   dsys011(%,+-)
 *  enddo
 *  do IDAT = 1,NDATBIN(IBINRAP)
 *    pt(min)   pt(max)   dsys012(%,+-)   dsys013(%,+-)   dsys014(%,+-)   dsys015(%,+-)   dsys016(%,+-)   dsys017(%,+-)
 *  enddo
 *  do IDAT = 1,NDATBIN(IBINRAP)
 *    pt(min)   pt(max)   dsys018(%,+-)   dsys019(%,+-)   dsys020(%,+-)   dsys021(%,+-)   dsys022(%,+-)   dsys023(%,+-)
 *  enddo
 * enddo
 *
 * Therefore we have 23 sources of correlated systematic
 * uncertainties with a single source of uncorrelated
 * systematic uncertainties
 *
 */

#include "D0.h"

void D0R2CONFilter::ReadData()
{
  // Opening files
  fstream f1, f2;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/D0-jets-cone-meas.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream sys("");
  sys << dataPath() << "rawdata/"
  << fSetName << "/D0-jets-cone-RunII-sys.data";
  f2.open(sys.str().c_str(), ios::in);

  if (f2.fail()) {
    cerr << "Error opening data file " << sys.str() << endl;
    exit(-1);
  }

  // Read data file
  string line;
  double etamin,etamax,eta,ptmin,ptmax,ptavg;
  int nrapbin, ndatbin[6];
  double s = 1960;

  int idat = 0;
  getline(f1,line);
  istringstream lstream(line);
  lstream >> nrapbin;
  for (int i = 0; i < nrapbin; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    lstream >> etamin >> etamax;
    eta = (etamax+etamin)*0.5;

    getline(f1,line);
    istringstream lstream2(line);
    lstream2 >> ndatbin[i];
    for (int j = 0; j < ndatbin[i]; j++)
    {
      getline(f1,line);
      istringstream lstream(line);

      fKin1[idat] = eta;                   //eta
      lstream >> ptmin >> ptmax >> ptavg;
      fKin2[idat] = ptavg*ptavg;    //pt2
      fKin3[idat] = s;              //sqrt(s)

      lstream >> fData[idat];             //obs
      lstream >> fStat[idat];             //stat

      idat++;
    }
  }

  // Read experimental systematic correlation matrix
  double sysread[fNData][24];           // array to store systematics temporarily
  double up,down,shift[fNData],dtmp,stmp;

  for (int i = 0; i < fNData; i++)
    shift[i]=0;

  int isys = 0;
  int jdat = 0;
  idat = 0;
  getline(f2,line);
  for (int i = 0; i < nrapbin; i++)
  {
    getline(f2,line);
    getline(f2,line);
    for (int k = 0; k < 4; k++)          //4 sets of 6 systematics
    {
      isys = k*6;
      jdat = idat;
      for (int j = 0; j < ndatbin[i]; j++)
      {
        getline(f2,line);
        istringstream lstream(line);

        lstream >> ptmin >> ptmax;

        for (int l = 0; l < 6; l++)      //6 systematics per line
        {
          lstream >> up >> down;
          symmetriseErrors(up, down, &stmp, &dtmp);
          sysread[jdat][isys+l] = stmp;
          shift[jdat] += dtmp;
        }
        jdat++;
      }
    }
    idat+= ndatbin[i];
  }

  // Filter systematics, apply shift to data
  for (int i = 0; i < fNData; i++)
  {
    fSys[i][1].mult = sysread[i][0];    //Uncorrelated systematic
    fSys[i][1].type = MULT;
    fSys[i][1].name = "UNCORR";

    //Need to avoid including luminosity twice in systematics
    for (int j = 2; j < 7; j++)
      fSys[i][j].mult = sysread[i][j-1];       // First 5 systematics (%)

    for (int j = 7; j < fNSys; j++)
      fSys[i][j].mult = sysread[i][j];        // Remaining systematics (%)

    for (int j = 2; j < fNSys; j++)
    {
      fSys[i][j].type = MULT;
      fSys[i][j].name = "CORR";
    }

    fSys[i][0].mult = 6.1;              // Luminosity: 6.1%
    fSys[i][0].type = MULT;
    fSys[i][0].name = "D0LUMI";

    // Convert additive uncertainties
    fStat[i] *= fData[i]*1e-2;
    for (int j = 0; j < fNSys; j++)
      fSys[i][j].add = fSys[i][j].mult*fData[i]*1e-2;

    fData[i] *= (1.0 + shift[i]*0.01);       // Shift due to asymmetric uncertainty

  }

  f1.close();
  f2.close();
}


/**
 *
 *     Z rapidity distribution measurements at the Tevatron
 *     by the D0 collaboration
 *
 *     Published data from D0 can be found in the following reference:
 *
 *     Published data can be found in
 *     Subjects: 	High Energy Physics - Experiment (hep-ex)
 *     Journal reference: 	Phys.Rev.D76:012003,2007
 *     DOI: 	10.1103/PhysRevD.76.012003
 *     Report number: 	Fermilab-Pub-07-040-E
 *       Cite as: 	arXiv:hep-ex/0702025v1
 *
 *     which does not includes the full experimental systematic covariance
 *     matrix -> Systematic uncertainties are to be added in quadrature
 *     with the statistical uncertainties
 *     (Checked with the experimentalists)
 *
 *
 *      The measured spectra has been folder to positive rapidities
 *     due to symmetry
 *
 *     Experimental data is given as (1/sigma) dsigma/dy
 *     on different rapidity bins, fNData = 28
 *
 *     Electroweak convenor at D0: Heidi Schelmann, schellman@fnal.gov
 *
 */

void D0ZRAPFilter::ReadData()
{
  // Opening files
  fstream f1;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/D0-Zrapdist.data";
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

    lstream >> fKin1[i];   //Z boson rapidity y
    fKin2[i] = MZ2;        //Mass Z squared
    fKin3[i] = s;          //sqrt(s)

    lstream >> fData[i];
    lstream >> fStat[i];            //stat

    double up = 0;double down = 0;
    lstream >> up >> down;          // Asymmetric systematic

    double sys = 0;
    double shift = 0;
    symmetriseErrors(up,down,&sys,&shift);

    fSys[i][0].add = sys;
    fSys[i][0].mult = (sys*100.0)/fData[i];

    fSys[i][0].type = MULT;
    fSys[i][0].name = "UNCORR";     //treat sys as uncorrelated

    fData[i]+= shift;     // Shift due to asymmetric error
  }

  f1.close();

}

/**
 *
 *     W asymmetry rapidity distribution measurements with muon production
 *     at the Tevatron by the D0 collaboration.
 *
 *     Published data from D0 can be found in the following reference:
 *
 *     Published data can be found in
 *     Subjects: 	High Energy Physics - Experiment (hep-ex)
 *     Journal reference:	Phys. Rev. D 88, 091102 (2013)
 *     DOI:	10.1103/PhysRevD.88.091102
 *     Report number:	FERMILAB-PUB-13-361-E
 *     Cite as:	arXiv:1309.2591 [hep-ex]
 *
 *     The muon charge asymmetry is presented in the kinematic region where 
 *     the muon transverse momentum is 'pμT>25 GeV' and missing transverse 
 *     energy 'EmissingT>25 GeV'.
 *
 *     There is a total of 7 systematic uncertainties
 *     Description of systematic sources:
 *
 *     EW bkg
 *     MJ bkg 
 *     Chardge mis-id
 *     Relative charge efficiency
 *     Magnet polarity weighting
 *     Momentum resolution
 *     Trigger isolation (uncorrelated)
 *
 *     All uncertianties in the paper are mulitplied by 100.
 *     A total 7 systematics with 6 of which are correlated and 1 unccorelated.
 *
 */

void D0WMASYFilter::ReadData()
{
  // Opening files
  fstream f1;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/D0_wmAsym_25_dat.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Read data file
  string line;
  double eta1, eta2, etaAverage;
  double sysread;
  double s = 1960;
  for(int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> eta1 >> eta2;
    lstream >> etaAverage;
    fKin1[i] = etaAverage; // Average bin to bin rapidity
    
    fKin2[i] = MW * MW;    // Mass of the W-boson squared
    fKin3[i] = s;          // sqrt(s) Center of mass energy
    
    lstream >> fData[i];
    fData[i] /= 100.0;   // fData divided by 100
    
    lstream >> fStat[i];
    fStat[i] /= 100.0;   // fStat divided by 100
    
    for(int j = 0; j < fNSys; j++)
    {
      lstream >> sysread;
      fSys[i][j].add = sysread / 100.0;  // The systematics are divided by 100
      fSys[i][j].type = ADD;
      fSys[i][j].name = (j == 6) ? "UNCORR" : "CORR";
      fSys[i][j].mult = fSys[i][j].add * 1e2 / fData[i];

    }
  }
  f1.close();
}

/**
 *
 *     W asymmetry rapidity distribution measurements with electron production
 *     at the Tevatron by the D0 collaboration.
 *
 *     Published data from D0 can be found in the following reference:
 *
 *     Published data can be found in
 *     Subjects:	High Energy Physics - Experiment (hep-ex)
 *     Journal reference:	Phys. Rev. D 91, 032007 (2015)
 *     DOI:	10.1103/PhysRevD.91.032007
 *     Report number:	FERMILAB-PUB-14-514-E
 *     Cite as:	arXiv:1412.2862 [hep-ex]
 *
 *     The electron charge asymmetry is presented in the kinematic region where 
 *     the electron transverse momentum is 'peT>25 GeV' and missing transverse 
 *     energy 'EmissingT>25 GeV'.
 *
 *     There is a total of 9 systematic uncertainties
 *     Description of systematic sources:
 *
 *     Gen
 *     EMID      (uncorrelated)
 *     K_eff     (uncorrelated)
 *     Energy
 *     Recoil
 *     Model
 *     Bkgs
 *     Q_mis     (uncorrelated)
 *     Unfolding (uncorrelated)
 *
 *     All uncertianties in the raw-data are mulitplied by 100. (In the paper by 1000)
 *     A total 9 systematics with 5 of which are correlated and 4 unccorelated.
 *
 */

void D0WEASYFilter::ReadData()
{
  // Opening files
  fstream f1;
    
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/D0_Wel_pt25_asym.dat";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Read data file
  string line;
  double eta1, eta2, etaAverage;
  double sysread;
  double s = 1960;
  for(int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> eta1 >> eta2; // Redundant step
    lstream >> etaAverage;   
    fKin1[i] = etaAverage;   // Average bin to bin rapidity
    
    fKin2[i] = MW * MW;      // Mass of W-boson squared
    fKin3[i] = s;            // sqrt(s) Center of mass energy
    
    lstream >> fData[i];
    fData[i] /= 100.0;   // fData divided by 100
    
    lstream >> fStat[i];
    fStat[i] /= 100.0;   // fStat divided by 100
    
    for(int j = 0; j < fNSys; j++)
    {
        lstream >> sysread;
        fSys[i][j].add = sysread / 100.0;  // Systematics are divided by 100
        fSys[i][j].mult = fSys[i][j].add * 1e2 / fData[i];
        // If EMID, K_eff, Q_mis or Unfolding, treat as uncorrelated
	// In the paper they say that unfolding (j=8) is totally correlated
	// however in the HERApaper they say that exps say it must be treated as correlated
        fSys[i][j].name  = (j == 1 || j == 2 || j == 7 || j == 8) ? "UNCORR" :  "CORR";
        fSys[i][j].type = ADD;
    }
  }
  f1.close();
}


