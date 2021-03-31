#include "D0.h"

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

// Fix inconsistency between add and mult columns
void D0ZRAP_40Filter::ReadData()
{
  // Opening files
  fstream f1;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/D0ZRAP_SF/D0-Zrapdist.data";
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
    fData[i]+= shift; // Shift due to asymmetric error (before calculating mult)
    fSys[i][0].mult = (sys*100.0)/fData[i];

    fSys[i][0].type = MULT;
    fSys[i][0].name = "UNCORR";     //treat sys as uncorrelated

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


