/****************************************************************************
*
*     h1-hera2.f
*
*     H1 HERA-II inclusive NC and CC structure function data
*
*     Four datasets
*
*     H1-HERA2-CCep.data    29 datapoints
*     H1-HERA2-NCep.data    136 datapoints
*     H1-HERA2-CCem.data    29 datapoints
*     H1-HERA2-NCem.data    139 datapoints
*
*     where the data with positive and negative polarity has been combined
*
*     The format of all the NC files is the same
*
*  Q2  x σred  sigtot(%) sigstat(%) siguncorr(%) siguncorrE(%) siguncorrN(%)
*     -> sigcorr(%)  (sigcorr_i(%), i=1,5)
*
*     For the  total uncorrelated systematic errors, two of its contributions
*     from the electron energy error  and the hadronic energy error
*     are also shown,
*     The effect of the other uncorrelated
*     systematic errors is included in siguncorr
*     In addition the correlated systematic (δcor) and its contributions
*     from a positive variation of one standard deviation of the five
*     different sources
*
*     Note that the inelasticity can be taken from E, Q and x
*
*     Note that normalization and polarization errors are not included
*
*     Instead the format for all the CC files is
*
*     Q2  x y σred  sigtot(%) sigstat(%) siguncorr(%) siguncorrh(%)
*     -> sigcorr(%)  (sigcorr_i(%), i=1,4)
*
*     Check inelasticiy with the numerical formula
*
*     On the other hand, the list of systematics is different in the two cases
*
*     NC: E+, theta+, h+, N+, B+
*
*     CC: V+, h+, N+, B+
*
*     so in total we have 6 systematics - > Check correct correlation
*
****************************************************************************/

#include "H1HERA2.h"

void H1HERA2NCEPFilter::ReadData()
{
  fstream f1;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
           << fSetName << "/d12-107.table26.txt"; //"/H1-HERA2-NCep.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  double stat, stattot, sysuncor, systot, proc1, proc2;

  for (int i = 0; i < fNData; i++)
    {
      f1 >> fKin2[i] >> fKin1[i] >> fData[i]
         >> stattot >> stat >> sysuncor >> proc1 >> proc2 >> systot;

      // inelasticity y
      fKin3[i] = fKin2[i]/fKin1[i]/pow(319.0,2.0);

      // statistical uncertainties
      fStat[i] = stat*fData[i]*1e-2;

      // "uncorrelated" systematics are correlated
      fSys[i][0].mult = proc1;
      fSys[i][0].type = ADD;
      fSys[i][0].name = "CORR";

      fSys[i][1].mult = proc2;
      fSys[i][1].type = ADD;
      fSys[i][1].name = "CORR";

      fSys[i][2].mult = sqrt(sysuncor*sysuncor - proc1*proc1 - proc2*proc2);
      fSys[i][2].type = ADD;
      fSys[i][2].name = "CORR";

      fSys[i][3].mult = 2.3; // absolute 2.3% luminosity uncertainty
      fSys[i][3].type = MULT;
      fSys[i][3].name = "H1HERA2LUMI";

      fSys[i][4].mult = 1.5; // independent 1.5% luminosity uncertainty
      fSys[i][4].type = MULT;
      fSys[i][4].name = "CORR";

      // correlated systematics
      //delta E:
      f1 >> fSys[i][5].mult;
      fSys[i][5].type = ADD;
      fSys[i][5].name = "H1HERA2SYSE";

      //delta theta:
      f1 >> fSys[i][6].mult;
      fSys[i][6].type = ADD;
      fSys[i][6].name = "H1HERA2SYSTheta";

      //delta h:
      f1 >> fSys[i][7].mult;
      fSys[i][7].type = ADD;
      fSys[i][7].name = "H1HERA2SYSh";

      //delta N:
      f1 >> fSys[i][8].mult;
      fSys[i][8].type = ADD;
      fSys[i][8].name = "H1HERA2SYSN_NC";

      //delta B:
      f1 >> fSys[i][9].mult;
      fSys[i][9].type = ADD;
      fSys[i][9].name = "H1HERA2SYSB";

      //extra corsys uncertainty
      fSys[i][10].mult = systot*systot;
      for (int j = 5; j < 10; j++) fSys[i][10].mult -= fSys[i][j].mult*fSys[i][j].mult;
      fSys[i][10].mult = sqrt(max(fSys[i][10].mult,0.0));
      fSys[i][10].type = ADD;
      fSys[i][10].name = "CORR";

      for (int j = 0; j < fNSys; j++)
        fSys[i][j].add = fSys[i][j].mult*fData[i]*1e-2;

      // check of filtering data
      double cor_tot = 0;
      for (int j = 5; j < fNSys; j++)
        cor_tot += fSys[i][j].mult*fSys[i][j].mult;

      cor_tot = sqrt(cor_tot);

      // check sum of correlated systematics is total sys
      if (fabs(cor_tot - systot) > 2e-2)
        {
          cout << "Mismatch in filtering H1 HERA-II data" << endl;
          cout << "ii, cor_tot, tot_sys" << endl;
          cout << i << "\t" << cor_tot << "\t" << systot << endl;
          exit(-1);
        }

      // check total uncertainty adds up
      double exp_quad = sqrt(stat*stat + systot*systot + sysuncor*sysuncor);
      if (fabs(stattot - exp_quad) > 2e-2)
        {
          cout << "Mismatch in filtering H1 HERA-II data" << endl;
          cout << "tot  exp_quad" << endl;
          cout << i << "\t" << stattot << "\t" << exp_quad << endl;
          exit(-1);
        }

      // check uncorrelated uncertainties are not larger than total uncorrelated
      double uncor_tot = sqrt(proc1*proc1 + proc2*proc2);
      if (uncor_tot > sysuncor)
        {
          cout << "Mismatch in filtering H1 HERA-II data" << endl;
          cout << "uncor_tot, sysuncor = " << endl;
          cout << uncor_tot << "\t" << sysuncor << endl;
          exit(-1);
        }
    }

  f1.close();
}



void H1HERA2NCEMFilter::ReadData()
{
  fstream f1;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
           << fSetName << "/d12-107.table25.txt";//"/H1-HERA2-NCem.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  double stat, stattot, sysuncor, systot, proc1, proc2;

  for (int i = 0; i < fNData; i++)
    {
      f1 >> fKin2[i] >> fKin1[i] >> fData[i]
         >> stattot >> stat >> sysuncor >> proc1 >> proc2 >> systot;

      // inelasticity y
      fKin3[i] = fKin2[i]/fKin1[i]/pow(319.0,2.0);

      // statistical uncertainties
      fStat[i] = stat*fData[i]*1e-2;

      // "uncorrelated" systematics are correlated
      fSys[i][0].mult = proc1;
      fSys[i][0].type = ADD;
      fSys[i][0].name = "CORR";

      fSys[i][1].mult = proc2;
      fSys[i][1].type = ADD;
      fSys[i][1].name = "CORR";

      fSys[i][2].mult = sqrt(sysuncor*sysuncor - proc1*proc1 - proc2*proc2);
      fSys[i][2].type = ADD;
      fSys[i][2].name = "CORR";

      fSys[i][3].mult = 2.3; // absolute 2.3% luminosity uncertainty
      fSys[i][3].type = MULT;
      fSys[i][3].name = "H1HERA2LUMI";

      fSys[i][4].mult = 1.5; // independent 1.5% luminosity uncertainty
      fSys[i][4].type = MULT;
      fSys[i][4].name = "CORR";

      // correlated systematics
      //delta E:
      f1 >> fSys[i][5].mult;
      fSys[i][5].type = ADD;
      fSys[i][5].name = "H1HERA2SYSE";

      //delta theta:
      f1 >> fSys[i][6].mult;
      fSys[i][6].type = ADD;
      fSys[i][6].name = "H1HERA2SYSTheta";

      //delta h:
      f1 >> fSys[i][7].mult;
      fSys[i][7].type = ADD;
      fSys[i][7].name = "H1HERA2SYSh";

      //delta N:
      f1 >> fSys[i][8].mult;
      fSys[i][8].type = ADD;
      fSys[i][8].name = "H1HERA2SYSN_NC";

      //delta B:
      f1 >> fSys[i][9].mult;
      fSys[i][9].type = ADD;
      fSys[i][9].name = "H1HERA2SYSB";

      //extra corsys uncertainty
      fSys[i][10].mult = systot*systot;
      for (int j = 5; j < 10; j++) fSys[i][10].mult -= fSys[i][j].mult*fSys[i][j].mult;
      fSys[i][10].mult = sqrt(max(fSys[i][10].mult,0.0));
      fSys[i][10].type = ADD;
      fSys[i][10].name = "CORR";

      for (int j = 0; j < fNSys; j++)
        fSys[i][j].add = fSys[i][j].mult*fData[i]*1e-2;

      // check of filtering data
      double cor_tot = 0;
      for (int j = 5; j < fNSys; j++)
        cor_tot += fSys[i][j].mult*fSys[i][j].mult;

      cor_tot = sqrt(cor_tot);

      // check sum of correlated systematics is total sys
      if (fabs(cor_tot - systot) > 2e-2)
        {
          cout << "Mismatch in filtering H1 HERA-II data" << endl;
          cout << "ii, cor_tot, tot_sys" << endl;
          cout << i << "\t" << cor_tot << "\t" << systot << endl;
          exit(-1);
        }

      // check total uncertainty adds up
      double exp_quad = sqrt(stat*stat + systot*systot + sysuncor*sysuncor);
      if (fabs(stattot - exp_quad) > 2e-2)
        {
          cout << "Mismatch in filtering H1 HERA-II data" << endl;
          cout << "tot  exp_quad" << endl;
          cout << i << "\t" << stattot << "\t" << exp_quad << endl;
          exit(-1);
        }

      // check uncorrelated uncertainties are not larger than total uncorrelated
      double uncor_tot = sqrt(proc1*proc1 + proc2*proc2);
      if (uncor_tot > sysuncor)
        {
          cout << "Mismatch in filtering H1 HERA-II data" << endl;
          cout << "uncor_tot, sysuncor = " << endl;
          cout << uncor_tot << "\t" << sysuncor << endl;
          exit(-1);
        }
    }

  f1.close();
}

  // The raw data for the HERA-II H1 measurement
  // Extracted from Tables 27 and 28 of arxiv:1206.7007
  // are given in terms of the double differential cross section
  // in units of pb/GeV2
  // Need to define a proper conversion factor to write them
  // in terms of the reduced cross section, which is what
  // fkgenerator outputs
  double const Gfermi = 0.0000116637; // GeV^{-2}
  //double const MW = 80.385; // GeV, PDG 2012
  double const conv_pb_gev = 0.389*pow(10.0,9.0); // GeV^2 * pb in natural units (dimensionless)

void H1HERA2CCEPFilter::ReadData()
{
  fstream f1;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
           << fSetName << "/H1-HERA2-CCep.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  double stat, stattot, sysuncor, systot, proc1;

  for (int i = 0; i < fNData; i++)
    {
      f1 >> fKin2[i] >> fKin1[i] >> fKin3[i] >> fData[i]
         >> stattot >> stat >> sysuncor >> proc1 >> systot;

  // Conversion factir from CC 2D xsec -> CC reduced xsec
      double conv_fact = 4.0 * M_PI * fKin1[i] / ( Gfermi * Gfermi );
      conv_fact *= pow( ( ( MW * MW + fKin2[i] )/ ( MW * MW )),2.0);
      conv_fact /= conv_pb_gev;
      // Extra factor 1/2 from the defnition of reduced cross sections
      // in fkgenerator
      conv_fact /= 2.0;
      // Now correct the data
      fData[i]*=conv_fact;

      // statistical uncertainties
      fStat[i] = stat*fData[i]*1e-2;

      // "uncorrelated" systematics are correlated
      fSys[i][0].mult = proc1;
      fSys[i][0].type = ADD;
      fSys[i][0].name = "CORR";

      fSys[i][1].mult = sqrt(sysuncor*sysuncor - proc1*proc1);
      fSys[i][1].type = ADD;
      fSys[i][1].name = "CORR";

      fSys[i][2].mult = 2.3; // absolute 2.3% luminosity uncertainty
      fSys[i][2].type = MULT;
      fSys[i][2].name = "H1HERA2LUMI";

      fSys[i][3].mult = 1.5; // independent 1.5% luminosity uncertainty
      fSys[i][3].type = MULT;
      fSys[i][3].name = "CORR";

      // correlated systematics
      //delta V:
      f1 >> fSys[i][4].mult;
      fSys[i][4].type = ADD;
      fSys[i][4].name = "H1HERA2SYSV";

      //delta h:
      f1 >> fSys[i][5].mult;
      fSys[i][5].type = ADD;
      fSys[i][5].name = "H1HERA2SYSh";

      //delta N:
      f1 >> fSys[i][6].mult;
      fSys[i][6].type = ADD;
      fSys[i][6].name = "H1HERA2SYSN_CC";

      //delta B:
      f1 >> fSys[i][7].mult;
      fSys[i][7].type = ADD;
      fSys[i][7].name = "H1HERA2SYSB";

      //extra corsys uncertainty
      fSys[i][8].mult = systot*systot;
      for (int j = 4; j < 8; j++) fSys[i][8].mult -= fSys[i][j].mult*fSys[i][j].mult;
      fSys[i][8].mult = sqrt(max(fSys[i][8].mult,0.0));
      fSys[i][8].type = ADD;
      fSys[i][8].name = "CORR";

      for (int j = 0; j < fNSys; j++)
        fSys[i][j].add = fSys[i][j].mult*fData[i]*1e-2;

      // check of filtering data
      double cor_tot = 0;
      for (int j = 4; j < fNSys; j++)
        cor_tot += fSys[i][j].mult*fSys[i][j].mult;

      cor_tot = sqrt(cor_tot);

      // check sum of correlated systematics is total sys
      if (fabs(cor_tot - systot)/systot > 4e-1)
        {
          cout << "Mismatch in filtering H1 HERA-II data" << endl;
          cout << "ii, cor_tot, tot_sys" << endl;
          cout << i << "\t" << cor_tot << "\t" << systot << endl;
          exit(-1);
        }

      // check total uncertainty adds up
      double exp_quad = sqrt(stat*stat + systot*systot + sysuncor*sysuncor);
      if (fabs(stattot - exp_quad)/exp_quad > 5e-2)
        {
          cout << "Mismatch in filtering H1 HERA-II data" << endl;
          cout << "tot  exp_quad" << endl;
          cout << i << "\t" << stattot << "\t" << exp_quad << endl;
          exit(-1);
        }

      // check uncorrelated uncertainties are not larger than total uncorrelated
      if (proc1 > sysuncor)
        {
          cout << "Mismatch in filtering H1 HERA-II data" << endl;
          cout << "proc1, sysuncor = " << endl;
          cout << proc1 << "\t" << sysuncor << endl;
          exit(-1);
        }
    }

  f1.close();
}


void H1HERA2CCEMFilter::ReadData()
{
  fstream f1;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
           << fSetName << "/H1-HERA2-CCem.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  double stat, stattot, sysuncor, systot, proc1;

  for (int i = 0; i < fNData; i++)
    {
      f1 >> fKin2[i] >> fKin1[i] >> fKin3[i] >> fData[i]
         >> stattot >> stat >> sysuncor >> proc1 >> systot;

      // Conversion factir from CC 2D xsec -> CC reduced xsec
      double conv_fact = 4.0 * M_PI * fKin1[i] / ( Gfermi * Gfermi );
      conv_fact *= pow( ( ( MW * MW + fKin2[i] )/ ( MW * MW )),2.0);
      conv_fact /= conv_pb_gev;
      // Extra factor 1/2 from the defnition of reduced cross sections
      // in fkgenerator
      conv_fact /= 2.0;
      // Now correct the data
      fData[i]*=conv_fact;

      // statistical uncertainties
      fStat[i] = stat*fData[i]*1e-2;

      // "uncorrelated" systematics are correlated
      fSys[i][0].mult = proc1;
      fSys[i][0].type = ADD;
      fSys[i][0].name = "CORR";

      fSys[i][1].mult = sqrt(sysuncor*sysuncor - proc1*proc1);
      fSys[i][1].type = ADD;
      fSys[i][1].name = "CORR";

      fSys[i][2].mult = 2.3; // absolute 2.3% luminosity uncertainty
      fSys[i][2].type = MULT;
      fSys[i][2].name = "H1HERA2LUMI";

      fSys[i][3].mult = 1.5; // independent 1.5% luminosity uncertainty
      fSys[i][3].type = MULT;
      fSys[i][3].name = "CORR";

      // correlated systematics
      //delta V:
      f1 >> fSys[i][4].mult;
      fSys[i][4].type = ADD;
      fSys[i][4].name = "H1HERA2SYSV";

      //delta h:
      f1 >> fSys[i][5].mult;
      fSys[i][5].type = ADD;
      fSys[i][5].name = "H1HERA2SYSh";

      //delta N:
      f1 >> fSys[i][6].mult;
      fSys[i][6].type = ADD;
      fSys[i][6].name = "H1HERA2SYSN_CC";

      //delta B:
      f1 >> fSys[i][7].mult;
      fSys[i][7].type = ADD;
      fSys[i][7].name = "H1HERA2SYSB";

      //extra corsys uncertainty
      fSys[i][8].mult = systot*systot;
      for (int j = 4; j < 8; j++) fSys[i][8].mult -= fSys[i][j].mult*fSys[i][j].mult;
      fSys[i][8].mult = sqrt(max(fSys[i][8].mult,0.0));
      fSys[i][8].type = ADD;
      fSys[i][8].name = "CORR";

      for (int j = 0; j < fNSys; j++)
        fSys[i][j].add = fSys[i][j].mult*fData[i]*1e-2;

      // check of filtering data
      double cor_tot = 0;
      for (int j = 4; j < fNSys; j++)
        cor_tot += fSys[i][j].mult*fSys[i][j].mult;

      cor_tot = sqrt(cor_tot);

      // check sum of correlated systematics is total sys
      if (fabs(cor_tot - systot)/systot > 4e-1)
        {
          cout << "Mismatch in filtering H1 HERA-II data" << endl;
          cout << "ii, cor_tot, tot_sys" << endl;
          cout << i << "\t" << cor_tot << "\t" << systot << endl;
          exit(-1);
        }

      // check total uncertainty adds up
      double exp_quad = sqrt(stat*stat + systot*systot + sysuncor*sysuncor);
      if (fabs(stattot - exp_quad)/exp_quad > 5e-2)
        {
          cout << "Mismatch in filtering H1 HERA-II data" << endl;
          cout << "tot  exp_quad" << endl;
          cout << i << "\t" << stattot << "\t" << exp_quad << endl;
          exit(-1);
        }

      // check uncorrelated uncertainties are not larger than total uncorrelated
      if (proc1 > sysuncor)
        {
          cout << "Mismatch in filtering H1 HERA-II data" << endl;
          cout << "proc1, sysuncor = " << endl;
          cout << proc1 << "\t" << sysuncor << endl;
          exit(-1);
        }
    }

  f1.close();
}



//////////////////////////////////////////////////////////////////
// Now the HERA-II H1 low Q2 data
//////////////////////////////////////////////////////////////////

/*

The H1 HERA-II low-Q2 data has been obtained from
http://www-h1.desy.de/psfiles/figures/d10-228.syst_460_575.txt
The data is for reduced NC DIS cross sections
Since the data is well below the MZ threshold, we can use for the theory computations either to assume electrons or to assume positrons

Measurement of the Inclusive e{\pm}p Scattering Cross Section at High Ineelasticity y and of the Structure Function $F_L$.
 F.D. Aaron et al. DESY-10-228, Dec 2010. 71pp.
 Published in Eur.Phys.J.C71:1579,2011.
 e-Print: arXiv:1012.4355 [hep-ex]

d10-228.syst_460_575.dat  -- combination of the Ep=460 GeV and Ep=575 GeV HERA-II data

Q^2 values are given in GeV^2. x,y are Bjorken x and inelasticity. Sr stands for the reduced cross section.
F2 stands for the structure function F2. CME stands for centre-of-mass energy.
Errors are quoted in % of the reduced cross sections.

stat           stands for the statistical uncertainty.
unc            stands for the uncorrelated systematic uncertainty.
corr           is the sum in quadrature of sysL1-sysL8 for 460-575
tot            is the sum of stat, unc and cor in quadrature
SysL1-SysL8    are the correlated systematic ystematic error sources.

 */


void H1HERA2LOWQ2Filter::ReadData()
{

  fstream f1;
  stringstream datafile("");
  // Read war data including all the list of correlated systematics
  datafile << dataPath() << "rawdata/"
           << fSetName << "/d10-228.syst_460_575.txt";
  f1.open(datafile.str().c_str(), ios::in);

  // Check that the file can the properly opened
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Note that fNSys includes the uncorrelated systematic errors
  // and the luminosity uncertainty
  // So in this case we have fNSys = 8 + 1 + 1

  for (int i = 0; i < fNData; i++)
    {

      double adum=0, stat=0, sysuncor=0, systot=0, ertot=0;
      int ibin=0;
      int const nsys_cor=8;
      double syscor[nsys_cor]={0.0};
      // The format of the data is:
      //  Bin    Q2  X   y     Sr F2 Stat Unc Corr   Tot  SysL1  ... SysL8  CME
      f1 >> ibin>>fKin2[i] >> fKin1[i] >>fKin3[i] >> fData[i]>>adum>>
	stat >> sysuncor >> systot >> ertot;
      for(int isys=0;isys<nsys_cor;isys++){
	f1>>syscor[isys];
	//cout<<isys<<"   "<<syscor[isys];
      }
      //cout<<endl;
      f1>> adum; // Reading the CM energy

      // Check bin ordering
      if(ibin!=(i+1)){ cerr << "Mismatch" << endl;exit(-1);}

      // statistical uncertainties
      // Uncorrelated systematics are included below, different treatment
      // Saved in abolsute value, not in percentage as default in the raw data
      fStat[i] = stat*fData[i]*1e-2;

      // uncorrelated systematics
      // Saved in percent, to be able to use T0 afterwards
      fSys[i][0].type = ADD;
      fSys[i][0].name = "UNCORR";
      fSys[i][0].mult = sysuncor; // Save as percentage
      fSys[i][0].add = fSys[i][0].mult*fData[i]*1e-2; // Save in absolute

      // Normalization uncertainty
      // Uncorrelated to the luminosity uncertainty of the other H1 HERA-II data
      // See Sect. 3.7 of 1012.4355v1.pdf
      fSys[i][1].mult = 4.0; // absolute luminosity uncertainty, in percent
      fSys[i][1].type = MULT;
      fSys[i][1].name = "H1HERA2LOWQ2LUMI"; // Save as percentage
      fSys[i][1].add = fSys[i][1].mult*fData[i]*1e-2; // Save in absolute

      // correlated systematics
      for (int j = 2; j < fNSys; j++)
        {

	  // Define type of systematic error
          fSys[i][j].type = ADD;
          fSys[i][j].name = "CORR";

	  // Saved as multiplicative - percentage value
          fSys[i][j].mult = syscor[j-2] ; // offset between fSys and syscor

	  // Saved as additive  - absolute values
	  fSys[i][j].add = fSys[i][j].mult*fData[i]*1e-2;

        }

      //-----------------------
      // checks of filtering data

      // Check that the total correlated systematic error is the
      // sum in quadrature of the eight sources of systematic error

      double cor_tot = 0;
      for (int j = 2; j < fNSys; j++){
	//	cout<<fSys[i][j].mult<<endl;
        cor_tot += fSys[i][j].mult*fSys[i][j].mult;
      }
      cor_tot = sqrt(cor_tot);

      // due to rounding errors, this checked performed with 1% accuracy
      double diff = 1e2*fabs( (cor_tot - systot) /systot) ;
      double const diff_max=1.0;
      //   cout<<cor_tot<<" "<<systot<<endl;
      //cout<<"diff = "<<diff<<endl;
      if ( diff > diff_max)
        {
          cout << "Mismatch in filtering H1 low Q2 HERA-II data" << endl;
          cout << "idat, cor_tot, tot_sys, diff" << endl;
          cout << i << "\t" << cor_tot << "\t" << systot <<"\t"<<diff<< endl;
          exit(-1);
        }

      // Now check that the total error is the sum of
      // i) Correlated systematics
      // ii) Uncorrelated systematics
      // iii) Statistical uncertainties
      double exp_tot = sqrt(stat*stat + cor_tot*cor_tot + fSys[i][0].mult*fSys[i][0].mult);
      double diff_tot = 1e2*fabs( (exp_tot - ertot) /ertot ) ;
      double const diff_tot_max=1; // Checked with 1% accuracy due to rounding errors
      if (diff_tot > diff_tot_max )
        {
          cout << "Mismatch in filtering H1 low-Q2 HERA-II data" << endl;
          cout << "idat   ertot  exp_tot diff" << endl;
          cout << i << "\t" <<ertot << "\t" << exp_tot<<"\t"<<diff << endl;
          exit(-1);
        }


    } // End loop over data points

  f1.close();
}


//////////////////////////////////////////////////////////////////
// Now the HERA-II H1 high-y data
//////////////////////////////////////////////////////////////////

/*

The high-y data is taken from the HERA-II H1 publication

Measurement of the Inclusive e{\pm}p Scattering Cross Section at High Inelasticity y and of the Structure Function FL

Comments:	70 pages, 29 figures, 23 tables
Subjects:	High Energy Physics - Experiment (hep-ex)
Journal reference:	Eur.Phys.J.C71:1579,2011
DOI:	10.1140/epjc/s10052-011-1579-4
Report number:	DESY 10-228

And the recommendation from the H1 management is to add Tables 10
and 11 of this paper.

Since the papers are not available in .txt format, they have
been extracted from the Latex source of the paper

These tables provide the NC reduced cross sections (since it is at low Q2,
no difference between electrons and positros). Uncertainties are quoted
in percent. There is on top of that a global normalisation uncertainty of 3%

Threre are 5 sources of correlated systematic errors, which
are common in table 10 and Table 11


 */

void H1HERA2HGHYFilter::ReadData()
{

  fstream f1;
  stringstream datafile("");
  // Read war data including all the list of correlated systematics
  datafile << dataPath() << "rawdata/"
           << fSetName << "/h1_hera2_highy_raw.data";
  f1.open(datafile.str().c_str(), ios::in);

  // Check that the file can the properly opened
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Note that fNSys includes the uncorrelated systematic errors
  // and the luminosity uncertainty
  // So in this case we have fNSys = 5 + 1 + 1

  for (int i = 0; i < fNData; i++)
    {
      double stat=0, sysuncor=0,  ertot=0;
      int const nsys_cor=5;
      double syscor[nsys_cor]={0.0};
      // The format of the data is:
      // Q2   x   y   sigma_red  stat unc tot  sys1 ... sys5
      f1 >>fKin2[i] >> fKin1[i] >>fKin3[i] >> fData[i]>> stat >> sysuncor >> ertot;
      for(int isys=0;isys<nsys_cor;isys++)
          f1>>syscor[isys];

      // statistical uncertainties
      // Uncorrelated systematics are included below, different treatment
      // Saved in abolsute value, not in percentage as default in the raw data
      fStat[i] = stat*fData[i]*1e-2;

      // uncorrelated systematics
      // Saved in percent, to be able to use T0 afterwards
      fSys[i][0].type = ADD;
      fSys[i][0].name = "UNCORR";
      fSys[i][0].mult = sysuncor; // Save as percentage
      fSys[i][0].add = fSys[i][0].mult*fData[i]*1e-2; // Save in absolute

      // Normalization uncertainty
      // Uncorrelated to the luminosity uncertainty of the other H1 HERA-II data
      // caption of Table 10 in arxiv:1012.4355
      fSys[i][1].mult = 3.0; // absolute luminosity uncertainty, in percent
      fSys[i][1].type = MULT;
      fSys[i][1].name = "H1HERA2HGHYLUMI"; // Save as percentage
      fSys[i][1].add = fSys[i][1].mult*fData[i]*1e-2; // Save in absolute

      // correlated systematics
      for (int j = 2; j < fNSys; j++)
        {

	  // Define type of systematic error
          fSys[i][j].type = ADD;
          fSys[i][j].name = "CORR";

	  // Saved as multiplicative - percentage value
          fSys[i][j].mult = syscor[j-2] ; // offset between fSys and syscor

	  // Saved as additive  - absolute values
	  fSys[i][j].add = fSys[i][j].mult*fData[i]*1e-2;

        }

      //-----------------------
      // checks of filtering data

      // Compute total correlated systematic error
      double cor_tot = 0;
      for (int j = 2; j < fNSys; j++){
	//	cout<<fSys[i][j].mult<<endl;
        cor_tot += fSys[i][j].mult*fSys[i][j].mult;
      }
      cor_tot = sqrt(cor_tot);

      // Now check that the total error is the sum of
      // i) Correlated systematics
      // ii) Uncorrelated systematics
      // iii) Statistical uncertainties
      double exp_tot = sqrt(stat*stat + cor_tot*cor_tot + fSys[i][0].mult*fSys[i][0].mult);
      double diff_tot = 1e2*fabs( (exp_tot - ertot) /ertot ) ;
      double const diff_tot_max=0.5; // Checked with 0.5% accuracy due to rounding errors
      if (diff_tot > diff_tot_max )
        {
          cout << "Mismatch in filtering H1 high-y HERA-II data" << endl;
          cout << "idat   ertot  exp_tot diff_tot" << endl;
          cout << i << "\t" <<ertot << "\t" << exp_tot<<"\t"<<diff_tot << endl;
          exit(-1);
        }


    } // End loop over data points

  f1.close();
}

/////////////////////////////////////////////////////
