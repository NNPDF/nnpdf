/**
 * heraf2charm.f
 *
 * Combined H1+ZEUS HERA I+II data on reduced charm production
 *     cross sections in NC DIS
 *
 *     Combined measurement of the charm production cross section in $e+-p$ scattering at HERA
 *
 *                               --  DESY 12-172 --
 *
 *     Data files contain:
 *
 *     d12-172.charm-ep.dat   -- Reduced cross sections of charm production at HERA
 *
 *     The bin numbers correspond to those in Table 4 of DESY 12-172
 *
 *     Q^2 values are given in GeV^2. x,y are Bjorken x and inelasticity. s_r stands for the reduced cross section.
 *
 *     Errors are quoted in % of the reduced cross sections.
 *
 *-----------------------------------------------------------
 *     stat           stands for the statistical uncertainty.
 *     unc            stands for the uncorrelated systematic uncertainty.
 *     tot            Tot stands for the total uncertainty of the measurement.
 *     d1-d42      are the correlated systematic experimental systematic error sources, they are correlated
 *               also across the data tables.
 *     p1-p2       uncertainties arising in the combination procedure (procedural uncertainties), they are
 *               correlated also across the data tables.
 *
 *
 */

#include "HERA2-C.h"

void HERAF2CHARMFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/d12-172.charm-ep.dat";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  int index = 0;
  
  double datain[fNData][52];
  double syscor[fNData][fNSys];
  double absnorexp = 0.0;
  
  for (int i = 0; i < fNData; i++)
  {
    for (int l = 0; l < fNSys; l++)
      syscor[i][l] = 0.0;
  }

  // Reading data
  index = 0;
  string tmp;
  getline(f1,tmp);
  for (int i = 0; i < fNData; i++)
  {
    for (int j = 0; j < 52; j++)
      f1 >> datain[i][j];
    if (datain[i][0] != i+1) {
      cerr << "Mismatch" << endl;
      exit(-1);
    }
    index++;
  }
  
  if (index != fNData)
  {
    cerr << "Mismatch in the number of data points" << endl;
    exit(-1);
  }
  
  // Filtering data
  double stat[fNData], x[fNData], y[fNData];
  double q2[fNData], xsec_charm[fNData];
  
  for (int i = 0; i < fNData; i++)
  {
    x[i]   = datain[i][2];
    y[i]   = datain[i][3];
    q2[i]  = datain[i][1];
    
    // Observable
    xsec_charm[i] = datain[i][4];
    
    // Statistical errors - percentage with respect the observable
    stat[i] = datain[i][5];
    
    // Uncorrelated systematics
    fSys[i][1].mult = datain[i][6];
    fSys[i][1].type = ADD;
    fSys[i][1].name = "UNCORR";
    
    // Systematic errors
    double sys_tot = 0.0;
    for (int l = 0; l < fNSys-2; l++)
    {
      syscor[i][l] = datain[i][l+8];
      sys_tot += syscor[i][l]*syscor[i][l];
    }
    sys_tot = sqrt(sys_tot);
    
    // Normalization uncertainty
    fSys[i][0].mult = absnorexp;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CORR";
    
    // Check total error
    double tot = datain[i][7];
    double tot_tmp = sqrt(sys_tot*sys_tot + stat[i]*stat[i] + fSys[i][1].mult*fSys[i][1].mult);
    if (fabs(tot-tot_tmp) > 1.4)
    {
      cerr << "Mismatch " << i << "\t" << tot << "\t" << tot_tmp << endl;
      exit(-1);
    }
    
    fKin1[i] = x[i];
    fKin2[i] = q2[i];
    fKin3[i] = y[i];
    
    fData[i] = xsec_charm[i];
    fStat[i] = stat[i]*fData[i]*1e-2;
    
    for (int l = 2; l < fNSys; l++)
    {
      fSys[i][l].mult = syscor[i][l];
      fSys[i][l].type = ADD;
      fSys[i][l].name = "CORR";
    }
    
    for (int l = 0; l < fNSys; l++)
      fSys[i][l].add = fSys[i][l].mult*fData[i]*1e-2;
      
  }
  
  f1.close();
}

