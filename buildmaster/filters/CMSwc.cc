/**
 * Check data and writes it in a common format
 *
 * This file implements in buildmaster the CMS 7 TeV 5 1/fb
 * data on W production in association with charm quarks
 * arXiv:1310.1138
 *
 * There are two uncorrelated datasets
 * - Absolute differential cross sections for W + charm production
 * (Need to add the two FK tables for WplusC and WminusC production
 * to compare with the experimental data)
 * - Cross section ratio data, WplusC/Wminuc
 * (Need to take the ratio of the corresponding FK tables)
 *
 * Systematici uncertainties are supplemented by a theoretical uncertainty
 * that takes into account missing NNLO corrections in the matrix element.
 * This uncertainty was estimated by means as the asymmetric envelope of the 
 * 3pt renormalisation scale variation (Eq. 4.16 in 1906.10698).
*/

#include "CMSwc.h"

void wc_covmat(vector<double> &);

void CMSWCHARMTOTFilter::ReadData()
{
  // Opening files
  fstream f1;
  fstream f2;

  // Read raw data
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/CMSWCHARMTOT.data";
  f1.open(datafile.str().c_str(), ios::in);
  // Check file properly open
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream theoryfile("");
  theoryfile << dataPath() << "rawdata/"
	     << fSetName << "/theory.txt";
  f2.open(theoryfile.str().c_str(), ios::in);
  // Check file properly open
  if (f2.fail()) {
    cerr << "Error opening data file " << theoryfile.str() << endl;
    exit(-1);
  }

  // Now generate the experimental covariance matrix
  // Correlation matrix provided by the analyzers
  // No separation between systematic and statistical uncertainties provided
  vector<double> covmatv;
  wc_covmat(covmatv);
  int const ndat=5;
  if(ndat!=fNData)exit(-10);
  //double covmat[ndat][ndat]={{0.0}};
  int k=0;

  double** covmat = new double*[fNData];
  for(int i = 0; i < fNData; i++)
  {
    covmat[i] = new double[fNData];
  }

  for (int i=0; i<fNData; i++) {
    for (int j=0; j<fNData; j++){
      covmat[i][j] = covmatv.at(k);
      k++;
    }
  }

  /*
  // Check the covariance matrix
  for (int i=0; i<fNData; i++) {
       for (int j=0; j<fNData; j++)  cout<<covmat[i][j]<<"  ";
    cout<<endl;
  }
  */

  // Now read the data
  string line;
  double etamin, etamax;
  double MW2 = pow(MW,2.0);
  // Center of mass energy
  double s = 7000;

  // Reading data
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);

    lstream >> etamin >> etamax;
    fKin1[i] = (etamax + etamin)*0.5;    // eta
    fKin2[i] = MW2;                      // Mass W squared
    fKin3[i] = s;                        // sqrt(s)

    lstream >> fData[i];
    fData[i] *= 1e3;                    // correct for units to match FK tables

    lstream >> fStat[i];
    fStat[i] *= 1e3;                    // correct for units to match FK tables

    double sysdum=0.0;
    lstream >> sysdum;

  }

  // correct for units to match FK tables
  for(int i = 0; i < fNData; i++){
    for(int j = 0; j < fNData; j++){
      covmat[i][j] *= 1e6;
    }
  }

  // Generating artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
   {
     cerr << " in " << fSetName << endl;
     exit(-1);
   }

  for (int i = 0; i < fNData; i++)
    for (int l = 0; l < fNSys-2; l++)
    {
      fSys[i][l].add = syscor[i][l];
      fSys[i][l].mult = fSys[i][l].add*100/fData[i];
      fSys[i][l].type = MULT;
      fSys[i][l].name = "CORR";
    }

  for(int i = 0 ; i<fNData; i++)
    {
      getline(f2,line);
      istringstream kstream(line);
      
      kstream >> fSys[i][5].mult >> fSys[i][6].mult;
      fSys[i][5].mult /= sqrt(2.);
      fSys[i][5].add = fSys[i][5].mult*fData[i]/100;
      fSys[i][5].type = MULT;
      fSys[i][5].name = "SKIP";
      
      fSys[i][6].mult /= sqrt(2.);
      fSys[i][6].add = fSys[i][6].mult*fData[i]/100;
      fSys[i][6].type = MULT;
      fSys[i][6].name = "SKIP";
    }
  
  f1.close();
  f2.close();
  
  for(int i = 0; i < fNData; i++) 
    {
      delete[] syscor[i];
      delete[] covmat[i];
    }
  delete[] syscor;
  delete[] covmat;
  
}

///////////////////////////

// Computation of the experimental covariance matrix
// for the absolute differential cross-sections
void wc_covmat(vector<double> & covmatv){

  int const ndat=5;

  //  std::cout<<"\n Computing covariance matrix \n"<<std::endl;

  // Results for the  pt> 25 GeV lepton cut
  // Total cross section (Section 6)
  double xs = 107.7; // [pb]
  double unc_xs_stat = 3.3; // statistical error [pb]
  double unc_xs_syst = 6.9; // systematic error [pb]
  // total experimental uncertainty in pb
  double unc_xs = sqrt(unc_xs_stat*unc_xs_stat + unc_xs_syst*unc_xs_syst);
  //std::cout<<"xsec = "<<xs<<" +- "<<unc_xs<<" pb "<<std::endl;

  //Normalized differential cross section (Table 4)
  double diff_n[ndat] = {0.638, 0.556, 0.527, 0.416, 0.326};

  // Uncertainties in the normalized differential cross section (Table 4)
  double unc_diff_n_stat[ndat] = {0.016, 0.016, 0.015, 0.012, 0.012};
  double unc_diff_n_syst[ndat] = {0.012, 0.012, 0.011, 0.009, 0.009};
  double unc_diff_n[ndat];
  // Total uncertainty in the normalized differential cross section
  // Add everything in quadrature
  for (int i=0; i<5; i++) {
    unc_diff_n[i] = sqrt(unc_diff_n_stat[i]*unc_diff_n_stat[i]
			 + unc_diff_n_syst[i]*unc_diff_n_syst[i]);
  }

  // Covariance and correlation matrices for the normalized differential cross section
  double Vcov_n[ndat][ndat]={{0.0}};
  double Vcor_n[ndat][ndat]={{0.0}};

  // Initializing the correlation matrix of the normalized cross section (Table 5)
  // Question: correlation of the statistical and systematic errors or only
  // of the latter
  for (int i=0; i<ndat; i++) {
    Vcor_n[i][i] = 1.;
  }
  Vcor_n[0][1] = -0.22;
  Vcor_n[0][2] = -0.24;
  Vcor_n[0][3] = -0.26;
  Vcor_n[0][4] = -0.24;
  Vcor_n[1][2] = -0.22;
  Vcor_n[1][3] = -0.26;
  Vcor_n[1][4] = -0.24;
  Vcor_n[2][3] = -0.28;
  Vcor_n[2][4] = -0.26;
  Vcor_n[3][4] = -0.26;
  for (int i=0; i<ndat; i++) {
    for (int j=i+1; j<ndat; j++){
      Vcor_n[j][i] = Vcor_n[i][j];
    }
  }

  // Building the covariance matrix  of the normalized data
  // This implies multiplying the correlation matrix with the total
  // experimental uncertainties in the normalized measurements
  for (int i=0; i<ndat; i++) {
    for (int j=0; j<ndat; j++){
      Vcov_n[i][j] = Vcor_n[i][j]*unc_diff_n[i]*unc_diff_n[j];
    }
  }

  //Absolute differential cross section (Table 6 of the paper)
  double diff[ndat];
  double unc_diff_stat[ndat];
  double unc_diff_syst[ndat];
  double unc_diff[ndat];

  // Covariance and correlation matrices for the absolute differential cross section
  double Vcov[ndat][ndat]={{0.0}};
  double Vcor[ndat][ndat]={{0.0}};

  // Absolute differential cross section
  // Multiply the differential cross-section by the measured absolute cross section
  for (int i=0; i<ndat; i++) {
    diff[i] = diff_n[i]*xs;
  }

  // Uncertainties in the absolute differential cross section
  // constructed from the uncertainties in the normalized measurement and
  // in the absolute total cross section
  for (int i=0; i<ndat; i++) {

    // Statistical uncertainty
    unc_diff_stat[i] = ( unc_diff_n_stat[i]/diff_n[i] )* ( unc_diff_n_stat[i]/diff_n[i] ) + ( unc_xs_stat/xs ) * (unc_xs_stat/xs );

    // multiply the value of the asbolute differentual cross section
    // by the sum in quadrature of the relative stat errors in the absolute total
    // cross section and in the relative differential cross section
    unc_diff_stat[i] = diff[i]*sqrt(unc_diff_stat[i]);

    // Systematic uncertainty
    unc_diff_syst[i] = ( unc_diff_n_syst[i]/diff_n[i] ) * (unc_diff_n_syst[i]/diff_n[i] ) + ( unc_xs_syst/xs ) * ( unc_xs_syst/xs );
    unc_diff_syst[i] = diff[i]*sqrt(unc_diff_syst[i]);
  }

  // multiply the value of the asbolute differentual cross section
  // by the sum in quadrature of the relative stat errors in the absolute total
  // cross section and in the relative differential cross section
  for (int i=0; i<ndat; i++) {
    unc_diff[i] = sqrt(unc_diff_stat[i]*unc_diff_stat[i] + unc_diff_syst[i]*unc_diff_syst[i]);
  }

  // Building the correlation matrix of the absolute differential cross section
  for (int i=0; i<ndat; i++) {
    for (int j=0; j<ndat; j++){
      double s1 = diff_n[i]/sqrt(Vcov_n[i][i])*unc_xs/xs;
      double s2 = diff_n[j]/sqrt(Vcov_n[j][j])*unc_xs/xs;
      Vcor[i][j] = (Vcor_n[i][j]+s1*s2)/sqrt( (1.+s1*s1)*(1.+s2*s2) );
    }
  }
  // Building the covariance matrix of the absolute differential cross section
  for (int i=0; i<ndat; i++) {
    for (int j=0; j<ndat; j++){
      Vcov[i][j] = Vcor[i][j]*unc_diff[i]*unc_diff[j];
    }
  }

  // Save
  for (int i=0; i<ndat; i++) {
    for (int j=0; j<ndat; j++){
      covmatv.push_back(Vcov[i][j]);
    }
  }

}


// Now the cross-section ratio data
// Defined as the ratio of WplusCbar and WminusC production cross-sections
// The covariance matrix for this measurement is diagonal
// And it is not correlated
void CMSWCHARMRATFilter::ReadData()
{
  // Opening files
  fstream f1;
  fstream f2;

  // Read raw data
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
	   << fSetName << "/Wcharm_ratio_Table7.txt";
  f1.open(datafile.str().c_str(), ios::in);
  // Check file properly open
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream theoryfile("");
  theoryfile << dataPath() << "rawdata/"
	     << fSetName << "/theory.txt";
  f2.open(theoryfile.str().c_str(), ios::in);
  // Check file properly open
  if (f2.fail()) {
    cerr << "Error opening data file " << theoryfile.str() << endl;
    exit(-1);
  }
  
  // Now read the data
  string line;
  double etamin, etamax;
  double MW2 = pow(MW,2.0);
  // Center of mass energy
  double s = 7000;

  // Reading data
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);

    lstream >> etamin >> etamax;
    fKin1[i] = (etamax + etamin)*0.5;    // eta
    fKin2[i] = MW2;                      // Mass W squared
    fKin3[i] = s;                        // sqrt(s)

    lstream >> fData[i];

    lstream >> fStat[i];

    double sys=0.0;
    lstream >> sys;
    // The covariance matrix for this observable is diagonal
    // so add in quadrature statistical and systematic uncertainties
    fStat[i] = sqrt ( fStat[i] * fStat[i] + sys * sys );

    getline(f2,line);
    istringstream kstream(line);

    kstream >> fSys[i][0].mult >> fSys[i][1].mult;
    fSys[i][0].mult /= sqrt(2.);
    fSys[i][0].add = fSys[i][0].mult*fData[i]/100;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "SKIP";
    
    fSys[i][1].mult /= sqrt(2.);
    fSys[i][1].add = fSys[i][1].mult*fData[i]/100;
    fSys[i][1].type = MULT;
    fSys[i][1].name = "SKIP";

  }

  f1.close();
  f2.close();

}
