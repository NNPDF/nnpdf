
#include "NMC.h"

/**
 * NMCpd: hep-ex/9611022, M. Arneodo et al., Nucl. Phys. B 483 (1997) 3
 *
 *     EXPeriment     = CERN-NA-037
 *     REACtion       = muon    p --> muon X
 *     REACtion       = muon deut --> muon X
 *     Plab           = 90+120+200+280 GeV
 *     Collaboration  = NMC
 *     Author         = Arneodo et al
 *     REFerence      = Nucl. Phys. B487 (1997) 3
 *     Additional info: The structure function ratio F2d/F2p measured by
 *     the NMC   Collaboration in muon p(deut)
 *     deep inelastic scattering from the merged data sets at
 *     incident momenta 90, 120, 200 and 280 GeV. Data cover
 *     the x range 0.001 to 0.8 and Q^2 range from 0.1 to
 *     145 GeV^2.
 *     This data is  the final full NMC data set and
 *     supersedes any previous NMC data.
 *
 *=====================================================================
 *
 *     The systematic error is the quadratic sum of the contributions given
 *     in the last 5 columns. These are:
 *
 *     VX - the error on the correction due to vertex resolution
 *     SM - the error due to kinematic resolution
 *     RC - the quadratic sum of errors from radiative corrections and the
 *     functional form of the ratio parametrisations
 *     E  - the error due to the uncertainty of the incident muon momentum
 *     E' - the error due to the uncertainty of the scattered muon momentum.
 *
 *     The signs of E and E' correspnd to an increase in the respective muon
 *     energy.
 *
 *     x       Q^2  F2d/F2p     Errors       VX    SM    RC     E     E'
 *             GeV^2             stat     sys     %     %     %     %     %
 *
 */
void NMCpdFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/nmc_f2df2p.data";
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
    lstream >> fKin1[i]; 
    lstream >> fKin2[i];

    lstream >> fData[i];
    
    lstream >> fStat[i];
    
    double tmp;
    lstream >> tmp;
    
    for (int l = 0; l < fNSys; l++)
    {
      lstream >> fSys[i][l].mult;
      fSys[i][l].add = fSys[i][l].mult*fData[i]*1e-2;
      fSys[i][l].type = ADD;
      fSys[i][l].name = "CORR";
    }
  }
  
  f1.close();
}

/**
 *     NMC: hep-ph/9610231, M. Arneodo et al., Nucl. Phys. B 483 (1997) 3
 *
 *     nmc_targetbeamenergy.data (Table 3 and 4):
 *     x, Q2, y, xsec, E, E', AC, RC, RE, R, F2, stat, systot
 *
 *     The correlated sys error are
 *
 *     E and E' are correlated for different targets and independent for different beam energies
 *     RC are correlated for different beam energies and independent for different targets
 *     AC and RE are fully correlated for all data sets
 *
 *     they are given as percentage, while stat and systot are absolute values.
 *
 *     The normalization uncertainty is 2% and it is correlated for data points
 *     of the same beam energy (even if the target is different, see sect 4.1),
 *     but uncorrelated between different beam energies.
 *
 *     In order to allow for loops over the total number sys and norm and to take
 *     into account correlations properly we have proceeded this way:
 *     for RC, say, we have correlation for data points of the same target, thus
 *     we use two sys for RC one of which is set to zero for the proton,
 *     and the other one is set to zero for the deuteron. Analogously for other sys.
 *
 *     From 02/2011, use reduced cross sections to fit the NMC data
 *     intead of their determination of the structure function
 *       Only central value of the observable modified, uncertainties
 *     (in percentage) are left unchanged
 *
 *     The observable is now DIS_NCE (muon mass and charge effects neglected)
 *     The reduced cross section as used by the HERA experiments
 *
 *     Note that the units of the NCM cross section measurement are
 *     in barns/GeV2
 *     The conversion factor is 1 GeV^-2 = 0.3894 mb
 *     -> 1 b/GeV2 = 2568.05 1/GeV4
 *     [convfac = 0.389379304D9    ! conversion from GeV to picobarn ]
 */
void NMCFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3, f4, d1, d2, d3, d4;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/nmc_p90.data";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  stringstream datafile2("");
  datafile2 << dataPath() << "rawdata/"
  << fSetName << "/nmc_p120.data";
  f2.open(datafile2.str().c_str(), ios::in);
  
  if (f2.fail()) {
    cerr << "Error opening data file " << datafile2.str() << endl;
    exit(-1);
  }
  
  stringstream datafile3("");
  datafile3 << dataPath() << "rawdata/"
  << fSetName << "/nmc_p200.data";
  f3.open(datafile3.str().c_str(), ios::in);
  
  if (f3.fail()) {
    cerr << "Error opening data file " << datafile3.str() << endl;
    exit(-1);
  }
  
  stringstream datafile4("");
  datafile4 << dataPath() << "rawdata/"
  << fSetName << "/nmc_p280.data";
  f4.open(datafile4.str().c_str(), ios::in);
  
  if (f4.fail()) {
    cerr << "Error opening data file " << datafile4.str() << endl;
    exit(-1);
  }
  
  stringstream datafiled("");
  datafiled << dataPath() << "rawdata/"
  << fSetName << "/nmc_d90.data";
  d1.open(datafiled.str().c_str(), ios::in);
  
  if (d1.fail()) {
    cerr << "Error opening data file " << datafiled.str() << endl;
    exit(-1);
  }
  
  stringstream datafiled2("");
  datafiled2 << dataPath() << "rawdata/"
  << fSetName << "/nmc_d120.data";
  d2.open(datafiled2.str().c_str(), ios::in);
  
  if (d2.fail()) {
    cerr << "Error opening data file " << datafiled2.str() << endl;
    exit(-1);
  }
  
  stringstream datafiled3("");
  datafiled3 << dataPath() << "rawdata/"
  << fSetName << "/nmc_d200.data";
  d3.open(datafile3.str().c_str(), ios::in);
  
  if (d3.fail()) {
    cerr << "Error opening data file " << datafiled3.str() << endl;
    exit(-1);
  }
  
  stringstream datafiled4("");
  datafiled4 << dataPath() << "rawdata/"
  << fSetName << "/nmc_d280.data";
  d4.open(datafiled4.str().c_str(), ios::in);
  
  if (d4.fail()) {
    cerr << "Error opening data file " << datafiled4.str() << endl;
    exit(-1);
  }
  
  
  // Starting filter  
  double datain[584][13];
  double relnorbeam = 2.0;
    
  // Reading data
  string line;
  int nini = 0;
  
  for (int i = 0; i < 73; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    for (int j = 0; j < 13; j++)
      lstream >> datain[i][j];
    
    fSys[i][12].mult = relnorbeam;
    fSys[i][0].mult = datain[i][4];
    fSys[i][1].mult = datain[i][5];
    
  }
  nini += 73;  

  for (int i = nini; i < nini+65; i++)
  {
    getline(f2,line);
    istringstream lstream(line);
    for (int j = 0; j < 13; j++)
      lstream >> datain[i][j];
    
    fSys[i][13].mult = relnorbeam;
    fSys[i][2].mult = datain[i][4];
    fSys[i][3].mult = datain[i][5];
  }
  nini += 65;  

  for (int i = nini; i < nini+75; i++)
  {
    getline(f3,line);
    istringstream lstream(line);
    for (int j = 0; j < 13; j++)
      lstream >> datain[i][j];
    
    fSys[i][14].mult = relnorbeam;
    fSys[i][4].mult = datain[i][4];
    fSys[i][5].mult = datain[i][5];
  }
  nini += 75;  

  for (int i = nini; i < nini+79; i++)
  {
    getline(f4,line);
    istringstream lstream(line);
    for (int j = 0; j < 13; j++)
      lstream >> datain[i][j];
    
    fSys[i][15].mult = relnorbeam;
    fSys[i][6].mult = datain[i][4];
    fSys[i][7].mult = datain[i][5];
  }
  nini +=79;

  for (int i = 0; i < 292; i++)
    fSys[i][8].mult = datain[i][7];
  if (fNData > 292 )
  {
  // reading deuteron data
  for (int i = nini; i < nini+73; i++)
  {
    getline(d1,line);
    istringstream lstream(line);
    for (int j = 0; j < 13; j++)
      lstream >> datain[i][j];
    
    fSys[i][12].mult = relnorbeam;
    fSys[i][0].mult = datain[i][4];
    fSys[i][1].mult = datain[i][5];
  }
  nini += 73;  

  for (int i = nini; i < nini+65; i++)
  {
    getline(d2,line);
    istringstream lstream(line);
    for (int j = 0; j < 13; j++)
      lstream >> datain[i][j];
    
    fSys[i][13].mult = relnorbeam;
    fSys[i][2].mult = datain[i][4];
    fSys[i][3].mult = datain[i][5];
  }
  nini += 65;  

  for (int i = nini; i < nini+75; i++)
  {
    getline(d3,line);
    istringstream lstream(line);
    for (int j = 0; j < 13; j++)
      lstream >> datain[i][j];
    
    fSys[i][14].mult = relnorbeam;
    fSys[i][4].mult = datain[i][4];
    fSys[i][5].mult = datain[i][5];
  }
  nini += 75;  

  for (int i = nini; i < nini+79; i++)
  {
    getline(d4,line);
    istringstream lstream(line);
    for (int j = 0; j < 13; j++)
      lstream >> datain[i][j];
    
    fSys[i][15].mult = relnorbeam;
    fSys[i][6].mult = datain[i][4];
    fSys[i][7].mult = datain[i][5];
  }
  nini +=79;

  for (int i = 292; i < fNData; i++)
    fSys[i][9].mult = datain[i][7];
  }
  // Filtering data
  double obs, x, y, rr;
  double q2, rxsec_qed;
  double mp = 0.93799999999999994;
  
  for (int i = 0; i < fNData; i++)
  {
    x  = datain[i][0];
    y  = datain[i][2];
    q2 = datain[i][1];
    
    // Compute reduced cross sections - from measured cross section
    // uncorrected by QCD effects
    double yplus = 1.0 + pow(1.0 - y, 2.0);
    //rxsec_noqed[i] = ( (x[i] * q2[i]*q2[i])
    //                   / (2.0 * M_PI * alphaEM*alphaEM * yplus))
    //                  * datain[i][3] * convfact;
    
    // Set structure function as determined by NMC collaboration
    obs = datain[i][10];
    
    // Compute reduced cross sections - from the published value
    // of f2 and the used value of R in each bin
    // Account for target mass effects, but neglect
    // muon mass effects
    // This should lead (?) to the QED corrected reduced cross section
    // to be used in the PDF fit
    rr = datain[i][9];
    double tmc = mp*mp * x*x * y*y / q2;
    rxsec_qed = (2.0 * obs / yplus) *
    (1.0 - y - tmc + (y*y + 4.0*tmc) / (2.0 * (1.0 + rr) ) );
    
    // Check
    // Reduced cross section always smaller than F2 due to
    // the negative contribution of the longitudinal structure function
    if (rxsec_qed >= obs)
    {
      cerr << "Problem with NMC data" << endl;
      exit(-1);
    }
    
    fSys[i][10].mult = datain[i][6];
    fSys[i][11].mult = datain[i][8];
    
    fKin1[i] = x;
    fKin2[i] = q2;
    fKin3[i] = y;
    fData[i] = rxsec_qed;
    fStat[i] = datain[i][11] / obs * rxsec_qed;
    
    for (int l = 0; l < fNSys; l++)
    {
      fSys[i][l].add = fSys[i][l].mult * fData[i] * 1e-2;
      fSys[i][l].type = ADD;
      fSys[i][l].name = "CORR";
    }
    for (int l = 12; l < fNSys; l++)
      fSys[i][l].type = MULT;
    
      
  }
  
  f1.close();
  f2.close();
  f3.close();
  f4.close();
  
  d1.close();
  d2.close();
  d3.close();
  d4.close();
}

