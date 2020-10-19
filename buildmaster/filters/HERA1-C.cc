/**
 *     hera1-av
 *
 *     Combined HERA-I (H1+ZEUS) NC and CC inclusive
 *     structure functions data set
 *
 *     As described in
 *
 *     "Combined Measurement and QCD Analysis of the
 *     Inclusive ep Scattering Cross Sections at HERA"
 *     arXiv:0911.0884v1 [hep-ex]
 *
 *     The combined data set is divided into four sets
 *
 *     1.- d09-158.nce+p.txt -> NC e+ cross sections, 528 data points
 *     File structure:
 *     Bin    Q2          x       y    s_r    F2    sta   unc   cor   exp   rel    gp   had    tot    (d_i, i=1,100)
 *
 *     2.- d09-158.nce-p.txt -> NC e- cross sections, 145 data points
 *     Identical file structure as above
 *
 *     3.-  d09-158.cce+p.txt -> CC e+ cross sections, 34 data points
 *     File structure:
 *     Bin      Q2          x       y     d2s/dxdQ2   s_r    sta   unc   cor   exp   rel    gp   had    tot  (d_i, i=1,100)
 *
 *     4.- d09-158.cce-p.txt -> CC e- cross sections, 34 data points
 *     Identical file structure as above
 *
 *     So we have four sets for this experiment
 *
 *     The overall normalization uncertainty, common to this
 *     dataset, is 0.5%
 *
 *     Note that all of the 110 correlated systematic uncertainties are
 *     symmetric
 */

#include "HERA1-C.h"

void HERA1NCEPFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/d09-158.nce+p.txt";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Filtering data
  double tmp;
  string line;
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    lstream >> tmp;
    lstream >> fKin2[i];   //q2
    lstream >> fKin1[i];   //x
    lstream >> fKin3[i];   //y
    
    // Observable - Neutral current reduced cross-section
    lstream >> fData[i];
    lstream >> tmp;
    
    // Statistical errors - Percentage with respect the observable
    lstream >> fStat[i];
    fStat[i] *= fData[i]*1e-2;
    
    // Uncorrelated systematics
    lstream >> fSys[i][0].mult;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    
    // Total correlated uncertainty
    lstream >> tmp;
    
    // Total experimental uncertainty
    lstream >> tmp;
    
    // Correlated procedural uncertainties
    for (int j = 1; j < 4; j++)
    {
      lstream >> fSys[i][j].mult;
      fSys[i][j].type = ADD;
      ostringstream sysname;
      sysname << "HERA1PROC" << j;
      fSys[i][j].name = sysname.str();
    }   
        
    // Total uncertainty
    lstream >> tmp;
    
    // Normalization uncertainty
    fSys[i][4].mult = 0.5; // Absolute -> 0.5% in combined HERAI set
    fSys[i][4].type = MULT;
    fSys[i][4].name = "HERA1NORM";   
    
    // 110 experimental correlated systematic uncertainties
    for (int icor = 5; icor < fNSys; icor++)
    {
      lstream >> fSys[i][icor].mult;
      fSys[i][icor].type = ADD;
      ostringstream sysname;
      sysname << "HERA1CORRSYS" << icor-4;
      fSys[i][icor].name = sysname.str();
    }
    
    // Convert from percentage to absolute
    for (int isys = 0; isys < fNSys; isys++)
      fSys[i][isys].add = fSys[i][isys].mult*fData[i]*1e-2;
  }
  
  f1.close();
}

/**
 *     hera1-av.f
 *
 *     Combined HERA-I (H1+ZEUS) NC and CC inclusive
 *     structure functions data set
 *
 *     As described in
 *
 *     "Combined Measurement and QCD Analysis of the
 *     Inclusive ep Scattering Cross Sections at HERA"
 *     arXiv:0911.0884v1 [hep-ex]
 *
 *     The combined data set is divided into four sets
 *
 *     1.- d09-158.nce+p.txt -> NC e+ cross sections, 528 data points
 *     File structure:
 *     Bin    Q2          x       y    s_r    F2    sta   unc   cor   exp   rel    gp   had    tot    (d_i, i=1,100)
 *
 *     2.- d09-158.nce-p.txt -> NC e- cross sections, 145 data points
 *     Identical file structure as above
 *
 *     3.-  d09-158.cce+p.txt -> CC e+ cross sections, 34 data points
 *     File structure:
 *     Bin      Q2          x       y     d2s/dxdQ2   s_r    sta   unc   cor   exp   rel    gp   had    tot  (d_i, i=1,100)
 *
 *     4.- d09-158.cce-p.txt -> CC e- cross sections, 34 data points
 *     Identical file structure as above
 *
 *     So we have four sets for this experiment
 *
 *     The overall normalization uncertainty, common to this
 *     dataset, is 0.5%
 *
 *     Note that all of the 110 correlated systematic uncertainties are
 *     symmetric
 */
void HERA1NCEMFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/d09-158.nce-p.txt";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Filtering data
  double tmp;
  string line;
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    lstream >> tmp;
    lstream >> fKin2[i];   //q2
    lstream >> fKin1[i];   //x
    lstream >> fKin3[i];   //y
    
    // Observable - Neutral current reduced cross-section
    lstream >> fData[i];
    lstream >> tmp;
    
    // Statistical errors - Percentage with respect the observable
    lstream >> fStat[i];
    fStat[i] *= fData[i]*1e-2;
    
    // Uncorrelated systematics
    lstream >> fSys[i][0].mult;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    
    // Total correlated uncertainty
    lstream >> tmp;
    
    // Total experimental uncertainty
    lstream >> tmp;
    
    // Correlated procedural uncertainties
    for (int j = 1; j < 4; j++)
    {
      lstream >> fSys[i][j].mult;
      fSys[i][j].type = ADD;
      ostringstream sysname;
      sysname << "HERA1PROC" << j;
      fSys[i][j].name = sysname.str();
    }   
        
    // Total uncertainty
    lstream >> tmp;
    
    // Normalization uncertainty
    fSys[i][4].mult = 0.5; // Absolute -> 0.5% in combined HERAI set
    fSys[i][4].type = MULT;
    fSys[i][4].name = "HERA1NORM";   
    
    // 110 experimental correlated systematic uncertainties
    for (int icor = 5; icor < fNSys; icor++)
    {
      lstream >> fSys[i][icor].mult;
      fSys[i][icor].type = ADD;
      ostringstream sysname;
      sysname << "HERA1CORRSYS" << icor-4;
      fSys[i][icor].name = sysname.str();
    }
    
    // Convert from percentage to absolute
    for (int isys = 0; isys < fNSys; isys++)
      fSys[i][isys].add = fSys[i][isys].mult*fData[i]*1e-2;
  }
  
  f1.close();
}

/**
 *     hera1-av.f
 *
 *     Combined HERA-I (H1+ZEUS) NC and CC inclusive
 *     structure functions data set
 *
 *     As described in
 *
 *     "Combined Measurement and QCD Analysis of the
 *     Inclusive ep Scattering Cross Sections at HERA"
 *     arXiv:0911.0884v1 [hep-ex]
 *
 *     The combined data set is divided into four sets
 *
 *     1.- d09-158.nce+p.txt -> NC e+ cross sections, 528 data points
 *     File structure:
 *     Bin    Q2          x       y    s_r    F2    sta   unc   cor   exp   rel    gp   had    tot    (d_i, i=1,100)
 *
 *     2.- d09-158.nce-p.txt -> NC e- cross sections, 145 data points
 *     Identical file structure as above
 *
 *     3.-  d09-158.cce+p.txt -> CC e+ cross sections, 34 data points
 *     File structure:
 *     Bin      Q2          x       y     d2s/dxdQ2   s_r    sta   unc   cor   exp   rel    gp   had    tot  (d_i, i=1,100)
 *
 *     4.- d09-158.cce-p.txt -> CC e- cross sections, 34 data points
 *     Identical file structure as above
 *
 *     So we have four sets for this experiment
 *
 *     The overall normalization uncertainty, common to this
 *     dataset, is 0.5%
 *
 *     Note that all of the 110 correlated systematic uncertainties are
 *     symmetric
 */
void HERA1CCEPFilter::ReadData()
{
  // Opening files
  fstream f1, g1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/d09-158.cce+p.txt";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Filtering data
  string line;
  double tmp;
  
  for (int i = 0; i < fNData; i++)
    string line;
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    lstream >> tmp;
    lstream >> fKin2[i];   //q2
    lstream >> fKin1[i];   //x
    lstream >> fKin3[i];   //y
    
    // Observable - Charged current reduced cross-section
    lstream >> tmp;
    lstream >> fData[i];
    
    // Statistical errors - Percentage with respect the observable
    lstream >> fStat[i];
    fStat[i] *= fData[i]*1e-2;
    
    // Uncorrelated systematics
    lstream >> fSys[i][0].mult;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    
    // Total correlated uncertainty
    lstream >> tmp;
    
    // Total experimental uncertainty
    lstream >> tmp;
    
    // Correlated procedural uncertainties
    for (int j = 1; j < 4; j++)
    {
      lstream >> fSys[i][j].mult;
      fSys[i][j].type = ADD;
      ostringstream sysname;
      sysname << "HERA1PROC" << j;
      fSys[i][j].name = sysname.str();
    }   
        
    // Total uncertainty
    lstream >> tmp;
    
    // Normalization uncertainty
    fSys[i][4].mult = 0.5; // Absolute -> 0.5% in combined HERAI set
    fSys[i][4].type = MULT;
    fSys[i][4].name = "HERA1NORM";   
    
    // 110 experimental correlated systematic uncertainties
    for (int icor = 5; icor < fNSys; icor++)
    {
      lstream >> fSys[i][icor].mult;
      fSys[i][icor].type = ADD;
      ostringstream sysname;
      sysname << "HERA1CORRSYS" << icor-4;
      fSys[i][icor].name = sysname.str();
    }
    
    // Convert from percentage to absolute
    for (int isys = 0; isys < fNSys; isys++)
      fSys[i][isys].add = fSys[i][isys].mult*fData[i]*1e-2;
  }
  
  f1.close();
}

/**
 *     hera1-av.f
 *
 *     Combined HERA-I (H1+ZEUS) NC and CC inclusive
 *     structure functions data set
 *
 *     As described in
 *
 *     "Combined Measurement and QCD Analysis of the
 *     Inclusive ep Scattering Cross Sections at HERA"
 *     arXiv:0911.0884v1 [hep-ex]
 *
 *     The combined data set is divided into four sets
 *
 *     1.- d09-158.nce+p.txt -> NC e+ cross sections, 528 data points
 *     File structure:
 *     Bin    Q2          x       y    s_r    F2    sta   unc   cor   exp   rel    gp   had    tot    (d_i, i=1,100)
 *
 *     2.- d09-158.nce-p.txt -> NC e- cross sections, 145 data points
 *     Identical file structure as above
 *
 *     3.-  d09-158.cce+p.txt -> CC e+ cross sections, 34 data points
 *     File structure:
 *     Bin      Q2          x       y     d2s/dxdQ2   s_r    sta   unc   cor   exp   rel    gp   had    tot  (d_i, i=1,100)
 *
 *     4.- d09-158.cce-p.txt -> CC e- cross sections, 34 data points
 *     Identical file structure as above
 *
 *     So we have four sets for this experiment
 *
 *     The overall normalization uncertainty, common to this
 *     dataset, is 0.5%
 *
 *     Note that all of the 110 correlated systematic uncertainties are
 *     symmetric
 */
void HERA1CCEMFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/d09-158.cce-p.txt";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Filtering data
  string line;
  double tmp;
  
  for (int i = 0; i < fNData; i++)
    string line;
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    lstream >> tmp;
    lstream >> fKin2[i];   //q2
    lstream >> fKin1[i];   //x
    lstream >> fKin3[i];   //y
    
    // Observable - Charged current reduced cross-section
    lstream >> tmp;
    lstream >> fData[i];
    
    // Statistical errors - Percentage with respect the observable
    lstream >> fStat[i];
    fStat[i] *= fData[i]*1e-2;
    
    // Uncorrelated systematics
    lstream >> fSys[i][0].mult;
    fSys[i][0].type = ADD;
    fSys[i][0].name = "UNCORR";
    
    // Total correlated uncertainty
    lstream >> tmp;
    
    // Total experimental uncertainty
    lstream >> tmp;
    
    // Correlated procedural uncertainties
    for (int j = 1; j < 4; j++)
    {
      lstream >> fSys[i][j].mult;
      fSys[i][j].type = ADD;
      ostringstream sysname;
      sysname << "HERA1PROC" << j;
      fSys[i][j].name = sysname.str();
    }   
        
    // Total uncertainty
    lstream >> tmp;
    
    // Normalization uncertainty
    fSys[i][4].mult = 0.5; // Absolute -> 0.5% in combined HERAI set
    fSys[i][4].type = MULT;
    fSys[i][4].name = "HERA1NORM";   
    
    // 110 experimental correlated systematic uncertainties
    for (int icor = 5; icor < fNSys; icor++)
    {
      lstream >> fSys[i][icor].mult;
      fSys[i][icor].type = ADD;
      ostringstream sysname;
      sysname << "HERA1CORRSYS" << icor-4;
      fSys[i][icor].name = sysname.str();
    }
    
    // Convert from percentage to absolute
    for (int isys = 0; isys < fNSys; isys++)
      fSys[i][isys].add = fSys[i][isys].mult*fData[i]*1e-2;
  }
  
  f1.close();
}
