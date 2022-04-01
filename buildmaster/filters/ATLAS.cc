/**
 * ATLAS JETS 2010 - 1112.6297v1
 *
 * Inclusive jet and dijet data from the ATLAS 2010 dataset
 * Statistics from 36 pb^-1 dataset
 *
 * There are 18 sources of correlated systematics, that
 * for which each rapidity bin is fully correlated in pt,
 * but only some bins in rapidity might be correlated among them
 * So effectively there are 86 source of systematic uncertainties
 * fully correlated between all pt and eta
 *
 * This code is used for both the R=0.4 and R=0.6 analyses, with differing setnames
 *
 * See Note added on 27 Jun 2014 at http://hepdata.cedar.ac.uk/view/ins1082936
 * updated luminosity
 */

#include "ATLAS.h"
#include "NNPDF/utils.h"
#include <numeric>
#include <algorithm>

void ATLAS2010JETSFilter::ReadData()
{
  // Constants
  const int nrapbin  = 7;      // There are seven rapidity bins
  const int nsystype = 18;     // There are 18 systematics per rapidity bin
  const double lcorr = 1.0187; // Correction factor due to luminosity upgrade
  const std::string rawdata_path = dataPath() + "rawdata/" + fSetName + "/";
  
  ifstream hepdata_file(rawdata_path + "/ATLAS-jets-36pb-hepdata.data");
  if (!hepdata_file.is_open()) throw runtime_error("Cannot open file: ATLAS-jets-36pb-hepdata.data");

  ifstream kinematics_file(rawdata_path + "ATLAS-jets-36pb.data");
  if (!kinematics_file.is_open()) throw runtime_error("Cannot open file: ATLAS-jets-36pb.data");

  ifstream correlation_file(rawdata_path + "ATLAS-jets-36pb.sysmat");
  if (!correlation_file.is_open()) throw runtime_error("Cannot open file: ATLAS-jets-36pb.sysmat");

  ifstream systematics_file(rawdata_path + "ATLAS-jets-36pb.sys");
  if (!systematics_file.is_open()) throw runtime_error("Cannot open file: ATLAS-jets-36pb.sys");
    
  ifstream nptcorr_file(rawdata_path + "ATLAS-jets-36pb.npt");
  if (!nptcorr_file.is_open()) throw runtime_error("Cannot open file: ATLAS-jets-36pb.npt");

  // Reading data, non-perturbative corrections, and kinematics 
  for (int i = 0; i < fNData; i++)
  {
    double dummy;
    
    // Read data from hepdata file
    hepdata_file >> fKin2[i];           // pT
    hepdata_file >> dummy >> dummy;     // pTmax and pTmin (not needed)
    hepdata_file >> fData[i];           // Data point
    hepdata_file >> fStat[i] >> dummy;  // Statistical uncertainty (symmetric)

    // Read nonperturbative correction from file
    double nptcorr, npterr;
    nptcorr_file >> nptcorr >> npterr;

    // Convert pT to pT^2
    fKin2[i] *= fKin2[i];
    fKin3[i] = 7000;                // sqrt(s)
    
    // Apply luminosity and nonperturbative correction
    fData[i] *= lcorr/nptcorr;
    fStat[i] *= lcorr/nptcorr;

    // Uncertainty due to luminosity
    fSys[i][0].mult = 3.5;           //Luminosity uncertainty of 3.5% (updated in 27 JUN 2014)
    fSys[i][0].name = "ATLASLUMI10";

    // Uncertainty due to nonperturbative correction
    fSys[i][1].mult = npterr/nptcorr*100.;    // Error on nonperturbative correction (%)
    fSys[i][1].name = "CORR";

    hepdata_file >> dummy >> dummy;   //total (correlated) systematic uncertainty (asymmetric)
    hepdata_file >> dummy >> dummy;  //uncorrelated systematic uncertainty (symmetric)
  }

  // Read data file - need rapidity bins to read systematics file
  int ndatbin[nrapbin];

  // Test number of rapidity bins
  int test_nrapbin; kinematics_file >> test_nrapbin;
  if (test_nrapbin != nrapbin) throw runtime_error("Mismatch between number of rapidity bins");
    
  // Read pseudorapidity bin information (bin value, number of points in bin)
  int idat = 0;
  for (int i = 0; i < nrapbin; i++)
  {
    // Pesudorapidity values
    double etamin, etamax;
    kinematics_file >> etamin >> etamax;
    const double eta = (etamax+etamin)*0.5;

    // Number of datapoints in pseudorapidity bin
    kinematics_file >> ndatbin[i];

    // Catch the \n
    string dummy_line; getline(kinematics_file,dummy_line);
    // Set eta for points in bin
    for (int j = 0; j < ndatbin[i]; j++)
    {
      fKin1[idat] = eta; // Pseudorapidity value
      //Rest of values in file are already in the hepdata file
      getline(kinematics_file,dummy_line);
      idat++;
    }
  }

  // Test pseudorapidity bin information
  const int total_datapoints = std::accumulate(ndatbin, ndatbin+nrapbin, 0); 
  if (total_datapoints != GetNData()) throw runtime_error("Mismatch between number of data points: " + to_string(total_datapoints) + " vs " + to_string(GetNData()));

  // Reading information on how rapidity bin-by rapidity errors should be correlated 
  int sysnumber[nrapbin][nsystype];
  for (int isys = 0; isys < nsystype; isys++)
  {
    string line; getline(correlation_file,line);
    // Clear comment at beginning of each line
    line.erase(0, line.find_last_of("\"")+1);
    const vector<int> split_line = NNPDF::isplit(line);

    if (split_line.size() != nrapbin) 
        throw runtime_error("Mismatch between number of rapidity bins and systematic error designations: " + to_string(split_line.size())+" vs "+to_string(nrapbin));

    // Read in number identifing systematic
    // Save [0] for luminosity, [1] nonperturbative uncertainty, [2]-[4] for uncorrelated systematics
    for (int j = 0; j < nrapbin; j++)
      sysnumber[j][isys] = split_line[j] + 4;      
  }

  //Reading systematics
  for (int i = 0; i < nrapbin; i++)
  {
    // Skip first line for every rapidity bin
    string dummy_string; getline(systematics_file,dummy_string);
    for (int j = 0; j < ndatbin[i]; j++)
    {
      double dummy; systematics_file >> dummy >> dummy;          //ptmax and ptmin
      const int idat = std::accumulate(ndatbin, ndatbin + i, j); // Current datapoint

      // Firstly, let's set the inactive systypes for this rapidity bin to zero
      // Correlated systematics block starts at systematic 5
      for (int l = 5; l < fNSys; l++)
        if (std::find(std::begin(sysnumber[i]), std::end(sysnumber[i]), l)  == std::end(sysnumber[i]))
            fSys[idat][l].mult = 0;

      // Now read the relevant systematics
      for (int isystype = 0; isystype < nsystype; isystype++)
      {
        double up, down, symmetrised;
        systematics_file >> up >> down;
        symmetriseErrors(up,down,&symmetrised,&dummy);
        fSys[idat][sysnumber[i][isystype]].mult = symmetrised;   //systematics are in %
        fSys[idat][sysnumber[i][isystype]].name = "CORR";
      }

      // Last few errors are uncorrelated
      for (int j = 2; j < 5; j++)
      {
        systematics_file >> fSys[idat][j].mult >> dummy;   //uncorrelated systematics are symmetric (%)
        fSys[idat][j].name = "UNCORR";
      }
    }
      systematics_file >> dummy_string; // Eat the newline
  }

  for (int i=0; i<fNData; i++)
    for (int l = 0; l < fNSys; l++)
    {
      fSys[i][l].add = fSys[i][l].mult*fData[i]*1e-2;    // Calculate absolute systematics for additive part
      fSys[i][l].type = MULT;                            // All systematics multiplicative
    }

  hepdata_file.close();
  kinematics_file.close();
  correlation_file.close();
  systematics_file.close();
  nptcorr_file.close();
}


/**
 * ATLASR04JETS2P76TEV filter
 * 1304.4739v1
 *
 */
void ATLASR04JETS2P76TEVFilter::ReadData()
{
  // Constants
  const int nrapbin = 7;   // There are seven rapidity bins
  const int nsystype = 19; // There are 19 systematics per rapidity bin (ignoring luminosity for the moment)
  const std::string rawdata_path = dataPath() + "rawdata/" + fSetName + "/";
  
  ifstream rawdata_file(rawdata_path + "/ATLASR04JETS2P76TEV_raw.data");
  if (!rawdata_file.is_open()) throw runtime_error("/ATLASR04JETS2P76TEV_raw.data");

  ifstream correlation_file(rawdata_path + "/ATLASR04JETS2P76TEV_raw_SYSTYPE.data");
  if (!correlation_file.is_open()) throw runtime_error("/ATLASR04JETS2P76TEV_raw_SYSTYPE.data");

  // Reading information on how rapidity bin-by rapidity errors should be correlated 
  int sysnumber[nrapbin][nsystype];
  for (int isys = 0; isys < nsystype; isys++)
  {
    string line; getline(correlation_file,line);
    // Clear comment at beginning of each line
    line.erase(0, line.find_last_of("\"")+1);
    const vector<int> split_line = NNPDF::isplit(line);

    if (split_line.size() != nrapbin) 
        throw runtime_error("Mismatch between number of rapidity bins and systematic error designations: " + to_string(split_line.size())+" vs "+to_string(nrapbin));

    // Read in number identifing systematic
    //Save [0] for luminosity, [1] nonperturbative uncertainty, [2]-[3] for uncorrelated systematics
    for (int j = 0; j < nrapbin; j++)
      sysnumber[j][isys] = split_line[j] + 3;      
  }

  //Reading data
  string line;
  int ndatbin[nrapbin];
  int idat = 0;
  for (int i = 0; i < nrapbin; i++)
  {
    double etamin, etamax, dummy;
    rawdata_file >> etamin >> etamax;
    const double eta = (etamax+etamin)*0.5;

    rawdata_file >> ndatbin[i];
    string dummy_string; getline(rawdata_file, dummy_string); // Eat up newline
    for (int j = 0; j < ndatbin[i]; j++)
    {
      // Read pT bin limits
      double ptmin, ptmax;
      rawdata_file >> ptmin >> ptmax;

      // Set kinematics
      fKin1[idat] = eta;                      //eta
      fKin2[idat] = pow((ptmin+ptmax)*0.5,2); //pt2
      fKin3[idat] = 2760;                     //sqrt(s)

      fSys[idat][0].mult = 2.7;      //2.7% luminosity uncertainty (uncorrelated with 36PB luminosity uncertainty)
      fSys[idat][0].name = "CORR";

      // Nonperturbative correction and associated uncertainty
      double up, down, symmetrised, nonpert;
      rawdata_file >> nonpert >> up >> down;
      symmetriseErrors(up,down,&symmetrised,&dummy);
      fSys[idat][1].mult = symmetrised*100/nonpert;
      fSys[idat][1].name = "CORR";

      rawdata_file >> fData[idat];
      fData[idat] /= nonpert;           // apply nonperturbative correction
                                        // note: multiplicative correction so additive systematics should be calcuated with corrected data

      rawdata_file >> fStat[idat];      // statistical uncertainty in %
      fStat[idat] *= fData[idat]*0.01;  // convert to absolute
      
      // Firstly, let's set the inactive systypes for this rapidity bin to zero
      // Correlated systematics block starts at systematic 4
      for (int l = 4; l < fNSys; l++)
        if (std::find(std::begin(sysnumber[i]), std::end(sysnumber[i]), l)  == std::end(sysnumber[i]))
            fSys[idat][l].mult = 0;

      //First 5 systematics
      for (int isystype = 0; isystype < 5; isystype++)
      {
        double up, down, symmetrised;
        rawdata_file >> up >> down;
        symmetriseErrors(up,down,&symmetrised,&dummy);
        fSys[idat][sysnumber[i][isystype]].mult = symmetrised;   //systematics are in %
        fSys[idat][sysnumber[i][isystype]].name = "CORR";
      }

      //Systematics number 88 and 89 are correlated in quadrature with 31 from ATLASJETS36PB
      // NH NOTE: They don't seem to be being correlated!
      rawdata_file >> up >> down;
      symmetriseErrors(up, down, &symmetrised, &dummy);
      const double sys88 = symmetrised; //systematics are in %

      //Next 7 systematics
      for (int isystype = 6; isystype < 13; isystype++)
      {
        double up, down, symmetrised;
        rawdata_file >> up >> down;
        symmetriseErrors(up, down, &symmetrised, &dummy);
        fSys[idat][sysnumber[i][isystype]].mult = symmetrised;   //systematics are in %
        fSys[idat][sysnumber[i][isystype]].name = "CORR";
      }

      //Systematics number 88 and 89 are correlated in quadrature with 31 from ATLASJETS36PB
      // NH NOTE: They don't seem to be being correlated!
      rawdata_file >> up >> down;
      symmetriseErrors(up,down, &symmetrised, &dummy);
      const double sys89 = symmetrised; //systematics are in %

      fSys[idat][sysnumber[6][4]+1].mult=sqrt(sys88*sys88+sys89*sys89);  //in quadrature for sys 31 of this set
      fSys[idat][sysnumber[6][4]+1].name = "CORR";

      //Next 5 systematics - symmetric
      for (int isystype = 14; isystype < nsystype; isystype++)
      {
        rawdata_file >> fSys[idat][sysnumber[i][isystype]].mult;
        fSys[idat][sysnumber[i][isystype]].name = "CORR";
      }

      //Two sources of uncorrelated systematics (%)
      rawdata_file >> fSys[idat][2].mult;
      fSys[idat][2].name = "UNCORR";

      rawdata_file >> fSys[idat][3].mult;
      fSys[idat][3].name = "UNCORR";

      idat++;
    }
  }

  // Final processing
  for (int i=0; i<fNData; i++)
    for (int l = 0; l < fNSys; l++)
    {
      fSys[i][l].add = fSys[i][l].mult*fData[i]*1e-2;    // Calculate absolute systematics for additive part
      fSys[i][l].type = MULT;                            // All systematics multiplicative
    }

  rawdata_file.close();
  correlation_file.close();
}

/**
 * atlas-wz
 *
 * This data is taken from the paper arxiv:1109.5141
 *
 * It contains the combined lepton rapidity distributions
 * for W+, W- and Z production, with the full covariance matrix
 * available
 *
 * Fit separatelt W+ and W- with covariance matrix instead of
 * the asymmetry distribution
 *
 * The that the results are rather close to a direct
 * refitting with the asymmetry
 * distribution
 *
 * This data will be included using reweighting, with DYNNLO
 * for the computation of the theory predictions
 *
 * Put all data together into a single set
 * This is due to the strong correlation between the W+, W^- and
 * Z0 which makes not suitable to treat them separately
 * in the first stages of the minimization
 */
void ATLASWZRAP36PBFilter::ReadData()
{
  // Opening files
  fstream fWZ[3];

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/ATLAS-36pb-Wplrap.data";
  fWZ[0].open(datafile.str().c_str(), ios::in);

  if (fWZ[0].fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream datafile2("");
  datafile2 << dataPath() << "rawdata/"
  << fSetName << "/ATLAS-36pb-Wmlrap.data";
  fWZ[1].open(datafile2.str().c_str(), ios::in);

  if (fWZ[1].fail()) {
    cerr << "Error opening data file " << datafile2.str() << endl;
    exit(-1);
  }

  stringstream datafile3("");
  datafile3 << dataPath() << "rawdata/"
  << fSetName << "/ATLAS-36pb-Zrap.data";
  fWZ[2].open(datafile3.str().c_str(), ios::in);

  if (fWZ[2].fail()) {
    cerr << "Error opening data file " << datafile3.str() << endl;
    exit(-1);
  }

  // Starting filter
  const double lcorr = 1.0187; // correction factor due to luminosity upgrade
  const int ndataWZ[3] = {11,11,8};  //Number of data for W+, W- and Z respectively
  const double convfac = lcorr*1000.; // Must multiply from pb to fb
  const double MWZ2[3] = {pow(MW,2.0), pow(MW,2.0), pow(MZ,2.0)};   //Mass squared of W (+ and -) and Z

  string line;
  int idat = 0;
  double etamin,etamax,tmp;

  cout << "********** WARNING: Converting pb to fb to match ApplGrid output ********" << endl;

  for (int iWZ = 0; iWZ < 3; iWZ++)
  {
    // rapidity
    getline(fWZ[iWZ],line);
    istringstream lstream(line);
    for (int i = 0; i < ndataWZ[iWZ]; i++)
    {
      lstream >> etamin >> etamax;
      fKin1[idat+i] = etamin + (etamax-etamin)*0.5;
    }

    // M_W,Z
    for (int i = 0; i < ndataWZ[iWZ]; i++)
      fKin2[idat+i] = MWZ2[iWZ];

    // sqrt(s)
    for (int i = 0; i < ndataWZ[iWZ]; i++)
      fKin3[idat+i] = 7000;

    // obs
    getline(fWZ[iWZ],line);
    istringstream lstream2(line);
    for (int i = 0; i < ndataWZ[iWZ]; i++)
      {
        lstream2 >> fData[idat+i];
        fData[idat+i] *= convfac;
      }

    // stat (%, converted later)
    getline(fWZ[iWZ],line);
    istringstream lstream3(line);
    for (int i = 0; i < ndataWZ[iWZ]; i++)
      lstream3 >> fStat[idat+i];

    // uncorrelated sys
    getline(fWZ[iWZ],line);
    istringstream lstream4(line);
    for (int i = 0; i < ndataWZ[iWZ]; i++)
    {
      lstream4 >> fSys[idat+i][0].mult;
      fSys[idat+i][0].name = "UNCORR";
    }

    // total correlated sys (unused)
    getline(fWZ[iWZ],line);

    // total uncertainty (unused)
    getline(fWZ[iWZ],line);

    // correlated systematics
    for (int isys = 2; isys < fNSys; isys++)  //2 to skip uncorr and lumi
    {
      getline(fWZ[iWZ],line);
      istringstream lstream(line);
      lstream >> tmp;
      for (int i = 0; i < ndataWZ[iWZ]; i++)
      {
        lstream >> fSys[idat+i][isys].mult;
        fSys[idat+i][isys].name = "CORR";
      }
    }

    // luminosity: 3.4%
    for (int i = 0; i < ndataWZ[iWZ]; i++)
    {
      fSys[idat+i][1].mult = 3.5;
      fSys[idat+i][1].name = "ATLASLUMI10";
    }

    idat+=ndataWZ[iWZ];
  }

  // Convert additive uncertainties to absolute form
  for (int i = 0; i < fNData; i++)
  {
    fStat[i] *= fData[i]*1e-2;
    for(int l = 0; l < fNSys; l++)
    {
      fSys[i][l].type = MULT; // All systematics multiplicative
      fSys[i][l].add = fSys[i][l].mult*fData[i]*1e-2;
    }
  }


  fWZ[0].close();
  fWZ[1].close();
  fWZ[2].close();
}


/**
 * atlas-highmass-dy
 *
 * data taken from: http://hepdata.cedar.ac.uk/view/ins1234228
 * and from paper: http://arxiv.org/abs/ARXIV:1305.4192
 *
 */
void ATLASZHIGHMASS49FBFilter::ReadData()
{
  cout << "WARNING: kinematics are not implemented" << endl;
  // Opening files
  fstream f1, f2;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/ATLAS-49fb-Zhighmass.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream datafile2("");
  datafile2 << dataPath() << "rawdata/"
  << fSetName << "/ATLAS-49fb-Zhighmass.sys";
  f2.open(datafile2.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double mbin[fNData+1], stat;

  // Reading data
  string tmp;
  getline(f1,tmp);
  for (int i = 0; i < 6; i++) getline(f2,tmp);

  // Filtering data
  for (int i = 0; i < fNData; i++)
  {
    f1 >> mbin[i] >> mbin[i+1] >> fData[i] >> tmp >> tmp >> tmp;
    fData[i] *= 1E3; // converting pb to fb

    f2 >> tmp >> tmp >> stat;
    for (int j = 0; j < fNSys; j++) f2 >> fSys[i][j].mult;

    fStat[i] = stat*fData[i]*1e-2;
    fKin2[i] = pow( 0.5*(mbin[i] + mbin[i+1]) , 2.0);
    fKin1[i] = 0.5*(mbin[i] + mbin[i+1]);

    fKin3[i] = 7E3;

    for (int j = 0; j < fNSys-1; j++)
    {
      fSys[i][j].add = fSys[i][j].mult*fData[i]*1e-2;
      fSys[i][j].type = MULT;
      if (j < 2)
        fSys[i][j].name = "UNCORR"; // for the uncorrelated
      else
        fSys[i][j].name = "CORR"; // for the correlated
    }

    fSys[i][fNSys-1].add = fSys[i][fNSys-1].mult*fData[i]*1e-2;
    fSys[i][fNSys-1].type = MULT;
    fSys[i][fNSys-1].name = "ATLASLUMI11";

  }

  f1.close();
  f2.close();
}

/**
 * ATLASWPT31PB data
 *
 * data taken from: http://hepdata.cedar.ac.uk/view/ins941555
 * and from paper: http://arxiv.org/abs/ARXIV:1108.6308v2
 *
 */
void ATLASWPT31PBFilter::ReadData()
{
  // Opening files
  fstream f1, f2;
  int idum,jdum;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/ATLASWPT31PB.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream datafile2("");
  datafile2 << dataPath() << "rawdata/"
  << fSetName << "/ATLASWPT31PB.covmat";
  f2.open(datafile2.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double mbin[fNData+1], systot;

  // Reading data
  string tmp;
  getline(f1,tmp);

  // Filtering data
  for (int i = 0; i < fNData; i++)
  {

    // Data are normalized to the total cross section, Units are GeV-1
    f1 >> fKin1[i] >> mbin[i] >> mbin[i+1] >> fData[i] >> fStat[i] >> systot;
    fKin2[i] = pow(MZ,2.0);
    fKin3[i] = 7E3;

    // Question: statistical uncertainties are also given in the covariance
    // matrix. Shall I ignore them and just use the total stat uncertainty given
    // in tha above file for each bin or shall I ignore the above and use the
    // one given in the cov matrix file?

  }

  // Reading covmat
  string line;
  double statmat[fNData][fNData];
  double sysmat[fNData][fNData];
  double** covmat = new double*[fNData];
  for(int i = 0; i < fNData; i++)
  {
    covmat[i] = new double[fNData];
    for(int j = i; j < fNData; j++)
    {
      getline(f2,line);
      istringstream lstream(line);
      lstream >> idum >> jdum >> idum >> jdum >> statmat[i][j] >> sysmat[i][j];
      covmat[i][j] = sysmat[i][j] + statmat[i][j]; // Here I consistently do not add the Stat uncertainty
    }                                // Unless I do not consider the one above. Not sure!
  }

  // Make it symmetric
  for(int i = 0; i < fNData; i++)
  for(int j = 0; j < i; j++)
    covmat[i][j] = covmat[j][i];

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
    for (int l = 0; l < fNSys; l++)
    {
      fSys[i][l].add = syscor[i][l];
      fSys[i][l].mult = fSys[i][l].add*100/fData[i];
      fSys[i][l].type = MULT;
      fSys[i][l].name = "CORR";
    }

  for(int i = 0; i < fNData; i++) 
    {
      delete[] syscor[i];
      delete[] covmat[i];
    }
  delete[] syscor;
  delete[] covmat;
  
  f1.close();
  f2.close();
}
