/*

Experiment: CERN-LHC-CMS (CMS)
Preprinted as CERN-PH-EP-2015-059
Preprinted as CMS-SMP-13-013
Archived as: ARXIV:1504.03511
Auxiliary Material: https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSMP13013
Record in: INSPIRE
Record in: CERN Document Server
Record in: HEPData http://hepdata.cedar.ac.uk/view/ins1359450

Description of the measurement
CERN-LHC. Measurements of the double-differential Drell-Yan cross sections are presented using an integrated luminosity of 19.7 inverse femtobarns in the dimuon channel of proton-proton collision data recorded in 2012 with the CMS detector at the LHC at sqrt(s) = 8 TeV. Covariance matrices are provided in addition to the tables given in the paper.

Description of the buildmaster implementation
Measured double differential fiducial cross section normalised to the inclusive fiducial cross section. The uncertainty indicates the total experimental uncertainties (statistical and systematic added in quadrature).

Raw data and covariance matrix are from HepData:
http://hepdata.cedar.ac.uk/view/ins1359450

*/

#include "CMSZDIFF12.h"

void CMSZDIFF12Filter::ReadData()
{
  // Opening files
  fstream f1, f2;

  bool Table1 = false; // Set to true or false (Table1 or Table 2 respectively)
                       // false -> unnormalised, true -> normalised
  if (Table1)
  {
    stringstream datafile("");
    datafile << dataPath() << "rawdata/"
    << fSetName << "/CMS_Z_DIFF_TABLE1.data";
    f1.open(datafile.str().c_str(), ios::in);

    if (f1.fail()) {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

    stringstream covfile("");
    covfile << dataPath() << "rawdata/"
    	    << fSetName << "/CMS_Z_DIFF_covmat.data";
    f2.open(covfile.str().c_str(), ios::in);

    if (f2.fail()) {
      cerr << "Error opening data file " << covfile.str() << endl;
      exit(-1);
    }
    cout << "We are using Table 1 for CMSZDIFF12.h" << endl;
  }
  else
  {
    stringstream datafile("");
    datafile << dataPath() << "rawdata/"
    << fSetName << "/CMS_Z_DIFF_TABLE2.data";
    f1.open(datafile.str().c_str(), ios::in);

    if (f1.fail()) {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

    stringstream covfile("");
    covfile << dataPath() << "rawdata/"
    	    << fSetName << "/CMS_Z_DIFF_covmat2.data";
    f2.open(covfile.str().c_str(), ios::in);

    if (f2.fail()) {
      cerr << "Error opening data file " << covfile.str() << endl;
      exit(-1);
    }
    cout << "We are using Table 2 for CMSZDIFF12.h" << endl;
  }
  string line;

  int idat;
  double idum, jdum, kdum;
  double ptmax, ptmin;
  double etamax, etamin;

  const int fNDataRap = 5;                  // Number of rapidity bins
  const int fNDataPT = 10;		              // Number of pt bins
  const double s = 8000;

  // Reading data
  for (int i = 0; i < fNDataPT; i++)
  {
    getline(f1,line);
    istringstream lstream(line);

    lstream >> ptmin >> ptmax;

    for (int j = 0; j < fNDataRap; j++)
    {
      idat = i + j*fNDataPT;		       // Structure for Covariance matrix

      etamin = 0.0 + j * 0.4;
      etamax = 0.4 + j * 0.4;

      fKin1[idat] = (etamin + etamax) / 2.0;   // average rapidity
      fKin2[idat] = (ptmin + ptmax) / 2.0;     // average pt
      fKin2[idat] *= fKin2[idat];
      fKin3[idat] = s;                         // only eta and pt relevant

      lstream >> fData[idat];
      fData[idat] *= 1000.;                    // Make pb -> fb

      lstream >> idum;
      fStat[idat] = 0.0;
      fStat[idat] *= 1000.;                   // Make pb -> fb
    }
  }
  // Reading covmat
  double** covmat = new double*[fNData];
  for(int i = 0; i < fNData; i++)
  {
    covmat[i] = new double[fNData];

    getline(f2,line);
    istringstream lstream(line);
    lstream >> idum >> jdum >> kdum;

    for(int j = 0; j < fNData; j++)
    {
      lstream >> covmat[i][j];
      covmat[i][j] *= 1.0e6;
      //cout << i << " " << j << " " << covmat[i][j] << endl;
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
    {
      for (int l = 0; l < fNSys - 1; l++)
      {
        fSys[i][l].add  = syscor[i][l];
        fSys[i][l].mult = fSys[i][l].add*100/fData[i];
        fSys[i][l].type = ADD;
        fSys[i][l].name = "CORR";
      }

      // Luminosity Uncertainty
      // CMS Luminosity Uncertainty, 2012 data set: 2.6%
      // arXiv:1504.03511v2
        fSys[i][fNSys-1].mult = 2.6;
        fSys[i][fNSys-1].add  = fData[i]*fSys[i][fNSys-1].mult/100;
        fSys[i][fNSys-1].type = MULT;
        fSys[i][fNSys-1].name = "CMSLUMI12";
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
