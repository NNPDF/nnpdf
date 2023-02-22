/*
Reference: Phys.Rev.Lett. 121 (2018) 6, 062002
arXiv: [1805.04736]
HepData: https://www.hepdata.net/record/ins1672941
Description: The pseudorapidity distributions of dijets as a function of their 
average transverse momentum are measured in proton-proton (pp) collisions. 
The data samples were collected by the CMS experiment at the CERN LHC, at a 
nucleon-nucleon center-of-mass energy of 5.02 TeV. 
The data are taken from Tabs. 3a, 3b, 3c, 3d and 3e of the HepData entry. No
information on correlations is provided.
*/

#include "CMS_2JET_5TEV.h"

void CMS_2JET_5TEVFilter::ReadData()
{

  //These are the bin limits. The data is in the middle
  std::vector<double> pTavg_bins = {55., 75., 95., 115., 150., 400.};
  int n = 0;
  for (size_t ibin = 1; ibin < pTavg_bins.size(); ibin++)
  {
    // Opening files
    fstream f;

    //std::size_t found = fSetName.find("bin");
    //int bin = std::stoi(fSetName.substr(found+3, 1));

    // rapidity distribution
    stringstream DataFile("");
    DataFile << dataPath() << "rawdata/CMS_2JET_5TEV/CMS_2JET_5TEV_bin" + std::to_string(ibin) + ".csv";
    f.open(DataFile.str().c_str(), ios::in);

    if (f.fail())
    {
      cerr << "Error opening data file " << DataFile.str() << endl;
      exit(-1);
    }

    // Starting filter
    double sqrt_s = 5020; // GeV
    string line;

    std::vector<double> totsys(fNData);
    double sys, yavg, ymin, ymax, dum;

    // read the lines containing comments and the total cross section (not used in this distribution)
    for (int i = 0; i < 13; i++)
      getline(f, line);

    while (getline(f, line))
    {
      istringstream lstream(line); // the line is split in its columns, so that f reads a single value in each go

      //Already in CM (from HepData)
      lstream >> yavg >> ymin >> ymax;
      fKin1[n] = yavg;                                                                          // dijet rapidity
      fKin2[n] = pow(pTavg_bins[ibin - 1] + (pTavg_bins[ibin] - pTavg_bins[ibin - 1]) / 2., 2); // dijet <pTavg1,2>^2
      fKin3[n] = sqrt_s;                                                                        // sqrt(s)

      lstream >> fData[n];
      lstream >> fStat[n];
      lstream >> dum;
      lstream >> sys;

      /*
	The total systematic uncertainties in ηdijet and in the ratios of
	the pPb and pp spectra are evaluated by summing in quadrature over
	the contributions from the above sources: the uncertainties due to
	the JES, jet energy resolution, and jet angular resolution 
      */
      fSys[n][0].add = sys;
      fSys[n][0].mult = (fSys[n][0].add / fData[n]) * 100;
      fSys[n][0].type = ADD;
      fSys[n][0].name = "UNCORR";

      n++;
    }
  }
}
