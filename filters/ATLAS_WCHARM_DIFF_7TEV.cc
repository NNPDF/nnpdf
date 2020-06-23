/*
*   Experiment: ATLAS
*   Process: proton + proton -> W- + charm
*   Center of mass energy = 7 TeV
*   Intergrated luminosity  = 4.6 1/fb
*   arXiv:1402.6263
*   Published in JHEP 05 (2014) 068
*   HEPData: https://www.hepdata.net/record/ins1282447
*   Description: The production of a W boson in association with a single charm 
*   quark is studied using 4.6 fb^-1 of pp collision data at sqrt(s) = 7 TeV 
*   collected in 2011 with the ATLAS detector at the Large Hadron Collider. 
*   In events in which a W boson decays to an electron or muon, the charm quark *   is tagged either by its semileptonic decay to a muon or by the presence of 
*   a charmed meson. The integrated and differential cross sections as a 
*   function of the pseudorapidity of the lepton from the W-boson decay are 
*   measured. Results are compared to the predictions of next-to-leading-order 
*   QCD calculations obtained from various parton distribution function 
*   parameterisations. The ratio of the strange-to-down sea quark distributions 
*   is determined to be 0.96 +0.26,-0.30 at Q^2 = 1.9 GeV^2 which supports the 
*   hypothesis of an SU(3) symmetric composition of the light quark sea. 
*   Additionally, the cross-section ratio sigma(W^+ + cbar)/sigma(W^- + c) is 
*   compared to the predictions obtained using parton distribution function 
*   parameterisations with different assumptions on the s-sbar asymmetry.
* 
*   The relevant data is taken from Table 3 (central value and statistical
*   uncertainty) and Table 12 and Table 13 (breakdown of systematic 
*   uncertainties), respectively for W+ +cbar and W- + c.
*
*   Systematici uncertainties are supplemented by a theoretical uncertainty
*   that takes into account missing NNLO corrections in the matrix element.
*   This uncertainty was estimated by means as the asymmetric envelope of the 
*   3pt renormalisation scale variation (Eq. 4.16 in 1906.10698).
*/

#include "ATLAS_WCHARM_WP_DIFF_7TEV.h"
void ATLAS_WCHARM_WP_DIFF_7TEVFilter::ReadData()
{
  // Opening files
  fstream f1;
  stringstream datafileWP("");
  datafileWP << dataPath()
       << "rawdata/" << fSetName << "/" << fSetName << ".data";
  f1.open(datafileWP.str().c_str(), ios::in);
  if (f1.fail())
  {
    cerr << "Error opening data file " << datafileWP.str() << endl;
    exit(-1);
  }

  //Starting filter
  string line;
  double etamin, etamax;       //rapidity binning
  double MW2 = pow(MW, 2.0);   //W mass
  double s = 7000;             //LHC at 7TeV

  for(int i=0; i<fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    //Reading in an interpretation of each column
    lstream >> etamin >> etamax >> fData[i] >> fStat[i];

    fData[i] = fData[i]*1000; // changing pb to fb for APPLgrid
    fStat[i] = fStat[i]*1000; // changing pb to fb for APPLgrid
    //Defining the kinematic variables
    fKin1[i] = (etamin + etamax)*0.5;    // eta
    fKin2[i] = MW2;                      // Mass W squared
    fKin3[i] = s;                        // sqrt(s)

    //Reading in the systematics
    for(int k=0; k<fNSys; k++)
    {
      lstream >> fSys[i][k].mult;
      if(k==fNSys-2||k==fNSys-1)
	fSys[i][k].mult /= sqrt(2.);
      fSys[i][k].type = MULT;
      fSys[i][k].add = fSys[i][k].mult*fData[i]/100;
      if(k == 0)
        fSys[i][k].name = "UNCORR";
      else if(k == fNSys - 3)
        fSys[i][k].name = "ATLASLUMI11";
      else if(k == fNSys - 2)
        fSys[i][k].name = "SKIP";
      else if(k == fNSys - 1)
        fSys[i][k].name = "SKIP"; 
      else
        fSys[i][k].name = "CORR";
    }
  }
  f1.close();
}

#include "ATLAS_WCHARM_WM_DIFF_7TEV.h"
void ATLAS_WCHARM_WM_DIFF_7TEVFilter::ReadData()
{
  // Opening files
  fstream f1;
  stringstream datafileWM("");
  datafileWM << dataPath()
       << "rawdata/" << fSetName << "/" << fSetName << ".data";
  f1.open(datafileWM.str().c_str(), ios::in);
  if (f1.fail())
  {
    cerr << "Error opening data file " << datafileWM.str() << endl;
    exit(-1);
  }

  //Starting filter
  string line;
  double etamin, etamax;       //rapidity binning
  double MW2 = pow(MW, 2.0);   //W mass
  double s = 7000;             //LHC at 7TeV [GeV]

  for(int i=0; i<fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    //Reading in an interpretation of each column
    lstream >> etamin >> etamax >> fData[i] >> fStat[i];

    fData[i] = fData[i]*1000; // changing pb to fb for APPLgrid
    fStat[i] = fStat[i]*1000; // changing pb to fb for APPLgrid
    //Defining the kinematic variables
    fKin1[i] = (etamin + etamax)*0.5;    // eta
    fKin2[i] = MW2;                      // Mass W squared
    fKin3[i] = s;                        // sqrt(s)

    //Reading in the systematics
    for(int k=0; k<fNSys; k++)
    {
      lstream >> fSys[i][k].mult;
      if(k==fNSys-2||k==fNSys-1)
	fSys[i][k].mult /= sqrt(2.);
      fSys[i][k].type = MULT;
      fSys[i][k].add = fSys[i][k].mult*fData[i]/100;
      if(k == 0)
        fSys[i][k].name = "UNCORR";
      else if(k == fNSys - 3)
        fSys[i][k].name = "ATLASLUMI11";
      else if(k == fNSys - 2)
        fSys[i][k].name = "SKIP";
      else if(k == fNSys - 1)
        fSys[i][k].name = "SKIP"; 
      else
        fSys[i][k].name = "CORR";
    }
  }
  f1.close();
}

