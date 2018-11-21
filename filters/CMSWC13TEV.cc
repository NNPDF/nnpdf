/*
*   Experiment: CMS 
*   Process: proton + proton -> W + charm
*   Center of mass energy = 13 TeV
*   Intergrated luminosity  = 35.7 1/fb
*   Reference: https://cds.cern.ch/record/2314570/files/SMP-17-014-pas.pdf
*/

#include "CMSWC13TEV.h"

void CMSWC13TEVFilter::ReadData()
{
  // Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() 
       << "rawdata/" << fSetName << "/" << fSetName << ".data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  //Starting filter
  string line;
  double etamin, etamax;        //eta ranges
  double MW2 = pow(MW,2.0);     //W mass
  double s = 13000;             //LHC at 13TeV
  double stmp, dtmp;
  //Systematics
  double fLumiP, fLumiM;
  double fTrackingP, fTrackingM;
  double fBranchingP, fBranchingM;
  double fMuonsP, fMuonsM;
  double fNselP, fNselM;
  double fKinematicsP, fKinematicsM;
  double fNormalizationP, fNormalizationM;
  double fPmissP, fPmissM;
  double fPileUpP, fPileUpM;
  double fPDFP, fPDFM;
  double fSecP, fSecM;
  double fFragmentationP, fFragmentationM;
  double fMonteCarloP, fMonteCarloM;

  for(int i=0; i<fNData;i++)
    {
      getline(f1,line);
      istringstream lstream(line);

      //Reading in an interpretation of each column
      lstream >> etamin >> etamax >> fData[i] >> fStat[i]
      >> fLumiP >> fLumiM
      >> fTrackingP >> fTrackingM
      >> fBranchingP >> fBranchingM
      >> fMuonsP >> fMuonsM
      >> fNselP >> fNselM
      >> fKinematicsP >> fKinematicsM
      >> fNormalizationP >> fNormalizationM
      >> fPmissP >> fPmissM
      >> fPileUpP >> fPileUpM
      >> fPDFP >> fPDFM
      >> fSecP >> fSecM
      >> fFragmentationP >> fFragmentationM
      >> fMonteCarloP >> fMonteCarloM;

      //Defining the kinematic variables
      fKin1[i] = (etamax + etamin)*0.5;    // eta
      fKin2[i] = MW2;                      // Mass W squared
      fKin3[i] = s;                        // sqrt(s)

      //Declaring the systematic uncertainties and symmetrising

      //Luminosity
      symmetriseErrors(fLumiP,fLumiM,&stmp,&dtmp);
      fSys[i][0].mult=stmp;
      fSys[i][0].type = MULT;
      fSys[i][0].name = "CORR";

      //Tracking
      symmetriseErrors(fTrackingP,fTrackingM,&stmp,&dtmp);
      fSys[i][1].mult=stmp;
      fSys[i][1].type = MULT;
      fSys[i][1].name = "CORR";

      //Branching
      symmetriseErrors(fBranchingP,fBranchingM,&stmp,&dtmp);
      fSys[i][2].mult=stmp;
      fSys[i][2].type = MULT;
      fSys[i][2].name = "CORR";

      //Muons
      symmetriseErrors(fMuonsP,fMuonsM,&stmp,&dtmp);
      fSys[i][3].mult=stmp;
      fSys[i][3].type = MULT;
      fSys[i][3].name = "CORR";

      //Nsel determination
      symmetriseErrors(fNselP,fNselM,&stmp,&dtmp);
      fSys[i][4].mult=stmp;
      fSys[i][4].type = MULT;
      fSys[i][4].name = "CORR";

      //D* kinematics
      symmetriseErrors(fKinematicsP,fKinematicsM,&stmp,&dtmp);
      fSys[i][5].mult=stmp;
      fSys[i][5].type = MULT;
      fSys[i][5].name = "CORR";

      //Bg normalization
      symmetriseErrors(fNormalizationP,fNormalizationM,&stmp,&dtmp);
      fSys[i][6].mult=stmp;
      fSys[i][6].type = MULT;
      fSys[i][6].name = "CORR";

      //P_T^miss
      symmetriseErrors(fPmissP,fPmissM,&stmp,&dtmp);
      fSys[i][7].mult=stmp;
      fSys[i][7].type = MULT;
      fSys[i][7].name = "CORR";

      //Pile up
      symmetriseErrors(fPileUpP,fPileUpM,&stmp,&dtmp);
      fSys[i][8].mult=stmp;
      fSys[i][8].type = MULT;
      fSys[i][8].name = "CORR";

      //PDF
      symmetriseErrors(fPDFP,fPDFM,&stmp,&dtmp);
      fSys[i][9].mult=stmp;
      fSys[i][9].type = MULT;
      fSys[i][9].name = "CORR";


      //Sec. Vtx.
      symmetriseErrors(fSecP,fSecM,&stmp,&dtmp);
      fSys[i][10].mult=stmp;
      fSys[i][10].type = MULT;
      fSys[i][10].name = "CORR";

      //Sec. Vtx.
      symmetriseErrors(fFragmentationP,fFragmentationM,&stmp,&dtmp);
      fSys[i][11].mult=stmp;
      fSys[i][11].type = MULT;
      fSys[i][11].name = "CORR";

      //Sec. Vtx.
      symmetriseErrors(fMonteCarloP,fMonteCarloM,&stmp,&dtmp);
      fSys[i][12].mult=stmp;
      fSys[i][12].type = MULT;
      fSys[i][12].name = "CORR";
      
      // Convert mutltiplicative uncertainties to additive
      for (int l = 0; l < fNSys; l++)
        {
          fSys[i][l].add = fSys[i][l].mult*fData[i]*1e-2;
        }
      
    }  

  f1.close();

}