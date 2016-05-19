/**
 * EMC.cc
 *
 * F2P
 * A Detailed Study of the Proton Structure Functions in Deep Inelastic Muon - Proton Scattering
 * Nucl.Phys. B259 (1985) 189, 1985
 * http://dx.doi.org/10.17182/hepdata.13916
 * From Table 237 onwards
 *
 * F2D
 * Measurements of the Nucleon Structure Functions F(2)N in Deep Inelastic
 * Muon Scattering from Deuterium and Comparison with Those from Hydrogen and Iron
 * Nucl.Phys. B293 (1987) 740, 1987
 * http://dx.doi.org/10.17182/hepdata.33821
 *
 * EMC 1987
 * Charm Production in Deep Inelastic Muon-Iron Interactions at 200 GeV/c
 * Z. Phys. C - Particles and Fields 35, 1-6 (1987)
 * http://hepdata.cedar.ac.uk/view/ins230629
 * 
 */

#include "EMC.h"

void EMCF2PFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/EMCF2P.data";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Reading data
  string line;

  //skip first two lines
  getline(f1,line);
  getline(f1,line);
    

  //Filtering data
  for (int i = 0; i < fNData; i++)
  {
    f1 >> fKin1[i]; //x
    f1 >> fKin2[i]; //Q2
    
    fKin3[i] = 0.0; //y

    f1 >> fData[i]; //obs
 
    double stat = 0;
    double sist = 0;
    double dummy;

    f1 >> stat;
    f1 >> dummy;

    fStat[i] = stat;

    f1 >> sist;
    f1 >> dummy;

    //check here
    fSys[i][0].add = sist;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CORR";

    fSys[i][0].mult = fSys[i][0].add/(fData[i]*1e-2);

  }

  
  f1.close();

}

void EMCF2DFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/EMCF2D.data";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Reading data
  string line;

  //skip first two lines
  getline(f1,line);
  getline(f1,line);
    

  //Filtering data
  for (int i = 0; i < fNData; i++)
  {
    f1 >> fKin1[i]; //x
    f1 >> fKin2[i]; //Q2
    
    fKin3[i] = 0.0; //y

    f1 >> fData[i]; //obs

    double stat = 0;
    double sist = 0;
    double dummy;

    f1 >> stat;
    f1 >> dummy;

    fStat[i] = stat;

    f1 >> sist;
    f1 >> dummy;

    //check here
    fSys[i][0].add = sist;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CORR";

    fSys[i][0].mult = fSys[i][0].add/(fData[i]*1e-2);

  }

  
  f1.close();

}


void EMCFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/EMC.data";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Reading data
  string line;

  //skip first line
  getline(f1,line);
    
  //Filtering data
  for (int i = 0; i < fNData; i++)
  {
    double Q2l;
    double Q2h;
    double nu;
    double dummy;
    double Q2;

    f1 >> Q2l;
    f1 >> Q2h;


    Q2 = (Q2l+Q2h)/2.;
    //data are integrated over each bin
    double deltaQ2 = Q2h-Q2l;
    double deltanu = 20.; //20 GeV
    fKin2[i] = Q2; //Q2 

    f1 >> nu;
    f1 >> dummy;
    f1 >> dummy;

    double Mp = 0.938272046;  //Proton Mass (GeV)
    double x = Q2/(2.*Mp*nu);

    fKin1[i] = x; //x

    double Ein = 200.; //Ebeam = 200 GeV
    double y = nu/Ein;

    fKin3[i] = y; //y

    double data;
    f1 >> data; //obs

    double stat = 0;

    f1 >> stat;
    f1 >> dummy;

    double Yp = (1.+pow(1.-y,2));
    double Jac = Q2/(2*pow(x,2)*Mp);
    double as = 1./137.035999074;
    double resc = (pow(Q2,2)*x)/(2.*pow(as,2)*M_PI*Yp);
    double convfact = 1e-6*(1./0.3894); //nb to GeV^-2 conv 

    fData[i] = resc*Jac*convfact*data/deltanu/deltaQ2;
    fStat[i] = resc*Jac*convfact*stat/deltanu/deltaQ2;
   
    /*
    fSys[i][0].add = sist;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CORR";

    fSys[i][0].mult = fSys[i][0].add/(fData[i]*1e-2);
    */
  
  }

  
  f1.close();

}