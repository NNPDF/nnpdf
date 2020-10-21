/**
 * EMCF2c1987.cc
 *
 * EMC 1987
 * Charm Production in Deep Inelastic Muon-Iron Interactions at 200 GeV/c
 * Z. Phys. C - Particles and Fields 35, 1-6 (1987)
 * http://hepdata.cedar.ac.uk/view/ins230629
 * 
 */

#include "EMCF2c1987.h"

void EMCF2c1987Filter::ReadData()
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
    double alpha = 1./137.035999074;
    double resc = (pow(Q2,2)*x)/(2.*pow(alpha,2)*M_PI*Yp);
    double convfact = 1e-6*(1./0.3894); //nb to GeV^-2 conv 

    double eps = pow(1.+(Q2+nu*nu)/(2.*Ein*(Ein-nu)-Q2/2.),-1);
    double Gamma = alpha*(nu-Q2/(2.*Mp))/(2.*M_PI*Q2*Ein*Ein*(1-eps));

    fData[i] = Gamma*resc*Jac*convfact*data;
    fStat[i] = Gamma*resc*Jac*convfact*stat;
   
   //rescaling for BR - check
    fData[i] = fData[i]*0.82;
    double sist = fData[i]*0.15;
    fSys[i][0].add = sist;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "CORR_EMC";
    fSys[i][0].mult = fSys[i][0].add/(fData[i]*1e-2);
  
  }

  
  f1.close();

}