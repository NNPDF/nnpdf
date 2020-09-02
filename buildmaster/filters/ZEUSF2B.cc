/**
 * ZEUS F2B
 * 
 * 
 */

#include "ZEUSF2B.h"

void ZEUSF2BFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/ZEUSHERAF2B.data";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Reading data
  string line;
  getline(f1,line);

  //Filtering data
  for (int i = 0; i < fNData; i++)
  {
    double dummy;
    double data;

    f1 >> fKin2[i]; //Q2
    f1 >> dummy >> dummy;
    f1 >> fKin1[i]; //x
    f1 >> dummy >> dummy;
    
    double E= 27.5;    // GeV
    double mp = 0.938; //GeV
    fKin3[i] = fKin2[i] / ( 2.*E*mp*fKin1[i] ); //y

    f1 >> data; //obs
 
    double statr, statl, sistr, sistl;
    double stat, sist;
    double delta = 0;

    f1 >> statr;
    f1 >> statl;
    symmetriseErrors(statr,statl,&stat,&delta);

    fStat[i] = stat;    
    data += delta;

    for (int l = 0; l < fNSys; ++l)
    {
      f1 >> sistr;
      f1 >> sistl;

      symmetriseErrors(sistr,sistl,&sist,&delta);

      data += delta;

      fSys[i][l].add = sist;
      fSys[i][l].type = MULT;
      fSys[i][l].name = "CORR";
      //check whether the extr uncertainty is correlated - unclear in the paper

    }

    fData[i] = data;

    for (int l = 0; l < fNSys; l++)
      fSys[i][l].mult = fSys[i][l].add/(fData[i]*1e-2);


  }

  
  f1.close();

}
