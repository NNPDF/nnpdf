/*
   EIC pseudodata for the Yellow Report.
   The pseudodata correspond to the CC and NC cross sections, evaluated from the
   200609-ern-001 fit. The kinematic bins and the associated statistical and
   systematic uncertainties have been provided by the IR EIC WG and correspond,
   for both CC and NC data, to a pessimistic and an optimistic scenario. See
   https://github.com/JeffersonLab/txgrids
*/

#include "EIC.h"

void EICFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/" + fSetName + ".csv";
  f1.open(datafile.str().c_str(), ios::in);
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  string line;
  getline(f1,line);

  for(int i=0; i<fNData; i++)
    {
      getline(f1,line);
      istringstream lstream(line);
      double sqrts;
      int idum;
      char comma;
      lstream >> idum            >> comma
	      >> sqrts           >> comma
	      >> fKin1[i]        >> comma
	      >> fKin2[i]        >> comma // x, Q2
	      >> fData[i]        >> comma
	      >> fStat[i]        >> comma
	      >> fSys[i][0].mult >> comma
	      >> fSys[i][1].mult;

      fStat[i] = fStat[i]/100. * fData[i]; 
      
      fKin3[i] = fKin2[i] / fKin1[i] / sqrts / sqrts; //y
      fSys[i][0].add  = fSys[i][0].mult / 100. * fData[i];
      fSys[i][0].type = ADD;
      fSys[i][0].name = "UNCORR";
      fSys[i][1].add  = fSys[i][1].mult / 100. * fData[i];
      fSys[i][1].type = MULT;
      fSys[i][1].name = "CORR";      
    }

  f1.close();
  
}
