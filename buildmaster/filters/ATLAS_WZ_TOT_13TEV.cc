/*
This implementation supersedes that in ATLASWZTOT13TEV81PB
in that it includes the experimental correlation coefficients
between W+, W- and Z cross sections.
*/


#include "ATLAS_WZ_TOT_13TEV.h"

void ATLAS_WZ_TOT_13TEVFilter::ReadData()
{
  fstream f1, f2;
  
  stringstream datafile(""), corrfile("");
  datafile << dataPath()
	   << "rawdata/ATLASWZTOT13TEV81PB/data.txt";
  corrfile << dataPath()
	   << "rawdata/ATLASWZTOT13TEV81PB/corr.txt";

  f1.open(datafile.str().c_str(), ios::in);
  if (f1.fail())
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }

  f2.open(corrfile.str().c_str(), ios::in);
  if (f2.fail())
    {
      cerr << "Error opening data file " << corrfile.str() << endl;
      exit(-1);
    }

  //Central value, statistical and systematic uncertainties
   string line;
  for(int i=0; i<11; i++)
    {
      getline(f1,line);
    }

  vector<double> sys(fNData);
  vector<double> lumi(fNData);

  for(int i=0; i<fNData; i++)
    {
      string sdum;

      getline(f1,line);
      istringstream lstream(line);
      lstream >> sdum 
	      >> fData[i]
	      >> fStat[i]
	      >> sys[i]
	      >> lumi[i]
	      >> fKin1[i]
	      >> fKin2[i]
	      >> fKin3[i];
      //conversion to femtobarns
      fData[i] *= 1e6;
      fStat[i] *= 1e6;
      sys[i]   *= 1e6;
      lumi[i]  *= 1e6;
    }
  
  //Correlation coefficients
  for(int i=0; i<25; i++)
    {
       getline(f2,line);
    }

  double** covmat = new double*[fNData];
  for(int i=0; i<fNData; i++)
    {
      covmat[i] = new double[fNData];
      for(int j=0; j<fNData; j++)
	{
	  covmat[i][j] = 1.;
	}
    }

  char comma;
  double ddum;
  
  vector<double> corrcoeff(3);

  for(int k=0; k<3; k++)
    {
      getline(f2,line);
      istringstream lstream(line);
      lstream >> ddum >> comma >> ddum >> comma >> corrcoeff[k];
    }

  covmat[0][1] = corrcoeff[0];
  covmat[1][0] = covmat[0][1];
  covmat[0][2] = corrcoeff[1];
  covmat[2][0] = covmat[0][2];
  covmat[1][2] = corrcoeff[2];
  covmat[2][1] = covmat[1][2];

  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  covmat[i][j] = covmat[i][j]*sys[i]*sys[j];
	}
    }

  //Generate artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
    {
      throw runtime_error("Couldn't generate artificial systematics for " + fSetName);
    }
  
  for(int i=0; i<fNData; i++)
    {
      for(int j=0; j<fNData; j++)
	{
	  fSys[i][j].add  = syscor[i][j];
	  fSys[i][j].mult = fSys[i][j].add*1e2/fData[i];
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = "CORR";
	}
      fSys[i][fNSys-1].add = lumi[i];
      fSys[i][fNSys-1].mult = fSys[i][fNSys-1].add*1e2/fData[i];
      fSys[i][fNSys-1].type = MULT;
      fSys[i][fNSys-1].name = "ATLASLUMI13";
    }

  f1.close();
  f2.close();

  for(int i = 0; i < fNData; i++) 
    {
      delete[] syscor[i];
      delete[] covmat[i];
    }
  delete[] syscor;
  delete[] covmat; 
  
}
