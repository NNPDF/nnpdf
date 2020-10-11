/*
   EIC pseudodata for the Yellow Report.
   The pseudodata correspond to CC and NC electron-proton reduced cross 
   sections. The kinematic bins and the associated statistical and systematic 
   uncertainties have been provided by the IR EIC WG and correspond, for both 
   CC and NC data, to a pessimistic and an optimistic scenario. For details, see
   https://github.com/JeffersonLab/txgrids
   The central value of the pseudodata is geenrated by flucutating the
   theoretical predictions obtained from two different PDF set:
   - 200609-ern-001 for a NNPDF31-like proton fit
   - 150719-NNPDF31-pch-nonuc for a proton fit used as baseline for nuclear PDFs
*/

#include "EIC.h"


/*
PROTON PSEUDODATA
*/

//Charged-current measurements, optimistic scenario
void EIC_CC_140_OPTFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/CCemP-optimistic.dat";
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
      lstream >> idum >> sqrts >> fKin1[i] >> fKin2[i] // x, Q2
	      >> fData[i]
	      >> fStat[i]
	      >> fSys[i][0].mult
	      >> fSys[i][1].mult;

      fStat[i] = fStat[i]/100. * fData[i]; 
      
      fKin3[i] = fKin2[i] / fKin1[i] / sqrts / sqrts; //y
      fSys[i][0].add  = fSys[i][0].mult / 100. * fData[i];
      fSys[i][0].type = ADD;
      fSys[i][0].name = "CORR";
      fSys[i][1].add  = fSys[i][1].mult / 100. * fData[i];
      fSys[i][1].type = MULT;
      fSys[i][1].name = "CORR";      
    }

  f1.close();
  
}

//Charged-current measurements, pessimistic scenario
void EIC_CC_140_PESFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/CCemP-pessimistic.dat";
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
      lstream >> idum >> sqrts >> fKin1[i] >> fKin2[i] // x, Q2
	      >> fData[i]
	      >> fStat[i]
	      >> fSys[i][0].mult
	      >> fSys[i][1].mult;
      
      fStat[i] = fStat[i]/100. * fData[i]; 
      
      fKin3[i] = fKin2[i] / fKin1[i] / sqrts / sqrts; //y
      fSys[i][0].add  = fSys[i][0].mult / 100. * fData[i];
      fSys[i][0].type = ADD;
      fSys[i][0].name = "CORR";
      fSys[i][1].add  = fSys[i][1].mult / 100. * fData[i];
      fSys[i][1].type = MULT;
      fSys[i][1].name = "CORR";
    }

  f1.close();
  
}

//Neutral-current measurements, optimistic scenario
void EIC_NC_140_OPTFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemP-optimistic_140.dat";
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
      lstream >> idum >> sqrts >> fKin1[i] >> fKin2[i] // x, Q2
	      >> fData[i]
	      >> fStat[i]
	      >> fSys[i][0].mult
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

void EIC_NC_63_OPTFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemP-optimistic_63.dat";
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
      lstream >> idum >> sqrts >> fKin1[i] >> fKin2[i] // x, Q2
	      >> fData[i]
	      >> fStat[i]
	      >> fSys[i][0].mult
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

void EIC_NC_44_OPTFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemP-optimistic_44.dat";
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
      lstream >> idum >> sqrts >> fKin1[i] >> fKin2[i] // x, Q2
	      >> fData[i]
	      >> fStat[i]
	      >> fSys[i][0].mult
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

void EIC_NC_28_OPTFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemP-optimistic_28.dat";
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
      lstream >> idum >> sqrts >> fKin1[i] >> fKin2[i] // x, Q2
	      >> fData[i]
	      >> fStat[i]
	      >> fSys[i][0].mult
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

//Neutral-current measurements, pessimistic scenario
void EIC_NC_140_PESFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemP-pessimistic_140.dat";
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
      lstream >> idum >> sqrts >> fKin1[i] >> fKin2[i] // x, Q2
	      >> fData[i]
	      >> fStat[i]
	      >> fSys[i][0].mult
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

void EIC_NC_63_PESFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemP-pessimistic_63.dat";
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
      lstream >> idum >> sqrts >> fKin1[i] >> fKin2[i] // x, Q2
	      >> fData[i]
	      >> fStat[i]
	      >> fSys[i][0].mult
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

void EIC_NC_44_PESFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemP-pessimistic_44.dat";
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
      lstream >> idum >> sqrts >> fKin1[i] >> fKin2[i] // x, Q2
	      >> fData[i]
	      >> fStat[i]
	      >> fSys[i][0].mult
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

void EIC_NC_28_PESFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemP-pessimistic_28.dat";
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
      lstream >> idum >> sqrts >> fKin1[i] >> fKin2[i] // x, Q2
	      >> fData[i]
	      >> fStat[i]
	      >> fSys[i][0].mult
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

/*
NUCLEAR PSEUDODATA
*/


//Charged-current measurements, optimistic scenario
void EIC_CC_140_OPT_NUCLFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/CCemP-optimistic_nucl.dat";
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
      lstream >> idum >> sqrts >> fKin1[i] >> fKin2[i] // x, Q2
	      >> fData[i]
	      >> fStat[i]
	      >> fSys[i][0].mult
	      >> fSys[i][1].mult;

      fStat[i] = fStat[i]/100. * fData[i]; 
      
      fKin3[i] = fKin2[i] / fKin1[i] / sqrts / sqrts; //y
      fSys[i][0].add  = fSys[i][0].mult / 100. * fData[i];
      fSys[i][0].type = ADD;
      fSys[i][0].name = "CORR";
      fSys[i][1].add  = fSys[i][1].mult / 100. * fData[i];
      fSys[i][1].type = MULT;
      fSys[i][1].name = "CORR";      
    }

  f1.close();
  
}

//Charged-current measurements, pessimistic scenario
void EIC_CC_140_PES_NUCLFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/CCemP-pessimistic_nucl.dat";
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
      lstream >> idum >> sqrts >> fKin1[i] >> fKin2[i] // x, Q2
	      >> fData[i]
	      >> fStat[i]
	      >> fSys[i][0].mult
	      >> fSys[i][1].mult;
      
      fStat[i] = fStat[i]/100. * fData[i]; 
      
      fKin3[i] = fKin2[i] / fKin1[i] / sqrts / sqrts; //y
      fSys[i][0].add  = fSys[i][0].mult / 100. * fData[i];
      fSys[i][0].type = ADD;
      fSys[i][0].name = "CORR";
      fSys[i][1].add  = fSys[i][1].mult / 100. * fData[i];
      fSys[i][1].type = MULT;
      fSys[i][1].name = "CORR";
    }

  f1.close();
  
}

//Neutral-current measurements, optimistic scenario
void EIC_NC_140_OPT_NUCLFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemP-optimistic_nucl_140_u.dat";
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
      lstream >> idum >> sqrts >> fKin1[i] >> fKin2[i] // x, Q2
	      >> fData[i]
	      >> fStat[i]
	      >> fSys[i][0].mult
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

void EIC_NC_63_OPT_NUCLFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemP-optimistic_nucl_63_u.dat";
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
      lstream >> idum >> sqrts >> fKin1[i] >> fKin2[i] // x, Q2
	      >> fData[i]
	      >> fStat[i]
	      >> fSys[i][0].mult
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

void EIC_NC_44_OPT_NUCLFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemP-optimistic_nucl_44_u.dat";
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
      lstream >> idum >> sqrts >> fKin1[i] >> fKin2[i] // x, Q2
	      >> fData[i]
	      >> fStat[i]
	      >> fSys[i][0].mult
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

void EIC_NC_28_OPT_NUCLFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemP-optimistic_nucl_28_u.dat";
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
      lstream >> idum >> sqrts >> fKin1[i] >> fKin2[i] // x, Q2
	      >> fData[i]
	      >> fStat[i]
	      >> fSys[i][0].mult
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

//Neutral-current measurements, pessimistic scenario
void EIC_NC_140_PES_NUCLFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemP-pessimistic_nucl_140_u.dat";
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
      lstream >> idum >> sqrts >> fKin1[i] >> fKin2[i] // x, Q2
	      >> fData[i]
	      >> fStat[i]
	      >> fSys[i][0].mult
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

void EIC_NC_63_PES_NUCLFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemP-pessimistic_nucl_63_u.dat";
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
      lstream >> idum >> sqrts >> fKin1[i] >> fKin2[i] // x, Q2
	      >> fData[i]
	      >> fStat[i]
	      >> fSys[i][0].mult
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

void EIC_NC_44_PES_NUCLFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemP-pessimistic_nucl_44_u.dat";
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
      lstream >> idum >> sqrts >> fKin1[i] >> fKin2[i] // x, Q2
	      >> fData[i]
	      >> fStat[i]
	      >> fSys[i][0].mult
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

void EIC_NC_28_PES_NUCLFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemP-pessimistic_nucl_28_u.dat";
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
      lstream >> idum >> sqrts >> fKin1[i] >> fKin2[i] // x, Q2
	      >> fData[i]
	      >> fStat[i]
	      >> fSys[i][0].mult
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
