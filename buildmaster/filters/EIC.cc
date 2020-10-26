/*
   EIC pseudodata for the Yellow Report.
   The pseudodata correspond to the CC and NC cross sections, evaluated from the
   200609-ern-001 fit. The kinematic bins and the associated statistical and
   systematic uncertainties have been provided by the IR EIC WG and correspond,
   for both CC and NC data, to a pessimistic and an optimistic scenario. See
   https://github.com/JeffersonLab/txgrids
*/

#include "EIC.h"

//Charged-current measurements, optimistic scenario

//A) electron-proton
void EIC_CC_EMP_140_OPTFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/CCemP-optimistic.csv";
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

//B) positron-proton
void EIC_CC_EPP_140_OPTFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/CCepP-optimistic.csv";
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

//Charged-current measurements, pessimistic scenario

//A) electron-proton
void EIC_CC_EMP_140_PESFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/CCemP-pessimistic.csv";
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

//B) positron-proton
void EIC_CC_EPP_140_PESFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/CCepP-pessimistic.csv";
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

//Neutral-current measurements, optimistic scenario

//A) electron-proton
void EIC_NC_EMP_140_OPTFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemP-optimistic_140.csv";
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

void EIC_NC_EMP_63_OPTFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemP-optimistic_63.csv";
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

void EIC_NC_EMP_44_OPTFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemP-optimistic_44.csv";
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

void EIC_NC_EMP_28_OPTFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemP-optimistic_28.csv";
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

//B) positron-proton
void EIC_NC_EPP_140_OPTFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCepP-optimistic_140.csv";
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

void EIC_NC_EPP_63_OPTFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCepP-optimistic_63.csv";
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

void EIC_NC_EPP_44_OPTFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCepP-optimistic_44.csv";
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

void EIC_NC_EPP_28_OPTFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCepP-optimistic_28.csv";
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

//C) electron-deuteron

void EIC_NC_EMD_88_OPTFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemD-optimistic_88.csv";
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

void EIC_NC_EMD_66_OPTFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemD-optimistic_66.csv";
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

void EIC_NC_EMD_28_OPTFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemD-optimistic_28.csv";
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

//D) positron-deuteron

void EIC_NC_EPD_88_OPTFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCepD-optimistic_88.csv";
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

void EIC_NC_EPD_66_OPTFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCepD-optimistic_66.csv";
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

void EIC_NC_EPD_28_OPTFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCepD-optimistic_28.csv";
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

//Neutral-current measurements, pessimistic scenario

//A) electron-proton
void EIC_NC_EMP_140_PESFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemP-pessimistic_140.csv";
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

void EIC_NC_EMP_63_PESFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemP-pessimistic_63.csv";
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

void EIC_NC_EMP_44_PESFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemP-pessimistic_44.csv";
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

void EIC_NC_EMP_28_PESFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemP-pessimistic_28.csv";
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

//B) positron-proton
void EIC_NC_EPP_140_PESFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCepP-pessimistic_140.csv";
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

void EIC_NC_EPP_63_PESFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCepP-pessimistic_63.csv";
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

void EIC_NC_EPP_44_PESFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCepP-pessimistic_44.csv";
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

void EIC_NC_EPP_28_PESFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCepP-pessimistic_28.csv";
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

//C) electron-deuteron

void EIC_NC_EMD_88_PESFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemD-pessimistic_88.csv";
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

void EIC_NC_EMD_66_PESFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemD-pessimistic_66.csv";
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

void EIC_NC_EMD_28_PESFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCemD-pessimistic_28.csv";
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

//D) positron-deuteron

void EIC_NC_EPD_88_PESFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCepD-pessimistic_88.csv";
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

void EIC_NC_EPD_66_PESFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCepD-pessimistic_66.csv";
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

void EIC_NC_EPD_28_PESFilter::ReadData()
{
  //Opening file
  fstream f1;
  stringstream datafile("");
  datafile << dataPath() << "rawdata/EIC/NCepD-pessimistic_28.csv";
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




































































