/*
Insert here some description about the experiment



*/

#include "ATLAS_TTB_DIFF_8TEV_LJ.h"
#include <iomanip>

//Function that construct and decomposes the global covariance matrix of statistical correlations
//this function MUST have an input argument (NORM, ABS) to selct the relevant tables)
//this functions MUST have an output, which is the array of additional fake systematic uncertainties


void statcorrs(int dist, string norm, double** extrasys) 
{
  int ndist=4;
  int ndat[] ={8, 5, 5, 7};
  int totdata=0.;
  for(int i=0; i<ndist; i++)
    {
      totdata += ndat[i]; 
    }

  //Tabcorr numbering for covariance matrix blocks
  string** tabcorr = new string*[ndist];
  for(int i=0; i<ndist; i++)
    {
      tabcorr[i] = new string[ndist];
    }
  
  string* tabstat = new string[ndist];
    
    if(norm=="ABS")
    {
      tabcorr[0][0]="172";
      tabcorr[0][1]="168";
      tabcorr[0][2]="173";
      tabcorr[0][3]="174";
      
      tabcorr[1][1]="167";
      tabcorr[1][2]="169";
      tabcorr[1][3]="170";   

      tabcorr[2][2]="176";
      tabcorr[2][3]="177";

      tabcorr[3][3]="179";

      tabstat[0]="29";
      tabstat[1]="31";
      tabstat[2]="27";
      tabstat[3]="23";
    }

  else if(norm=="NORM")
    {
      tabcorr[0][0]="227";
      tabcorr[0][1]="223";
      tabcorr[0][2]="228";
      tabcorr[0][3]="229";
      
      tabcorr[1][1]="222";
      tabcorr[1][2]="224";
      tabcorr[1][3]="225";   

      tabcorr[2][2]="231";
      tabcorr[2][3]="232";

      tabcorr[3][3]="234";
      
      tabstat[0]="30";
      tabstat[1]="32";
      tabstat[2]="28";
      tabstat[3]="24";
    }
  
  else
    {
      cout << "Normalisation not defined" << endl;
      exit(-1);
    }

  //Read statistical uncertainty
  double* stat = new double[totdata];
  int p=0;

  for(int idist=0; idist<ndist; idist++)
    { 

      fstream f2, f3;
      stringstream datafile_stat("");
      datafile_stat << dataPath()
		    <<"rawdata/ATLAS_TTB_DIFF_8TEV_LJ/HEPData-ins1404878-v2-csv/Table" << tabstat[idist] << ".csv";
      f3.open(datafile_stat.str().c_str(), ios::in);
 
      if (f3.fail()) 
	{
	  cerr << "Error opening data file " << datafile_stat.str() << endl;
	  exit(-1);
	}
      
      string qline;
      
      for(int i=0; i<11; i++)
	{
	  getline(f3,qline);
	}
    
      double* data = new double[totdata];
	  
      for(int j=0; j<ndat[idist]; j++)
	{
	  char comma;
	  double ddum;
	  getline(f3,qline);
	  istringstream kstream(qline);

	  kstream >> ddum      >> comma
		  >> ddum      >> comma
		  >> ddum      >> comma
		  >> data[j+p] >> comma
		  >> stat[j+p] >> comma;

	  stat[j+p] = stat[j+p]/100.*data[j+p];
	  
	}
      
      p += ndat[idist];
      f3.close();

    }

  //Read correlation coefficients 
  double** corrcoeff = new double*[totdata];
  int n=0;

  for(int i=0; i<totdata; i++)
    {
      corrcoeff[i] = new double[totdata];
    }

  for(int i=0; i<totdata; i++)
    {
      for(int j=0; j<totdata; j++)
	{
	  corrcoeff[i][j]=0.;
	}
    }

  for(int idist=0; idist<ndist; idist++)
    {
      int m;

      if(idist==0)
	m=0;
      else if(idist==1)
	m=ndat[0];
      else if(idist==2)
	m=ndat[0]+ndat[1];
      else if(idist==3)
	m=ndat[0]+ndat[1]+ndat[2];
      else
	{
	  cout << "Error with distribution assignment" << endl;
	  exit(-1);
	}

      for(int jdist=idist; jdist<ndist; jdist++)
	{ 
	 
	  fstream f2;
	  stringstream datafile_corr("");
	  datafile_corr << dataPath()
			<<"rawdata/ATLAS_TTB_DIFF_8TEV_LJ/HEPData-ins1404878-v2-csv/Table" << tabcorr[idist][jdist] << ".csv";

	  f2.open(datafile_corr.str().c_str(), ios::in);
	  
	  if (f2.fail()) 
	    {
	      cerr << "Error opening data file " << datafile_corr.str() << endl;
	      exit(-1);
	    }
	  
	  string pline;
	  
	  for(int i=0; i<10; i++)
	    {
	      getline(f2,pline);
	    }

	  for(int i=0; i<ndat[idist]; i++)
	    {
	      
	      for(int j=0; j<ndat[jdist]; j++)
		{
		  getline(f2,pline);
		  istringstream lstream(pline);

		  int idum;
		  char comma;

		  lstream >> idum >> comma
			  >> idum >> comma
			  >> corrcoeff[i+n][j+m];
		  
		  corrcoeff[i+n][j+m] *= stat[i+n]*stat[j+m];
		  
		}

	    }

	  m += ndat[jdist];
	  f2.close();

	}

      n += ndat[idist];  

    }

  //Symmetrise correlation matrix
  for(int i=0; i<totdata; i++)
    {
      for(int j=i; j<totdata; j++)
	{
	  corrcoeff[j][i] = corrcoeff[i][j];
	}
    }

  //Generate artificial systematics
  double** syscor = new double*[totdata];
  for(int i=0; i<totdata; i++)
    syscor[i] = new double[totdata];
  
  if(!genArtSys(totdata,corrcoeff,syscor))
    {
      cerr << "An error was encountered in the generation of the artificial systematics " << endl;
      exit(-1);
    }

  for(int k=0; k<ndat[dist]; k++)
    { 
      int a;
      if(dist==0) a=0;
      else if (dist==1) a=ndat[0];
      else if (dist==2) a=ndat[0]+ndat[1];
      else if (dist==3) a=ndat[0]+ndat[1]+ndat[2];
      else exit(-1);
      
      for(int j=0; j<totdata; j++)
	{
	  extrasys[k][j]=syscor[a+k][j];
	}    
    }

  delete[] syscor;
  delete[] corrcoeff;
  delete[] tabcorr;
  
  return;
  
}

//A - Absolute distributions

//1) Distribution differential in top quark transverse momentum

void ATLAS_TTB_DIFF_8TEV_LJ_TPTFilter::ReadData()
{
  //Opening files
  fstream f1;
  
  //Central values, total statistic and systematic uncertainties
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/ATLAS_TTB_DIFF_8TEV_LJ/HEPData-ins1404878-v2-csv/Table29.csv";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  string line;
  istringstream lstream(line);
  
  for(int i=0; i<22; i++)
    {
      getline(f1,line);
    }
  
  int realsys=56;
  
  double** extrasys = new double*[25];
  for(int i=0; i<25; i++)
    extrasys[i] = new double[25];

  for(int i=0; i<fNData; i++)
    {
      char comma;
      double ddum;
      getline(f1,line);
      istringstream lstream(line);
      
      lstream >> fKin1[i] >> comma
	      >> ddum     >> comma
	      >> ddum     >> comma
	      >> fData[i] >> comma
	      >> ddum >> comma >> comma
	      >> ddum     >> comma >> comma;
      
      fStat[i] = 0.;
      fKin2[i] = Mt*Mt;       
      fKin3[i] = 8000;         //sqrt(s)
      fStat[i] = fStat[i]/100.*fData[i];
      
      //Systematic uncertainties
      for(int j=0; j<realsys; j++)
	{
	  double sysR, sysL;
	  
	  if(j<realsys-1)
	    lstream >> sysR >> comma >> comma >> sysL >> comma >> comma;
	  else
	    lstream >> sysR >> comma >> comma >> sysL >> comma;
	  
	  double tmp1 = sysR;
	  double tmp2 = sysL;
	  
	  //Case 1: sysR and sysL are both negative
	  if(tmp1<0.0 && tmp2<0.0)
	    {
	      if(tmp2<tmp1)
		{
		  sysR = 0.0;
		  sysL = tmp2;
		}
	      if(tmp2>tmp1)
		{
		  sysR = 0.0;
		  sysL = tmp1;
		}
	    }

	  //Case 2: sysR and sysL are both positive
	  if(tmp1>0.0 && tmp2>0.0)
	    {
	      if(tmp1>tmp2)
		{
		  sysR = tmp1;
		  sysL = 0.0;
		}
	      if(tmp1<tmp2)
		{
		  sysR = tmp2;
		  sysL = 0.0;
		}
	    }

	  //Case3: sysR is negative and sysL is positive
	  if(tmp1<0.0 && tmp2>0.0)
	    {
	      sysR = tmp2;
	      sysL = tmp1;
	    }

	  sysR = sysR/sqrt(2.);
	  sysL = sysL/sqrt(2.);
	  
	  if(j<realsys-1)

	    {
	      fSys[i][2*j].mult = sysR;
	      fSys[i][2*j].add  = fSys[i][2*j].mult*fData[i]/100;
	      fSys[i][2*j].type = MULT;
	      fSys[i][2*j].name = "SYSCOR";
	      
	      fSys[i][2*j+1].mult = sysR;
	      fSys[i][2*j+1].add  = fSys[i][2*j+1].mult*fData[i]/100;
	      fSys[i][2*j+1].type = MULT;
	      fSys[i][2*j+1].name = "SYSCOR";  
	    }
	  else
	    {
	      fSys[i][2*j].mult = sysR;
	      fSys[i][2*j].add  = fSys[i][2*j].mult*fData[i]/100;
	      fSys[i][2*j].type = MULT;
	      fSys[i][2*j].name = "ATLASLUMI12";
	    }
	    
	}
      
      //Define artificial systematics
      for(int j=2*realsys-1; j<fNSys; j++)
	{

	  statcorrs(0,"ABS",extrasys);
	  fSys[i][j].add = extrasys[i][j-2*realsys+1];
	  fSys[i][j].mult = fSys[i][j].add/fData[i]*100.;
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = "STATCOR";
	  
	}
    
    }
 
  f1.close();

}

//2) Distribution differential in top rapidity

void ATLAS_TTB_DIFF_8TEV_LJ_TRAPFilter::ReadData()
{
  //Opening files
  fstream f1;
  
  //Central values, total statistic and systematic uncertainties
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/ATLAS_TTB_DIFF_8TEV_LJ/HEPData-ins1404878-v2-csv/Table31.csv";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  string line;
  istringstream lstream(line);
  
  for(int i=0; i<19; i++)
    {
      getline(f1,line);
    }
  
  int realsys=56;
  
  double** extrasys = new double*[25];
  for(int i=0; i<25; i++)
    extrasys[i] = new double[25];

  for(int i=0; i<fNData; i++)
    {
      char comma;
      double ddum;
      getline(f1,line);
      istringstream lstream(line);
      
      lstream >> fKin1[i] >> comma
	      >> ddum     >> comma
	      >> ddum     >> comma
	      >> fData[i] >> comma
	      >> ddum     >> comma >> comma
	      >> ddum     >> comma >> comma;
      
      fStat[i] = 0.;
      fKin2[i] = Mt*Mt;       
      fKin3[i] = 8000;         //sqrt(s)
      fStat[i] = fStat[i]/100.*fData[i];
      
      //Systematic uncertainties
      for(int j=0; j<realsys; j++)
	{
	  double sysR, sysL;
	  
	  if(j<realsys-1)
	    lstream >> sysR >> comma >> comma >> sysL >> comma >> comma;
	  else
	    lstream >> sysR >> comma >> comma >> sysL >> comma;
	  
	  double tmp1 = sysR;
	  double tmp2 = sysL;
	  
	  //Case 1: sysR and sysL are both negative
	  if(tmp1<0.0 && tmp2<0.0)
	    {
	      if(tmp2<tmp1)
		{
		  sysR = 0.0;
		  sysL = tmp2;
		}
	      if(tmp2>tmp1)
		{
		  sysR = 0.0;
		  sysL = tmp1;
		}
	    }

	  //Case 2: sysR and sysL are both positive
	  if(tmp1>0.0 && tmp2>0.0)
	    {
	      if(tmp1>tmp2)
		{
		  sysR = tmp1;
		  sysL = 0.0;
		}
	      if(tmp1<tmp2)
		{
		  sysR = tmp2;
		  sysL = 0.0;
		}
	    }

	  //Case3: sysR is negative and sysL is positive
	  if(tmp1<0.0 && tmp2>0.0)
	    {
	      sysR = tmp2;
	      sysL = tmp1;
	    }

	  sysR = sysR/sqrt(2.);
	  sysL = sysL/sqrt(2.);
	  
	  if(j<realsys-1)

	    {
	      fSys[i][2*j].mult = sysR;
	      fSys[i][2*j].add  = fSys[i][2*j].mult*fData[i]/100;
	      fSys[i][2*j].type = MULT;
	      fSys[i][2*j].name = "SYSCOR";
	      
	      fSys[i][2*j+1].mult = sysR;
	      fSys[i][2*j+1].add  = fSys[i][2*j+1].mult*fData[i]/100;
	      fSys[i][2*j+1].type = MULT;
	      fSys[i][2*j+1].name = "SYSCOR";  
	    }
	  else
	    {
	      fSys[i][2*j].mult = sysR;
	      fSys[i][2*j].add  = fSys[i][2*j].mult*fData[i]/100;
	      fSys[i][2*j].type = MULT;
	      fSys[i][2*j].name = "ATLASLUMI12";
	    }
	    
	}
      
      //Define artificial systematics
      for(int j=2*realsys-1; j<fNSys; j++)
	{

	  statcorrs(1,"ABS",extrasys);
	  fSys[i][j].add = extrasys[i][j-2*realsys+1];
	  fSys[i][j].mult = fSys[i][j].add/fData[i]*100.;
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = "STATCOR";
	  
	}
    
    }
 
  f1.close();

}

//3) Distribution differential in top pair rapidity

void ATLAS_TTB_DIFF_8TEV_LJ_TTRAPFilter::ReadData()
{
  //Opening files
  fstream f1;
  
  //Central values, total statistic and systematic uncertainties
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/ATLAS_TTB_DIFF_8TEV_LJ/HEPData-ins1404878-v2-csv/Table27.csv";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  string line;
  istringstream lstream(line);
  
  for(int i=0; i<19; i++)
    {
      getline(f1,line);
    }
  
  int realsys=56;
  
  double** extrasys = new double*[25];
  for(int i=0; i<25; i++)
    extrasys[i] = new double[25];

  for(int i=0; i<fNData; i++)
    {
      char comma;
      double ddum;
      getline(f1,line);
      istringstream lstream(line);
      
      lstream >> fKin1[i] >> comma
	      >> ddum     >> comma
	      >> ddum     >> comma
	      >> fData[i] >> comma
	      >> ddum     >> comma >> comma
	      >> ddum     >> comma >> comma;
      
      fStat[i] = 0.;
      fKin2[i] = Mt*Mt;       
      fKin3[i] = 8000;         //sqrt(s)
      fStat[i] = fStat[i]/100.*fData[i];
      
      //Systematic uncertainties
      for(int j=0; j<realsys; j++)
	{
	  double sysR, sysL;
	  
	  if(j<realsys-1)
	    lstream >> sysR >> comma >> comma >> sysL >> comma >> comma;
	  else
	    lstream >> sysR >> comma >> comma >> sysL >> comma;
	  
	  double tmp1 = sysR;
	  double tmp2 = sysL;
	  
	  //Case 1: sysR and sysL are both negative
	  if(tmp1<0.0 && tmp2<0.0)
	    {
	      if(tmp2<tmp1)
		{
		  sysR = 0.0;
		  sysL = tmp2;
		}
	      if(tmp2>tmp1)
		{
		  sysR = 0.0;
		  sysL = tmp1;
		}
	    }

	  //Case 2: sysR and sysL are both positive
	  if(tmp1>0.0 && tmp2>0.0)
	    {
	      if(tmp1>tmp2)
		{
		  sysR = tmp1;
		  sysL = 0.0;
		}
	      if(tmp1<tmp2)
		{
		  sysR = tmp2;
		  sysL = 0.0;
		}
	    }

	  //Case3: sysR is negative and sysL is positive
	  if(tmp1<0.0 && tmp2>0.0)
	    {
	      sysR = tmp2;
	      sysL = tmp1;
	    }

	  sysR = sysR/sqrt(2.);
	  sysL = sysL/sqrt(2.);
	  
	  if(j<realsys-1)

	    {
	      fSys[i][2*j].mult = sysR;
	      fSys[i][2*j].add  = fSys[i][2*j].mult*fData[i]/100;
	      fSys[i][2*j].type = MULT;
	      fSys[i][2*j].name = "SYSCOR";
	      
	      fSys[i][2*j+1].mult = sysR;
	      fSys[i][2*j+1].add  = fSys[i][2*j+1].mult*fData[i]/100;
	      fSys[i][2*j+1].type = MULT;
	      fSys[i][2*j+1].name = "SYSCOR";  
	    }
	  else
	    {
	      fSys[i][2*j].mult = sysR;
	      fSys[i][2*j].add  = fSys[i][2*j].mult*fData[i]/100;
	      fSys[i][2*j].type = MULT;
	      fSys[i][2*j].name = "ATLASLUMI12";
	    }
	    
	}
      
      //Define artificial systematics
      for(int j=2*realsys-1; j<fNSys; j++)
	{

	  statcorrs(2,"ABS",extrasys);
	  fSys[i][j].add = extrasys[i][j-2*realsys+1];
	  fSys[i][j].mult = fSys[i][j].add/fData[i]*100.;
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = "STATCOR";
	  
	}
    
    }
 
  f1.close();

}

//4) Distribution differential in top pair mass

void ATLAS_TTB_DIFF_8TEV_LJ_TTMFilter::ReadData()
{
  //Opening files
  fstream f1;
  
  //Central values, total statistic and systematic uncertainties
  stringstream datafile("");
  datafile << dataPath()
	   << "rawdata/ATLAS_TTB_DIFF_8TEV_LJ/HEPData-ins1404878-v2-csv/Table23.csv";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) 
    {
      cerr << "Error opening data file " << datafile.str() << endl;
      exit(-1);
    }
  
  string line;
  istringstream lstream(line);
  
  for(int i=0; i<21; i++)
    {
      getline(f1,line);
    }
  
  int realsys=56;
  
  double** extrasys = new double*[25];
  for(int i=0; i<25; i++)
    extrasys[i] = new double[25];

  for(int i=0; i<fNData; i++)
    {
      char comma;
      double ddum;
      getline(f1,line);
      istringstream lstream(line);
      
      lstream >> fKin1[i] >> comma
	      >> ddum     >> comma
	      >> ddum     >> comma
	      >> fData[i] >> comma
	      >> ddum     >> comma >> comma
	      >> ddum     >> comma >> comma;
      
      fStat[i] = 0.;
      fKin2[i] = Mt*Mt;       
      fKin3[i] = 8000;         //sqrt(s)
      
      //Systematic uncertainties
      for(int j=0; j<realsys; j++)
	{
	  double sysR, sysL;
	  
	  if(j<realsys-1)
	    lstream >> sysR >> comma >> comma >> sysL >> comma >> comma;
	  else
	    lstream >> sysR >> comma >> comma >> sysL >> comma;
	  
	  double tmp1 = sysR;
	  double tmp2 = sysL;
	  
	  //Case 1: sysR and sysL are both negative
	  if(tmp1<0.0 && tmp2<0.0)
	    {
	      if(tmp2<tmp1)
		{
		  sysR = 0.0;
		  sysL = tmp2;
		}
	      if(tmp2>tmp1)
		{
		  sysR = 0.0;
		  sysL = tmp1;
		}
	    }

	  //Case 2: sysR and sysL are both positive
	  if(tmp1>0.0 && tmp2>0.0)
	    {
	      if(tmp1>tmp2)
		{
		  sysR = tmp1;
		  sysL = 0.0;
		}
	      if(tmp1<tmp2)
		{
		  sysR = tmp2;
		  sysL = 0.0;
		}
	    }

	  //Case3: sysR is negative and sysL is positive
	  if(tmp1<0.0 && tmp2>0.0)
	    {
	      sysR = tmp2;
	      sysL = tmp1;
	    }

	  sysR = sysR/sqrt(2.);
	  sysL = sysL/sqrt(2.);
	  
	  if(j<realsys-1)

	    {
	      fSys[i][2*j].mult = sysR;
	      fSys[i][2*j].add  = fSys[i][2*j].mult*fData[i]/100;
	      fSys[i][2*j].type = MULT;
	      fSys[i][2*j].name = "SYSCOR";
	      
	      fSys[i][2*j+1].mult = sysR;
	      fSys[i][2*j+1].add  = fSys[i][2*j+1].mult*fData[i]/100;
	      fSys[i][2*j+1].type = MULT;
	      fSys[i][2*j+1].name = "SYSCOR";  
	    }
	  else
	    {
	      fSys[i][2*j].mult = sysR;
	      fSys[i][2*j].add  = fSys[i][2*j].mult*fData[i]/100;
	      fSys[i][2*j].type = MULT;
	      fSys[i][2*j].name = "ATLASLUMI12";
	    }
	    
	}
      
      //Define artificial systematics
      for(int j=2*realsys-1; j<fNSys; j++)
	{

	  statcorrs(3,"ABS",extrasys);
	  fSys[i][j].add = extrasys[i][j-2*realsys+1];
	  fSys[i][j].mult = fSys[i][j].add/fData[i]*100.;
	  fSys[i][j].type = ADD;
	  fSys[i][j].name = "STATCOR";
	  
	}
    
    }
 
  f1.close();

}

