/*
Z pT normalised measurements at 8 TeV, 20 fb^(-1)  -- norm = true
OR
Z pT unnormalised measurements at 8 TeV, 20 fb^(-1) -- norm = false
arXiv:1512.02192
http://hepdata.cedar.ac.uk/view/ins1408516
 * ATLAS Z pT normalised measurement at 8 TeV 20 fb^{-1}.
 * Born level, muon & electron combined measurements
 * Info contained in (un)normalized/output/ZcombPt_born_mXXYY_yZZHH/tab.dat
 * Table 17/29 - 66 GeV <  M_{ll} < 116 GeV  - 0.0 < y_{ll} < 0.4  - 20 datapoints
 * Table 18/30 - 66 GeV <  M_{ll} < 116 GeV  - 0.4 < y_{ll} < 0.8  - 20 datapoints
 * Table 19/31 - 66 GeV <  M_{ll} < 116 GeV  - 0.8 < y_{ll} < 1.2  - 20 datapoints
 * Table 20/32 - 66 GeV <  M_{ll} < 116 GeV  - 1.2 < y_{ll} < 1.6  - 20 datapoints
 * Table 21/33 - 66 GeV <  M_{ll} < 116 GeV  - 1.6 < y_{ll} < 2.0  - 20 datapoints
 * Table 22/34 - 66 GeV <  M_{ll} < 116 GeV  - 2.0 < y_{ll} < 2.4  - 20 datapoints
 * Table 23/35 - 12 GeV <  M_{ll} < 20  GeV  - 0.0 < y_{ll} < 2.4  - 8 datapoints
 * Table 24/36 - 20 GeV <  M_{ll} < 30  GeV  - 0.0 < y_{ll} < 2.4  - 8 datapoints
 * Table 25/37 - 30 GeV <  M_{ll} < 46  GeV  - 0.0 < y_{ll} < 2.4  - 8 datapoints
 * Table 26/38 - 46 GeV <  M_{ll} < 66  GeV  - 0.0 < y_{ll} < 2.4  - 20 datapoints
 * Table 27/39 - 66 GeV <  M_{ll} < 116 GeV  - 0.0 < y_{ll} < 2.4  - 43 datapoints (discarded if YDIST data included)
 * Table 28/40 - 116GeV <  M_{ll} < 150 GeV  - 0.0 < y_{ll} < 2.4  - 20 datapoints

- All systematic uncertainties are assumed to be multiplicative.
- Custom uncertainty descriptions are assumed to allow for cross-correlations
among the five differential distributions. 
- There are 99 sources of correlated systematics
*/

#include "ATLASZPT8TEV.h"

/****************************************************/
/* Central M_{ll} region, Z peak, 6 bins in rapidity*/
/* Normalised distributions                         */

void  ATLASZPT8TEVYDISTNORMFilter::ReadData()
{
  // Opening files
  fstream y1,y2,y3,y4,y5,y6;
  
  /* Table 17 - 66 GeV <  M_{ll} < 116 GeV  - 0.0 < y_{ll} < 0.4  - 20 datapoints*/
  stringstream yfile1("");
  yfile1 << dataPath() 
	 << "rawdata/" << fSetName << "/normalized/output/ZcombPt_born_m66116_y0004/tab.dat";
  y1.open(yfile1.str().c_str(), ios::in);
  if (y1.fail()) 
    {
      cerr << "Error opening data file " << yfile1.str() << endl;
      exit(-1);
    }
  
  /* Table 18 - 66 GeV <  M_{ll} < 116 GeV  - 0.4 < y_{ll} < 0.8  - 20 datapoints*/
  stringstream yfile2("");
  yfile2 << dataPath() 
	 << "rawdata/" << fSetName << "/normalized/output/ZcombPt_born_m66116_y0408/tab.dat";
  y2.open(yfile2.str().c_str(), ios::in);
  if (y2.fail()) 
    {
      cerr << "Error opening data file " << yfile2.str() << endl;
      exit(-1);
    }
  
  /* Table 19 - 66 GeV <  M_{ll} < 116 GeV  - 0.8 < y_{ll} < 1.2  - 20 datapoints*/
  stringstream yfile3("");
  yfile3 << dataPath() 
	 << "rawdata/" << fSetName << "/normalized/output/ZcombPt_born_m66116_y0812/tab.dat";
  y3.open(yfile3.str().c_str(), ios::in);
  if (y3.fail()) 
    {
      cerr << "Error opening data file " << yfile3.str() << endl;
      exit(-1);
    }
  
  /* Table 20 - 66 GeV <  M_{ll} < 116 GeV  - 1.2 < y_{ll} < 1.6  - 20 datapoints*/
  stringstream yfile4("");
  yfile4 << dataPath() 
	 << "rawdata/" << fSetName << "/normalized/output/ZcombPt_born_m66116_y1216/tab.dat";
  y4.open(yfile4.str().c_str(), ios::in);
  if (y4.fail()) 
    {
      cerr << "Error opening data file " << yfile4.str() << endl;
      exit(-1);
    }
  
  /* Table 21 - 66 GeV <  M_{ll} < 116 GeV  - 1.6 < y_{ll} < 2.0  - 20 datapoints*/
  stringstream yfile5("");
  yfile5 << dataPath() 
	 << "rawdata/" << fSetName << "/normalized/output/ZcombPt_born_m66116_y1620/tab.dat";
  y5.open(yfile5.str().c_str(), ios::in);
  if (y5.fail()) 
    {
      cerr << "Error opening data file " << yfile5.str() << endl;
      exit(-1);
    }
  
  /* Table 22 - 66 GeV <  M_{ll} < 116 GeV  - 2.0 < y_{ll} < 2.4  - 20 datapoints*/
  stringstream yfile6("");
  yfile6 << dataPath() 
	 << "rawdata/" << fSetName << "/normalized/output/ZcombPt_born_m66116_y2024/tab.dat";
  y6.open(yfile6.str().c_str(), ios::in);
  if (y6.fail()) 
    {
      cerr << "Error opening data file " << yfile6.str() << endl;
      exit(-1);
    }
  
  //Starting filter
  // Data per distribution
  int ydata = 20;
  
  //Skip thre description lines in each file
  string head,line,linec;
  for(int i=0; i<3; i++)
    {
      getline(y1,head);
      getline(y2,head);
      getline(y3,head);
      getline(y4,head);
      getline(y5,head);
      getline(y6,head);
    }
  
  for(int i=0; i<fNData;i++)
    {
      double ptZ,yZ,ddum,fTot;

      if(i<ydata){
	yZ = 0.2;
 	getline(y1,line);
	getline(y1,linec);
      }
      else if(i>ydata-1 && i < 2*ydata){
	yZ = 0.6;
	getline(y2,line);
	getline(y2,linec);
      }
      else if(i>2*ydata-1 && i < 3*ydata){
	yZ = 1.0;
	getline(y3,line);
	getline(y3,linec);
      }
      else if(i>3*ydata-1 && i < 4*ydata){
	yZ = 1.4;
	getline(y4,line);
	getline(y4,linec);
      }
      else if(i>4*ydata-1 && i < 5*ydata){
	yZ = 1.8;
	getline(y5,line);
	getline(y5,linec);
      }
      else if(i>5*ydata-1 && i < 6*ydata){
	yZ = 2.2;
	getline(y6,line);
	getline(y6,linec);
      }
      else{
	cout << "Wrong indices in ATLASZPT8TEV.cc !!! " << endl;
	exit(-1);
      }

      istringstream lstream(line);
      //      const double pb2fb = 1e3; // Must multiply from pb to fb

      lstream >> ptZ >> ddum >> ddum; 
      fKin1[i] = yZ;         // y_Z
      fKin2[i] = ptZ*ptZ;     // pT_Z^2
      fKin3[i] = 8000;        // sqrt(S)
      lstream >> fData[i];     //differential distribution
      lstream >> fStat[i];     //statistical uncertainty
      lstream >> fSys[i][0].add >> fTot;
      
      fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
      fSys[i][0].type = MULT;
      fSys[i][0].name = "UNCORR";
      // Read Correlated systematics (given in percentage value!!!)
      istringstream lstreamcorr(linec);	
      for ( int k = 1; k < fNSys; k++ )
	{
	  lstreamcorr >> fSys[i][k].mult; 
	  fSys[i][k].add  = fSys[i][k].mult*fData[i]/1e2; 
	  fSys[i][k].type  = MULT;
	  fSys[i][k].name = "CORR";
	}
    }
  y1.close();
  y2.close();
  y3.close();
  y4.close();
  y5.close();
  y6.close();
}

/********************************/
/* Un-normailsed distributions  */
/********************************/

void  ATLASZPT8TEVYDISTFilter::ReadData()
{
  
  // Opening files
  fstream y1,y2,y3,y4,y5,y6;
  
  /* Table 29 - 66 GeV <  M_{ll} < 116 GeV  - 0.0 < y_{ll} < 0.4  - 20 datapoints*/
  stringstream yfile1("");
  yfile1 << dataPath() 
	 << "rawdata/" << fSetName << "/unnormalized/output/ZcombPt_born_m66116_y0004/tab.dat";
  y1.open(yfile1.str().c_str(), ios::in);
  if (y1.fail()) 
    {
      cerr << "Error opening data file " << yfile1.str() << endl;
      exit(-1);
    }
  
  /* Table 30 - 66 GeV <  M_{ll} < 116 GeV  - 0.4 < y_{ll} < 0.8  - 20 datapoints*/
  stringstream yfile2("");
  yfile2 << dataPath() 
	 << "rawdata/" << fSetName << "/unnormalized/output/ZcombPt_born_m66116_y0408/tab.dat";
  y2.open(yfile2.str().c_str(), ios::in);
  if (y2.fail()) 
    {
      cerr << "Error opening data file " << yfile2.str() << endl;
      exit(-1);
    }
  
  /* Table 31 - 66 GeV <  M_{ll} < 116 GeV  - 0.8 < y_{ll} < 1.2  - 20 datapoints*/
  stringstream yfile3("");
  yfile3 << dataPath() 
	 << "rawdata/" << fSetName << "/unnormalized/output/ZcombPt_born_m66116_y0812/tab.dat";
  y3.open(yfile3.str().c_str(), ios::in);
  if (y3.fail()) 
    {
      cerr << "Error opening data file " << yfile3.str() << endl;
      exit(-1);
    }
  
  /* Table 32 - 66 GeV <  M_{ll} < 116 GeV  - 1.2 < y_{ll} < 1.6  - 20 datapoints*/
  stringstream yfile4("");
  yfile4 << dataPath() 
	 << "rawdata/" << fSetName << "/unnormalized/output/ZcombPt_born_m66116_y1216/tab.dat";
  y4.open(yfile4.str().c_str(), ios::in);
  if (y4.fail()) 
    {
      cerr << "Error opening data file " << yfile4.str() << endl;
      exit(-1);
    }
  
  /* Table 33 - 66 GeV <  M_{ll} < 116 GeV  - 1.6 < y_{ll} < 2.0  - 20 datapoints*/
  stringstream yfile5("");
  yfile5 << dataPath() 
	 << "rawdata/" << fSetName << "/unnormalized/output/ZcombPt_born_m66116_y1620/tab.dat";
  y5.open(yfile5.str().c_str(), ios::in);
  if (y5.fail()) 
    {
      cerr << "Error opening data file " << yfile5.str() << endl;
      exit(-1);
    }
  
  /* Table 34 - 66 GeV <  M_{ll} < 116 GeV  - 2.0 < y_{ll} < 2.4  - 20 datapoints*/
  stringstream yfile6("");
  yfile6 << dataPath() 
	 << "rawdata/" << fSetName << "/unnormalized/output/ZcombPt_born_m66116_y2024/tab.dat";
  y6.open(yfile6.str().c_str(), ios::in);
  if (y6.fail()) 
    {
      cerr << "Error opening data file " << yfile6.str() << endl;
      exit(-1);
    }
  
  //Starting filter
  // Data per distribution
  int ydata = 20;
  
  //Skip thre description lines in each file
  string head,line,linec;
  for(int i=0; i<3; i++)
    {
      getline(y1,head);
      getline(y2,head);
      getline(y3,head);
      getline(y4,head);
      getline(y5,head);
      getline(y6,head);
    }

  for(int i=0; i<fNData;i++)
    {
      double ptZ,yZ,ddum,fTot;

      if(i<ydata){
	yZ = 0.2;
 	getline(y1,line);
	getline(y1,linec);
      }
      else if(i>ydata-1 && i < 2*ydata){
	yZ = 0.6;
	getline(y2,line);
	getline(y2,linec);
      }
      else if(i>2*ydata-1 && i < 3*ydata){
	yZ = 1.0;
	getline(y3,line);
	getline(y3,linec);
      }
      else if(i>3*ydata-1 && i < 4*ydata){
	yZ = 1.4;
	getline(y4,line);
	getline(y4,linec);
      }
      else if(i>4*ydata-1 && i < 5*ydata){
	yZ = 1.8;
	getline(y5,line);
	getline(y5,linec);
      }
      else if(i>5*ydata-1 && i < 6*ydata){
	yZ = 2.2;
	getline(y6,line);
	getline(y6,linec);
      }
      else{
	cout << "Wrong indices in ATLASZPT8TEV.cc !!! " << endl;
	exit(-1);
      }
      istringstream lstream(line);
      const double pb2fb = 1e3; // Must multiply from pb to fb

      lstream >> ptZ >> ddum >> ddum; 
      fKin1[i] = yZ;         // y_Z
      fKin2[i] = ptZ*ptZ;     // pT_Z^2
      fKin3[i] = 8000;        // sqrt(S)
      lstream >> fData[i];     //differential distribution
      lstream >> fStat[i];     //statistical uncertainty
      lstream >> fSys[i][0].add >> fTot;

      fData[i] *= pb2fb;
      fStat[i] *= pb2fb;
      fSys[i][0].add *= pb2fb;
      fTot *= pb2fb;

      fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
      fSys[i][0].type = MULT;
      fSys[i][0].name = "UNCORR";
      // Read Correlated systematics (given in percentage value!!!)
      istringstream lstreamcorr(linec);	
      for ( int k = 1; k < fNSys-1; k++ )
	{
	  lstreamcorr >> fSys[i][k].mult; //what if the errors in these files are given as relative errors and not percentage relative?
	  //	  fSys[i][k].mult *= 1e2 ; ADDED - was that the bug?? --- No, I do not think so, they are given in %, not as relative
	  //	  fSys[i][k].mult *= pb2fb; ADDED - was that the bug?? --- No, I do not think so, they are given in %, not as absolute
	  fSys[i][k].add  = fSys[i][k].mult*fData[i]/1e2; 
	  fSys[i][k].type  = MULT;
	  fSys[i][k].name = "CORR";
	}
      // Luminosity Uncertainty: 2.8%
      fSys[i][fNSys-1].mult = 2.8;
      fSys[i][fNSys-1].add  = fSys[i][fNSys-1].mult*fData[i]/1e2;
      fSys[i][fNSys-1].type = MULT;
      fSys[i][fNSys-1].name = "ATLASLUMI12";
    }
  y1.close();
  y2.close();
  y3.close();
  y4.close();
  y5.close();
  y6.close();
}

/****************************************************************************************/
/* Low mass region, inclusive in rapidity, 4 bins in lepton invariant mass below Z peak */
/* Higs mass region, inclusive in rapidity, 1 bin in lepton invariant mass above Z peak */
/****************************************************************************************/

void  ATLASZPT8TEVMDISTFilter::ReadData()
{
  // Opening files
  fstream m1,m2,m3,m4,m5;

  bool norm = false; //Set to true for normalised data, false for unnormalised

  if(norm)
    {  
      /* Table 23 - 12 GeV <  M_{ll} < 20 GeV  - 0.0 < y_{ll} < 2.4  - 8 datapoints*/  
      stringstream mfile1("");
      mfile1 << dataPath() 
	     << "rawdata/" << fSetName << "/normalized/output/ZcombPt_born_m1220_y0024/tab.dat";
      m1.open(mfile1.str().c_str(), ios::in);
      if (m1.fail()) 
	{
	  cerr << "Error opening data file " << mfile1.str() << endl;
	  exit(-1);
	}
      
      /* Table 24 - 20 GeV <  M_{ll} < 30 GeV  - 0.0 < y_{ll} < 2.4  - 8 datapoints*/  
      stringstream mfile2("");
      mfile2 << dataPath() 
	     << "rawdata/" << fSetName << "/normalized/output/ZcombPt_born_m2030_y0024/tab.dat";
      m2.open(mfile2.str().c_str(), ios::in);
      if (m2.fail()) 
	{
	  cerr << "Error opening data file " << mfile2.str() << endl;
	  exit(-1);
	}
      
      /* Table 25 - 30 GeV <  M_{ll} < 46 GeV  - 0.0 < y_{ll} < 2.4  - 8 datapoints*/  
      stringstream mfile3("");
      mfile3 << dataPath() 
	     << "rawdata/" << fSetName << "/normalized/output/ZcombPt_born_m3046_y0024/tab.dat";
      m3.open(mfile3.str().c_str(), ios::in);
      if (m3.fail()) 
	{
	  cerr << "Error opening data file " << mfile3.str() << endl;
	  exit(-1);
	}
      
      /* Table 26 - 46 GeV <  M_{ll} < 66 GeV  - 0.0 < y_{ll} < 2.4  - 20 datapoints*/  
      stringstream mfile4("");
      mfile4 << dataPath() 
	     << "rawdata/" << fSetName << "/normalized/output/ZcombPt_born_m4666_y0024/tab.dat";
      m4.open(mfile4.str().c_str(), ios::in);
      if (m4.fail()) 
	{
	  cerr << "Error opening data file " << mfile4.str() << endl;
	  exit(-1);
	}
      
      /* Table 28 - 116 GeV <  M_{ll} < 150 GeV  - 0.0 < y_{ll} < 2.4  - 20 datapoints*/  
      stringstream mfile5("");
      mfile5 << dataPath() 
	     << "rawdata/" << fSetName << "/normalized/output/ZcombPt_born_m116150_y0024/tab.dat";
      m5.open(mfile5.str().c_str(), ios::in);
      if (m5.fail()) 
	{
	  cerr << "Error opening data file " << mfile5.str() << endl;
	  exit(-1);
	}
    }
  else
    {
      /* Table 35 - 12 GeV <  M_{ll} < 20 GeV  - 0.0 < y_{ll} < 2.4  - 8 datapoints*/  
      stringstream mfile1("");
      mfile1 << dataPath() 
	     << "rawdata/" << fSetName << "/unnormalized/output/ZcombPt_born_m1220_y0024/tab.dat";
      m1.open(mfile1.str().c_str(), ios::in);
      if (m1.fail()) 
	{
	  cerr << "Error opening data file " << mfile1.str() << endl;
	  exit(-1);
	}
      
      /* Table 36 - 20 GeV <  M_{ll} < 30 GeV  - 0.0 < y_{ll} < 2.4  - 8 datapoints*/  
      stringstream mfile2("");
      mfile2 << dataPath() 
	     << "rawdata/" << fSetName << "/unnormalized/output/ZcombPt_born_m2030_y0024/tab.dat";
      m2.open(mfile2.str().c_str(), ios::in);
      if (m2.fail()) 
	{
	  cerr << "Error opening data file " << mfile2.str() << endl;
	  exit(-1);
	}
      
      /* Table 37 - 30 GeV <  M_{ll} < 46 GeV  - 0.0 < y_{ll} < 2.4  - 8 datapoints*/  
      stringstream mfile3("");
      mfile3 << dataPath() 
	     << "rawdata/" << fSetName << "/unnormalized/output/ZcombPt_born_m3046_y0024/tab.dat";
      m3.open(mfile3.str().c_str(), ios::in);
      if (m3.fail()) 
	{
	  cerr << "Error opening data file " << mfile3.str() << endl;
	  exit(-1);
	}
      
      /* Table 38 - 46 GeV <  M_{ll} < 66 GeV  - 0.0 < y_{ll} < 2.4  - 20 datapoints*/  
      stringstream mfile4("");
      mfile4 << dataPath() 
	     << "rawdata/" << fSetName << "/unnormalized/output/ZcombPt_born_m4666_y0024/tab.dat";
      m4.open(mfile4.str().c_str(), ios::in);
      if (m4.fail()) 
	{
	  cerr << "Error opening data file " << mfile4.str() << endl;
	  exit(-1);
	}
      
      /* Table 40 - 116 GeV <  M_{ll} < 150 GeV  - 0.0 < y_{ll} < 2.4  - 20 datapoints*/  
      stringstream mfile5("");
      mfile5 << dataPath() 
	     << "rawdata/" << fSetName << "/unnormalized/output/ZcombPt_born_m116150_y0024/tab.dat";
      m5.open(mfile5.str().c_str(), ios::in);
      if (m5.fail()) 
	{
	  cerr << "Error opening data file " << mfile5.str() << endl;
	  exit(-1);
	}
    }
  
  //Starting filter
  int nbin1 = 8;
  int nbin2 = 8;
  int nbin3 = 8;
  int nbin4 = 20;
  //  int nbin5 = 43;
  int nbin6 = 20;

  //Skip thre description lines in each file
  string head,line,linec;
  for(int i=0; i<3; i++)
    {
      getline(m1,head);
      getline(m2,head);
      getline(m3,head);
      getline(m4,head);
      getline(m5,head);
    }
  for(int i=0; i<fNData;i++)
    {
      double ptZ,mLL,ddum,fTot;
      
      if(i<nbin1){
	mLL = 16.;
	getline(m1,line);
 	getline(m1,linec);
      }
      else if(i>nbin1-1 && i < nbin1 + nbin2){
	mLL = 25.;
	getline(m2,line);
	getline(m2,linec);
      }
      else if(i> nbin1 + nbin2-1 && i < nbin1 + nbin2 + nbin3){
	mLL = 38.;
	getline(m3,line);
	getline(m3,linec);
      }
      else if(i>nbin1 + nbin2 + nbin3 - 1 && i < nbin1+nbin2+nbin3+nbin4){
	mLL = 56.;
	getline(m4,line);
	getline(m4,linec);
      }
      else if(i>nbin1+nbin2+nbin3+nbin4-1 && i < nbin1+nbin2+nbin3+nbin4+nbin6){
	mLL = 138.;
	getline(m5,line);
	getline(m5,linec);
      }
      else{
	cout << "Wrong indices in ATLASZPT8TEV.cc !!! " << endl;
	exit(-1);
      }

      const double pb2fb = 1000.; // Must multiply from pb to fb
      istringstream lstream(line);
      lstream >> ptZ >> ddum >> ddum; 
      fKin1[i] = ptZ;         // P_T^(Z)
      fKin2[i] = mLL*mLL;     // M_Z^2
      fKin3[i] = 8000;        // sqrt(S)
      lstream >> fData[i];    //differential distribution
      lstream >> fStat[i];    //statistical uncertainty
      lstream >> fSys[i][0].add >> fTot;
      fSys[i][0].mult = fSys[i][0].add/fData[i]*1e2;
      fSys[i][0].type = MULT;
      fSys[i][0].name = "UNCORR";
      
      if(norm)
	{ 
	  cout << "No need to normalise normalised xsec!" << endl;
	}
      else{
	fData[i] *= pb2fb;
	fStat[i] *= pb2fb;
	fSys[i][0].add *= pb2fb;
	fTot *= pb2fb;
      }
      
      // Read Correlated systematics (given in percentage in the data file - wrong README in ATLAS!)
      istringstream lstreamcorr(linec);	
      for ( int k = 1; k < fNSys-1; k++ )
	{
	  lstreamcorr >> fSys[i][k].mult;
	  fSys[i][k].add  = fSys[i][k].mult*fData[i]/1e2;
	  fSys[i][k].type  = MULT;
	  fSys[i][k].name = "CORR";
	}      
      // Luminosity Uncertainty: 2.8%
      fSys[i][fNSys-1].mult = 2.8;
      fSys[i][fNSys-1].add  = fData[i]*fSys[i][fNSys-1].mult/1e2;
      fSys[i][fNSys-1].type = MULT;
      fSys[i][fNSys-1].name = "ATLASLUMI12";
    }  
  m1.close();
  m2.close();
  m3.close();
  m4.close();
  m5.close();
}


