/*********** NEW DATA IN NNPDF4.0 ************

* Points implemented are those reported at https://www.hepdata.net/record/ins1487277, tables 1 to 6.
*
* Full breakdown of the sysstematics from xfitter files. It contains 26 syst to be simmetrized, 
* given as 52 column entries, plus the lumi sys and 3 columns realted to stats.
*
*

bin 1:   0 < |y| < 0.5
========================
points    46 
real sys  26 + lumi 
art sys   37
systot    64

bin 2:   0.5 < |y| < 1.0
========================
points    46 
real sys  26 + lumi 
art sys   37
systot    64

bin 3:   1.5 < |y| < 2.0
========================
points    45 
real sys  26 + lumi 
art sys   36
systot    63

bin 4:   2.0 < |y| < 2.5
========================
points    41 
real sys  26 + lumi 
art sys   32
systot    59

bin 5:   2.5 < |y| < 3.0
========================
points    34 
real sys  26 + lumi 
art sys   25
systot    52

bin 6:   3.0 < |y| < 3.5
========================
points    27 
real sys  26 + lumi 
art sys   18
systot    45


tot points         239
tot sys per point   64  (those missing will be set to ~0)


*********************************************/


#include "CMS.h"

void CMS_1JET_8TEVFilter::ReadData()
{

  // Opening files
  fstream rS, rCorr;

  //Bins specification
  int nbins = 6;
  double  y[]  = {0.25, 0.75, 1.25, 1.75, 2.25, 2.75 };
  int ndata[]  = {46, 46, 45, 41, 34, 27};    
  int nstat[] =  {37, 37, 36, 32, 25, 18};   

  int realsys = 27;  // number of real systematics
  int n = 0;         // count total number of datapoints


  for(int bin=0; bin < nbins; bin++ )
  {

    cout << " bin ok --------- " << bin << " ------ " << endl;   
    // start filter
    string data_file = "/CMS_8TeV_jets_Ybin" + to_string(bin+1) + ".dat";
    string cov_file =  "/CMS_8TeV_jets_Ybin" + to_string(bin+1) + "___CMS_8TeV_jets_Ybin" + to_string(bin+1) + ".dat";

    // data files
    stringstream DataFile("");
    DataFile << dataPath() << "rawdata/" << fSetName << data_file;
    rS.open(DataFile.str().c_str(), ios::in);
    if (rS.fail()) {
      cerr << "Error opening data file " << DataFile.str() << endl;
      exit(-1);
    }

    // statistics correlation matrix
    stringstream DataFileCorr("");
    DataFileCorr << dataPath() << "rawdata/" << fSetName << cov_file;
    rCorr.open(DataFileCorr.str().c_str(), ios::in);

    if (rCorr.fail()) {
      cerr << "Error opening data file " << DataFileCorr.str() << endl;
      exit(-1);
    }
  
    string line;
    double ptmin, ptmax, dum;

    // pT distribution in the rapidity bin
  
    // skip the first 42 lines
    getline(rS,line);
    for (int i = 0; i < 42; i++)  
      getline(rS,line);

    for (int i = n; i < n + ndata[bin]; i++)
    {
      cout << " bin  " << bin << " ------ " << endl;
      getline(rS,line);                      // rS reads all a line of the file in a single go
      istringstream lstream(line);          // the line is split in its columns, so that rS reads a single value in each go

      lstream >> dum >> dum >> dum;
      lstream >> ptmin >> ptmax;
      fKin1[i] = (ptmin + ptmax)*0.5;       // pt, central value of the bin
      fKin2[i] = y[bin];                    // y, central value of the bin
      fKin3[i] = 8000;                         // sqrt(s)

      lstream >> fData[i] >> dum;    
      cout << i << " " << fData[i] << endl;
    
      //** read first systematics here
      //** read statystic here

      //then read the other systematics in a loop
      // put the real systematics first. Rememeber they are expressed as percentage
      for (int l = 0; l < realsys -1; l++)     
      {

        cout << " sys ok --------- " << bin << " ------ " << endl;
        lstream >> fSys[i][l].add;  
        fSys[i][l].mult = fSys[i][l].add/fData[i]*1e2;
        fSys[i][l].type = MULT;  
        fSys[i][l].name = "CORR";
      }


      // Lumi uncertainty
      lstream >> fSys[i][realsys -1].add;  
      //fSys[i][fNSys-1].add *= pb2fb;
      fSys[i][realsys -1].mult = fSys[i][realsys -1].add/fData[i]*1e2;
      fSys[i][realsys -1].type = MULT;
      fSys[i][realsys -1].name = "CMSLUMI12";

      // Statistical error
      fStat[i] = 0;
      lstream >> fStat[i];      // no, the sist is the one called abs_sis read before. this is the percenatge value probably
      lstream >> dum;           // last column is a dummy variable....?  yes
      cout << " totstat ok --------- " << bin << " ------ " << endl;   
    }


    
    // Defining covariance matrix for the specific bin
    double** covmat = new double*[nstat[bin]];
    for (int i = 0; i < nstat[bin]; i++) 
      covmat[i] = new double[nstat[bin]];

    // Reading Covariance Matrix
    for (int i = 0; i < 16; i++)
      getline(rCorr,line);


    for (int i = 0; i < nstat[bin]; i++){ 
 
      for (int j = 0; j < nstat[bin]; j++) { 
        getline(rCorr,line);                       
        rCorr >> dum >> dum >> dum >> dum;
        rCorr >> covmat[i][j];
        //covmat[j][i]=covmat[j][i]*totstat[j]*totstat[i];    
      }
      getline(rCorr,line);                       
    }

    for (int i = 0; i < nstat[bin]; i++){ 
 
      for (int j = 0; j < nstat[bin]; j++) {    

        cout << covmat[i][j] << " ";

      }
      cout << endl;
    }

    if (bin==4) break;

    
    cout << " reading ok --------- " << bin << " ------ " << endl;   
    // Generate artificial systematics
    double** syscor = new double*[nstat[bin]];
    for(int i = 0; i < nstat[bin]; i++)
      syscor[i] = new double[nstat[bin]];

    cout << " reading ok --------- " << bin << " ------ " << endl;   
    if(!genArtSys(nstat[bin],covmat,syscor))
    {
       cerr << " in " << fSetName << " : cannot generate artificial systematics" << endl;
       exit(-1);
    }
    cout << " generated ok --------- " << bin << " ------ " << endl;   
  
    // Copy the artificial systematics in the fSys matrix
    for (int i = n; i < n + nstat[bin]; i++)
      for (int l = realsys; l < realsys + nstat[bin]; l++)     // -1 is due to the fact that i have already put the lumninosity uncertainty
      {
        cout << "qui " << endl;
        fSys[i][l].add  = syscor[i-n][l-n];
        fSys[i][l].mult = fSys[i][l].add/fData[i]*1e2;
        fSys[i][l].type = MULT;  
        fSys[i][l].name = "CORR";
      }

    for(int i = 0; i < nstat[bin]; i++) 
      delete[] syscor[i];
    delete[] syscor;

    
    for(int i = 0; i < nstat[bin]; i++) 
      delete[] covmat[i];
    delete[] covmat;
    

    rS.close();
    rCorr.close();

    n += ndata[bin];
  } 
   
}

