/**
 * ZEUS experiment -> HERA-II run
 *
 * The ZEUS HERA-II experiment is divided into 2 sets
 *
 * - ZEUS06NC
 *
 *     ZEUS06.f
 *
 *     Reduced NC and CC cross sections from HERA-II ZEUS
 *
 *     Remember that in fortran, these files with "include" are
 *     not automatically recompiled when they are modified
 *
 *     For the NC dataset, the absolute normalization uncertainty is
 *     2.6% (See Pag. 11 of DESY-08-202)
 *
 *     For the CC dataset, the same absolute normalization due to
 *     the luminosity is quoted, 2.6%, and should be completely
 *     correlated between the two data sets (check!)
 */

#include "ZEUS2.h"

void Z06NCFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile1("");
  datafile1 << dataPath() << "rawdata/"
  << fSetName << "/ZEUS06NC.sys";
  f1.open(datafile1.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile1.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  double up, down, stmp, dtmp, shift, tmp, uncsys;
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin2[i];   //Q2
    lstream >> fKin1[i];   //x
    fKin3[i] = fKin2[i]/fKin1[i]/(318.0*318.0);   //y
    
    lstream >> fData[i];

    shift = 0;    
    lstream >> up >> down;
    symmetriseErrors(up, -down, &stmp, &dtmp);
    fStat[i] = stmp;
    shift += dtmp;
    fStat[i]*= fData[i]*1e-2;
        
    lstream >> tmp >> tmp; //systot

    fSys[i][0].mult = 2.6;  //Luminosity
    fSys[i][0].type = MULT;
    fSys[i][0].name = "Z06LUMI";

    fSys[i][1].mult = 4.0;  //4% polarization uncertainty
    fSys[i][1].type = MULT;
    fSys[i][1].name = "CORR";
    
    for (int j = 2; j < 9; j++)              //Systematics (asymmetric)
    {
      lstream >> up >> down;
      symmetriseErrors(up, -down, &stmp, &dtmp);
      shift += dtmp;
      fSys[i][j].mult = stmp;
      fSys[i][j].type = ADD;
      fSys[i][j].name = "CORR";
    }

    lstream >> up >> down;                     //Uncorrelated systematic
    symmetriseErrors(up, -down, &uncsys, &dtmp);    
    shift += dtmp;    
    fSys[i][9].mult = uncsys;
    fSys[i][9].type = ADD;
    fSys[i][9].name = "UNCORR";

    for (int j = 0; j < fNSys; j++)    
     fSys[i][j].add = fSys[i][j].mult*fData[i]*1e-2;   

    fData[i]*=(1.0 + shift*0.01);  //Shift from asymmetic errors
  }
  
  f1.close();
}

/**
 * ZEUS experiment -> HERA-II run
 *
 * The ZEUS HERA-II experiment is divided into 2 sets
 *
 * - ZEUS06CC
 */
void Z06CCFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/ZEUS06CC.data";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  string line;
  double up, down, shift, systot;
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> fKin2[i];   //Q2
    lstream >> fKin1[i];   //x
    fKin3[i] = fKin2[i]/fKin1[i]/(318.0*318.0);   //y
    
    lstream >> fData[i];
    
    lstream >> fStat[i];  

    fSys[i][0].mult = 2.6;  //Luminosity
    fSys[i][0].add = fSys[i][0].mult*fData[i]*1e-2;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "Z06LUMI";
    
    lstream >> up >> down;                     //systematics
    symmetriseErrors(up, -down, &systot, &shift);
    fSys[i][1].add = systot;
    fSys[i][1].mult = fSys[i][1].add *100/fData[i];
    fSys[i][1].type = ADD;
    fSys[i][1].name = "CORR";
      
    fData[i]+= shift;  //Shift from asymmetic errors
  }
  
  f1.close();
}

/**
 *     zeus-hera2-NCp.f
 *
 *     This is the data described in DESY-12-145
 *     arXiv:1208.6138
 *
 *     It includes the inclusive positron NC cross sections
 *     measured by ZEUS in the HERA-II data taking period between
 *     2006 and 2007 with proton beam energies Ep of 920 GeV
 *
 *     The measurement are provided both by posively and negatively
 *     polarized positrons beams, as well as the average over
 *     polarizations, which is what we use here.
 *
 *     For the time being we consider the data as uncorrelated to
 *     the existing ZEUS HERA-II data (2006 paper) in order
 *     to compare with the existing results.
 *     The experimentalists themselves treat those datasets as uncorrelated
 *     when they extract F3 (see top of page 13).
 *
 *     The neutral current data is available from
 *     HERAII-ZEUS-NC-ep.data
 *     Taken from the latex source of the paper.
 *     ( 22.10. 2012 Cross-checked by JR-MU )
 *
 *     The format of the file is
 *     Q2min Q2max Q2c xmin xmax xc  sigma_ref  stat
 *
 *     where the statistical error is given in absolute values,not in percent
 *
 *     The data has been taken with the center of mass energy of
 *     318 GeV, thus y = Q2(GeV)/(x*318^2)
 *
 *     The uncertainty on the total integrated luminosity is 1.8%
 *
 *     The full list of correlated systematic errors are available from
 *     http://www-zeus.desy.de/zeus_papers/ZEUS_PAPERS/DESY-12-145_sys_tables.ps
 *
 *     The systematic errors  have been tabulated from the file
 *     zeus-hera2/DESY-12-145_sys_tables.txt
 *
 *     which has been taken from the .eps source of the addendum to the paper,
 *     so it should be taken with great care. Here errors are already in %
 *     ( 22.10. 2012 Cross-checked by JR-MU )
 *
 *     We also received the original list of systematics from Trevor Stewart,
 *     the analyzer of the data, with the systematic uncertainties
 *     already symmetrized
 *     This info has been saved in the file
 *     HERAII-ZEUS-NC-ep.sys.paper
 *
 *     Which has the following format:
 *     column 1: Q2
 *     column 2: x
 *     column 3: xsec
 *     column 4: stat.
 *     column 5: the symmetrised systematic uncertainties added in quadrature.
 *     column 6-20: the symmetrised systematic uncertainties for each
 *     systematic source.
 *     column 21: the lumi uncertainty.
 *
 *
 *     Official recommendation from Mandy Cooper-Sarkar, 3/12/2012
 *
 *     YES there are some issues.
 *     In the paper itself we just say to use them uncorrelated.
 *     But now our work on the combination leads us to say that  sources 1,2,3,4,6,8,10,12,13 are point to point correlated- but not 5,7,9,11,14 and 15.
 *     Also we are now working on quantiying what are the correlations between these e+NC sources and the e-NC (HERA_Ii) sources which were published earlier. For example source 3 is correlated to source 7 from the e-p paper but with the opposite sign definition (for the table which you are using). I don't think we are finished on this yet and it may be better if you just used them uncorrelated until this is sorted out.
 *
--------------------------------------------------------
----- UPDATE MARCH 2013 --------------------------------
--------------------------------------------------------

The raw data files

HERAII-ZEUS-NC-ep.data      
HERAII-ZEUS-NC-ep.sys.paper

were the ones in the original version of the paper, published
in August 2013

In March 2013, these were updated with a corrected list of systematics,
which can be obtained from the ZEUS webpage

http://www-zeus.desy.de/zeus_papers/ZEUS_PAPERS/tables_txt_final_1304.tar

We only keep the data in the ALL folder, that is, the cross-sections
averaged over the different beam polarizations

The experimental data has been taked from Table 4 of DESY-12-145.pdf,
which gives the positron NC reduced cross setcions from ZEUS in the 
HERA-II run, corrected to zero beam polarization

We are interested in the correlated systematics of the reduced
cross section, file ALL/sys_redxsecUNPOL.txt 

The systematics are described in
http://www-zeus.desy.de/zeus_papers/ZEUS_PAPERS/DESY-12-145_sys_tables.pdf
which provides the individual systematic uncertainties,
δ1 – δ15 as deﬁned in Section 7 of DESY-12-145, separately. The ﬁrst(upper)
uncertainty listed in the tables always follows from an upward variation of a
parameter or cut. Bin-to-bin correlations were found for 
delta 1,2,3,4,6,8,10,12 and 13. Tables were updated March 2013.
In particular we need the asymmetric errors, Table 8, file
sys_redxsecUNPOL_asy.txt 

So in the final version of the buildmaster, the original
correlation matrix
HERAII-ZEUS-NC-ep.sys.paper
is replaced with the final one 
sys_redxsecUNPOL_asy.txt
 
 *
 */
void ZEUSHERA2NCPFilter::ReadData()
{
  // Opening files
  fstream f1;
 
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/sys_redxsecUNPOL_asy.txt";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Filtering data
  string line;
  double up, down, shtmp, untmp, shift;
  double systotu, systotd, tmp, sysup, sysdn;
  for (int i = 0; i < fNData; i++)
  {
    shift=0;
    getline(f1,line);
    istringstream lstream(line);
    
    lstream >> tmp >> tmp; //Q2 range
    lstream >> fKin2[i];   //Q2
    lstream >> tmp >> tmp; //x range
    lstream >> fKin1[i];   //x
    fKin3[i] = fKin2[i]/fKin1[i]/(318.0*318.0);   //y  
  
    lstream >> fData[i];   //obs
    
    lstream >> up >> down;   //stat (asymmetric)
    symmetriseErrors(up, down, &untmp, &shift);
    fStat[i] = untmp*fData[i]*1e-2;
    
    lstream >> systotu >> systotd;  //total systematic
    
    sysup=0; sysdn=0;
  // There are 15 systematic errors
  // Of which only 9 are bin by bin correlated
  // This corresponds to sources 1,2,3,4,6,8,10,12,13
    for (int j = 0; j < 15; j++)
    {
      lstream >> up >> down;        //systematics (asymmetric)
      sysup+= pow(max(up,down),2); 
      sysdn+= pow(min(up,down),2);
      untmp=0; shtmp=0;
      symmetriseErrors(up, down, &untmp, &shtmp);
      shift += shtmp;
      fSys[i][j].mult = untmp;
      fSys[i][j].add = fSys[i][j].mult*fData[i]*1e-2;
      fSys[i][j].type = ADD;
      fSys[i][j].name = "CORR";
    }
    
    // some systematics are uncorrelated (5,7,11,14,15)
    fSys[i][4].name = "UNCORR";
    fSys[i][6].name = "UNCORR";
    fSys[i][8].name = "UNCORR";
    fSys[i][10].name = "UNCORR";
    fSys[i][13].name = "UNCORR";
    fSys[i][14].name = "UNCORR";
    
    // shift datapoints
    fData[i]*=(1.0+shift*1e-2);
    
    lstream >> fSys[i][15].mult;  //lumi
    fSys[i][15].add = fSys[i][15].mult*fData[i]*1e-2;
    fSys[i][15].type = MULT;
    fSys[i][15].name = "ZEUS2LUMI2";
    
    // Check total systematic error
    sysup=sqrt(sysup); sysdn=sqrt(sysdn);
    if(fabs(1-sysup/systotu) > 1e-1)
 	{  
 	  cerr << "ZEUSHERA2NCP: Mismatch in total systematic uncertainty (+) "
	   << i+1 << "\t" << sysup << "\t" << systotu << endl;
	} 
	if(fabs(1+sysdn/systotd) > 1e-1)
 	{  
 	  cerr << "ZEUSHERA2NCP: Mismatch in total systematic uncertainty (-) "
	   << i+1 << "\t" << sysdn << "\t" << -systotd << endl;
	}  
    
  }
  
  f1.close();
}

/**
 *
 *     zeus-hera2-CCp.f
 *
 *     EPJ C 70, Issue 4 (2010) 945-963
 *     DESY-10-129
 *     This is the data described in DESY-12-145
 *     arXiv:1208.6138
 *
 *     It includes the inclusive positron CC cross sections
 *     measured by ZEUS in the HERA-II data taking period between
 *     2006 and 2007 with proton beam energies Ep of 920 GeV
 *
 *     The measurement are provided both by posively and negatively
 *     polarized positrons beams, as well as the average over
 *     polarizations, which is what we use here
 *
 *     For the time being we consider the data as uncorrelated to
 *     the existing ZEUS HERA-II data (2006 paper) in order
 *     to compare with the existing results
 *
 *     The charged current data is available from
 *      HERAII-ZEUS-CC-ep.data
 *     (Taken from the latex source of the paper, to cross-check)
 *
 *     The format of the file is
 *     Q2  x (adum,i=1,12)  sigma_ref +stat -stat +sys -sys factor
 *
 *     where the overall factor multiplies everything,
 *     and the uncertainties are given in absolute value
 *     not in percent
 *
 *     The data has been taken with the center of mass energy of
 *     318 GeV, thus y = Q2(GeV)/(x*318^2)
 *
 *     Note that for the helicity averaged measurement the separated
 *     correlated systematic errors are not available
 *     so we are forced to add in quadrature the statistical and
 *     systematic error.
 *     This should be corrected once the combined HERA-II dataset is available
 *
 *     The uncertainty on the total integrated luminosity is 2.6%
 *
 */
void ZEUSHERA2CCPFilter::ReadData()
{
  // Opening files
  fstream f1;
  
  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/HERAII-ZEUS-CCp.data";
  f1.open(datafile.str().c_str(), ios::in);
  
  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }
  
  // Starting filter
  double datain[fNData][20];
  
  // Reading data
  for (int i = 0; i < fNData; i++)
    for (int j = 0; j < 20; j++)
      f1 >> datain[i][j];
  
  // Filtering data
  double stat[fNData], x[fNData], y[fNData];
  double q2[fNData], rxsec[fNData], shift[fNData];
  double Ecm = 318.0; // GeV, center of mass energy
  
  for (int i = 0; i < fNData; i++)
  {
    x[i]   = datain[i][1];
    q2[i]  = datain[i][0];
    y[i]   = q2[i]/x[i]/pow(Ecm,2.0);
    
    if (y[i] >= 1.0)
    {
      cerr << "Invalid value of y" << endl;
      exit(-1);
    }
    
    // Observable
    rxsec[i] = datain[i][14]*datain[i][19]; // Overall factor
    shift[i] = 0.0;
    
    // Statistical errors
    double up = datain[i][15]*datain[i][19];
    double dn = datain[i][16]*datain[i][19];
    double stmp = 0, dtmp = 0;
    symmetriseErrors(up,dn,&stmp,&dtmp);
    stat[i] = stmp;
    shift[i] += dtmp;
    
    // Uncorrelated systematics
    up = datain[i][17]*datain[i][19];
    dn = datain[i][18]*datain[i][19];
    stmp = dtmp = 0;
    symmetriseErrors(up,dn,&stmp,&dtmp);
    double systot = stmp;
    shift[i] += dtmp;
    
    fKin1[i] = x[i];
    fKin2[i] = q2[i];
    fKin3[i] = y[i];
    
    fData[i] = rxsec[i];
    fStat[i] = stat[i];

    // Normalisation uncertainty
    fSys[i][0].mult = 2.6;
    fSys[i][0].add = fSys[i][0].mult*fData[i]*1e-2;
    fSys[i][0].type = MULT;
    fSys[i][0].name = "ZEUS2LUMI2";
    
    // Systematics
    fSys[i][1].add = systot;
    fSys[i][1].mult = fSys[i][1].add*100/fData[i]; 
    fSys[i][1].type = ADD;
    fSys[i][1].name = "UNCORR"; 
    
    // Shift in the data due to asymmetric errors
    fData[i] += shift[i]; 
  }
  
  f1.close();
}

