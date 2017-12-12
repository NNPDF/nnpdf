/**
 * ATLAS JETS 2010 - 1112.6297v1
 *
 * Inclusive jet and dijet data from the ATLAS 2010 dataset
 * Statistics from 36 pb^-1 dataset
 *
 * There are 18 sources of correlated systematics, that
 * for which each rapidity bin is fully correlated in pt,
 * but only some bins in rapidity might be correlated among them
 * So effectively there are 86 source of systematic uncertainties
 * fully correlated between all pt and eta
 *
 * See Note added on 27 Jun 2014 at http://hepdata.cedar.ac.uk/view/ins1082936
 * updated luminosity
 */

#include "ATLAS.h"

void ATLASR04JETS36PBFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3, f4, f5;

  stringstream hepdata("");
  hepdata << dataPath() << "rawdata/"
  << fSetName << "/ATLAS-jets-R04-36pb-hepdata.data";
  f1.open(hepdata.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << hepdata.str() << endl;
    exit(-1);
  }

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/ATLAS-jets-R04-36pb.data";
  f2.open(datafile.str().c_str(), ios::in);

  if (f2.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream covfile("");
  covfile << dataPath() << "rawdata/"
  << fSetName << "/ATLAS-jets-36pb.sysmat";
  f3.open(covfile.str().c_str(), ios::in);

  if (f3.fail()) {
    cerr << "Error opening data file " << covfile.str() << endl;
    exit(-1);
  }

  stringstream covfile2("");
  covfile2 << dataPath() << "rawdata/"
  << fSetName << "/ATLAS-jets-R04-36pb.sys";
  f4.open(covfile2.str().c_str(), ios::in);

  if (f4.fail()) {
    cerr << "Error opening data file " << covfile.str() << endl;
    exit(-1);
  }

  stringstream nptcorr("");
  nptcorr << dataPath() << "rawdata/"
  << fSetName << "/ATLAS-jets-R04-36pb.npt";
  f5.open(nptcorr.str().c_str(), ios::in);

  if (f5.fail()) {
    cerr << "Error opening data file " << nptcorr.str() << endl;
    exit(-1);
  }

  // Reading hebdata file
  string line;
  double tmp,pt,sysup,sysdown,sysuncor;
  double s = 7000;
  double lcorr = 1.0187; // correction factor due to luminosity upgrade
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);

    lstream >> pt;
    lstream >> tmp >> tmp;         //ptmax and ptmin

    fKin2[i] = pt*pt;              //pt2
    fKin3[i] = s;                  //sqrt(s)

    lstream >> fData[i];
    fData[i] *= lcorr; // apply lumi correction

    lstream >> fStat[i] >> tmp;    //statistical uncertainty (symmetric)
    fStat[i] *= lcorr; // apply lumi correction

    lstream >> sysup >> sysdown;   //total (correlated) systematic uncertainty (asymmetric)
    lstream >> sysuncor >> tmp;    //uncorrelated systematic uncertainty (symmetric)
  }

  // Read data file - need rapidity bins to read systematics file
  const int nrapbin = 7;   // There are seven rapidity bins
  double etamin, etamax, eta;
  int ndatbin[nrapbin];

  int idat = 0;
  getline(f2,line);
  for (int i = 0; i < nrapbin; i++)
  {
    getline(f2,line);
    istringstream lstream(line);
    lstream >> etamin >> etamax;
    eta = (etamax+etamin)*0.5;

    getline(f2,line);
    istringstream lstream2(line);
    lstream2 >> ndatbin[i];
    for (int j = 0; j < ndatbin[i]; j++)
    {
      fKin1[idat] = eta;                   //eta
      getline(f2,line);                   //Rest of values in file are already in the hepdata file
      idat++;
    }
  }

  // Reading systematic correlations
  const int nsystype = 18; // There are 18 systematics per rapidity bin
  int sysnumber[nrapbin][nsystype];
  string tmps;

  for (int isys = 0; isys < nsystype; isys++)
  {
    getline(f3,line);
    istringstream lstream(line);

    //Need to avoid the strings at front of each line
    bool atnumbers = false;
    while(atnumbers != true) {
      lstream >> tmps;
      if (tmps[tmps.size()-1] == '\"') atnumbers = true;
    }

    for (int j = 0; j < nrapbin; j++)
    {
      lstream >> sysnumber[j][isys];      //Read in number identifing systematic
      sysnumber[j][isys] += 4;            //Save [0] for luminosity, [1] nonperturbative uncertainty, [2]-[4] for uncorrelated systematics
    }
  }

  //Reading systematics
  double up,down,shift[fNData],stmp,dtmp;
  idat = 0;
  for (int i = 0; i < nrapbin; i++)
  {
    getline(f4,line);
    for (int j = 0; j < ndatbin[i]; j++)
    {
      getline(f4,line);
      istringstream lstream(line);
      lstream >> tmp >> tmp;      //ptmax and ptmin

      shift[idat] = 0;
      for (int isystype = 0; isystype < nsystype; isystype++)
      {
        lstream >> up >> down;
        symmetriseErrors(up,down,&stmp,&dtmp);
        shift[idat]+=dtmp;
        fSys[idat][sysnumber[i][isystype]].mult=stmp;   //systematics are in %
        fSys[idat][sysnumber[i][isystype]].name = "CORR";
      }

      for (int j = 2; j < 5; j++)
      {
        lstream >> fSys[idat][j].mult >> tmp;   //uncorrelated systematics are symmetric (%)
        fSys[idat][j].name = "UNCORR";
      }
      idat++;
    }
  }

  //Reading nonperturbative correction and associated uncertainty
  double nonper, error;
  for (int i = 0; i < fNData; i++)
  {
    getline(f5,line);
    istringstream lstream(line);

    lstream >> nonper >> error;

    fData[i]/= nonper;            //Multiplicative correction to theory so instead divide out of data
    fStat[i]/= nonper;            //Do the same to statistical uncertainty stored in absolute form

    fSys[i][0].mult = 3.5;           //Luminosity uncertainty of 3.5% (updated in 27 JUN 2014)
    fSys[i][0].name = "ATLASLUMI10";

    fSys[i][1].mult = error/nonper*100;    // Error on nonperturbative correction (%)
    fSys[i][1].name = "CORR";
  }

  for (int i=0; i<fNData; i++)
    for (int l = 0; l < fNSys; l++)
    {
      fSys[i][l].add = fSys[i][l].mult*fData[i]*1e-2;    // Calculate absolute systematics for additive part
      fSys[i][l].type = MULT;                            // All systematics multiplicative
    }

  f1.close();
  f2.close();
  f3.close();
  f4.close();
  f5.close();
}

/**
 * ATLAS JETS 2010 - 1112.6297v1
 *
 * Inclusive jet and dijet data from the ATLAS 2010 dataset
 * Statistics from 36 pb^-1 dataset, R = 0.6
 *
 * See Note added on 27 Jun 2014 at http://hepdata.cedar.ac.uk/view/ins1082936
 * updated luminosity
 */
void ATLASR06JETS36PBFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3, f4, f5;

  stringstream hepdata("");
  hepdata << dataPath() << "rawdata/"
  << fSetName << "/ATLAS-jets-R06-36pb-hepdata.data";
  f1.open(hepdata.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << hepdata.str() << endl;
    exit(-1);
  }

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/ATLAS-jets-R06-36pb.data";
  f2.open(datafile.str().c_str(), ios::in);

  if (f2.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream covfile("");
  covfile << dataPath() << "rawdata/"
  << fSetName << "/ATLAS-jets-36pb.sysmat";
  f3.open(covfile.str().c_str(), ios::in);

  if (f3.fail()) {
    cerr << "Error opening data file " << covfile.str() << endl;
    exit(-1);
  }

  stringstream covfile2("");
  covfile2 << dataPath() << "rawdata/"
  << fSetName << "/ATLAS-jets-R06-36pb.sys";
  f4.open(covfile2.str().c_str(), ios::in);

  if (f4.fail()) {
    cerr << "Error opening data file " << covfile.str() << endl;
    exit(-1);
  }

  stringstream nptcorr("");
  nptcorr << dataPath() << "rawdata/"
  << fSetName << "/ATLAS-jets-R06-36pb.npt";
  f5.open(nptcorr.str().c_str(), ios::in);

  if (f5.fail()) {
    cerr << "Error opening data file " << nptcorr.str() << endl;
    exit(-1);
  }

  // Reading hebdata file
  string line;
  double tmp,pt,sysup,sysdown,sysuncor;
  double s = 7000;
  double lcorr = 1.0187; // correction factor due to luminosity upgrade
  for (int i = 0; i < fNData; i++)
  {
    getline(f1,line);
    istringstream lstream(line);

    lstream >> pt;
    fKin2[i] = pt*pt;              //pt2
    lstream >> tmp >> tmp;         //ptmax and ptmin

    fKin3[i] = s;                  //sqrt(s)

    lstream >> fData[i];
    fData[i] *= lcorr; // apply lumi correction

    lstream >> fStat[i] >> tmp;    //statistical uncertainty (symmetric)
    fStat[i] *= lcorr; // apply lumi correction

    lstream >> sysup >> sysdown;   //total (correlated) systematic uncertainty (asymmetric)

    lstream >> sysuncor >> tmp;    //uncorrelated systematic uncertainty (symmetric)
  }

  // Read data file - need rapidity bins to read systematics file
  const int nrapbin = 7;   // There are seven rapidity bins
  double etamin, etamax, eta;
  int ndatbin[nrapbin];

  int idat = 0;
  getline(f2,line);
  for (int i = 0; i < nrapbin; i++)
  {
    getline(f2,line);
    istringstream lstream(line);
    lstream >> etamin >> etamax;
    eta = (etamax+etamin)*0.5;

    getline(f2,line);
    istringstream lstream2(line);
    lstream2 >> ndatbin[i];
    for (int j = 0; j < ndatbin[i]; j++)
    {
      fKin1[idat] = eta;                   //eta
      getline(f2,line);                   //Rest of values in file are already in the hepdata file
      idat++;
    }
  }

  // Reading systematic correlations
  const int nsystype = 18; // There are 18 systematics per rapidity bin
  int sysnumber[nrapbin][nsystype];
  string tmps;

  for (int isys = 0; isys < nsystype; isys++)
  {
    getline(f3,line);
    istringstream lstream(line);

    //Need to avoid the strings at front of each line
    bool atnumbers = false;
    while(atnumbers != true) {
      lstream >> tmps;
      if (tmps[tmps.size()-1] == '\"') atnumbers = true;
    }

    for (int j = 0; j < nrapbin; j++)
    {
      lstream >> sysnumber[j][isys];      //Read in number identifing systematic
      sysnumber[j][isys] += 4;            //Save [0] for luminosity, [1] nonperturbative uncertainty, [2]-[4] for uncorrelated systematics
    }
  }

  //Reading systematics
  double up,down,shift[fNData],stmp,dtmp;
  idat = 0;
  for (int i = 0; i < nrapbin; i++)
  {
    getline(f4,line);
    for (int j = 0; j < ndatbin[i]; j++)
    {
      getline(f4,line);
      istringstream lstream(line);
      lstream >> tmp >> tmp;      //ptmax and ptmin

      shift[idat] = 0;
      for (int isystype = 0; isystype < nsystype; isystype++)
      {
        lstream >> up >> down;
        symmetriseErrors(up,down,&stmp,&dtmp);
        shift[idat]+=dtmp;
        fSys[idat][sysnumber[i][isystype]].mult=stmp;   //systematics are in %
        fSys[idat][sysnumber[i][isystype]].name = "CORR";
      }

      for (int j = 2; j < 5; j++)
      {
        lstream >> fSys[idat][j].mult >> tmp;   //uncorrelated systematics are symmetric (%)
        fSys[idat][j].name = "UNCORR";
      }
      idat++;
    }
  }

  //Reading nonperturbative correction and associated uncertainty
  double nonper, error;
  for (int i = 0; i < fNData; i++)
  {
    getline(f5,line);
    istringstream lstream(line);

    lstream >> nonper >> error;

    fData[i]/= nonper;            //Multiplicative correction to theory so instead divide out of data
    fStat[i]/= nonper;            //Do the same to statistical uncertainty stored in absolute form

    fSys[i][0].mult = 3.5;           //Luminosity uncertainty of 3.5% (updated in 27 JUN 2014)
    fSys[i][0].name = "ATLASLUMI10";

    fSys[i][1].mult = error/nonper*100;    //Error on nonperturbative correction (%)
    fSys[i][1].name = "CORR";
  }

  for (int i=0; i<fNData; i++)
    for (int l = 0; l < fNSys; l++)
    {
      fSys[i][l].add = fSys[i][l].mult*fData[i]*1e-2;    // Calculate absolute systematics for additive part
      fSys[i][l].type = MULT;                            // All systematics multiplicative
    }

  f1.close();
  f2.close();
  f3.close();
  f4.close();
  f5.close();
}

/**
 * ATLASR04JETS2P76TEV filter
 * 1304.4739v1
 *
 */
void ATLASR04JETS2P76TEVFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3, f4, f5;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/ATLASR04JETS2P76TEV_raw.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream covfile("");
  covfile << dataPath() << "rawdata/"
  << fSetName << "/ATLASR04JETS2P76TEV_raw_SYSTYPE.data";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail()) {
    cerr << "Error opening data file " << covfile.str() << endl;
    exit(-1);
  }

  // Reading hebdata file
  string line;
  double s = 2760;

  // Reading systematic correlations
  const int nrapbin = 7;   // There are seven rapidity bins
  const int nsystype = 19; // There are 19 systematics per rapidity bin (ignoring luminosity for the moment)
  int sysnumber[nrapbin][nsystype];
  string tmps;

  for (int isys = 0; isys < nsystype; isys++)
  {
    getline(f2,line);
    istringstream lstream(line);

    //Need to avoid the strings at front of each line
    bool atnumbers = false;
    while(atnumbers != true) {
      lstream >> tmps;
      if (tmps[tmps.size()-1] == '\"') atnumbers = true;
    }

    for (int j = 0; j < nrapbin; j++)
    {
      lstream >> sysnumber[j][isys];      //Read in number identifing systematic
      sysnumber[j][isys] += 3;            //Save [0] for luminosity, [1] nonperturbative uncertainty, [2]-[3] for uncorrelated systematics
    }
  }

  //Reading data
  double etamin, etamax, eta, ptmax, ptmin, nonpert;
  double up,down,shift,stmp,dtmp,sys88,sys89;
  int ndatbin[nrapbin];
  int idat = 0;
  for (int i = 0; i < nrapbin; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    lstream >> etamin >> etamax;
    eta = (etamax+etamin)*0.5;

    getline(f1,line);
    istringstream lstream2(line);
    lstream2 >> ndatbin[i];
    for (int j = 0; j < ndatbin[i]; j++)
    {
      getline(f1,line);
      istringstream lstream(line);

      fKin1[idat] = eta;                   //eta

      lstream >> ptmin >> ptmax;           //ptmin and ptmax
      fKin2[idat] = (ptmin+ptmax)*0.5;     //pt
      fKin2[idat] *= fKin2[idat];          //pt2

      fKin3[idat] = s;                     //sqrt(s)

      fSys[idat][0].mult = 2.7;      //2.7% luminosity uncertainty (uncorrelated with 36PB luminosity uncertainty)
      fSys[idat][0].name = "CORR";

      shift = 0;

      lstream >> nonpert >> up >> down;   //nonperturbative correction and associated uncertainty
      up *= 100/nonpert;                  //convert to percentage
      down *= 100/nonpert;
      symmetriseErrors(up,down,&stmp,&dtmp);
      shift=dtmp;
      fSys[idat][1].mult=stmp;
      fSys[idat][1].name = "CORR";

      lstream >> fData[idat];
      fData[idat] /= nonpert;               //apply nonperturbative correction
                                            //note: multiplicative correction so additive systematics should be calcuated with corrected data

      lstream >> fStat[idat];                //statistical uncertainty in %
      fStat[idat] *= fData[idat]*0.01;  //convert to absolute

      //First 5 systematics
      for (int isystype = 0; isystype < 5; isystype++)
      {
        lstream >> up >> down;
        symmetriseErrors(up,down,&stmp,&dtmp);
        shift+=dtmp;
        fSys[idat][sysnumber[i][isystype]].mult=stmp;   //systematics are in %
        fSys[idat][sysnumber[i][isystype]].name = "CORR";
      }

      //Systematics number 88 and 89 are correlated in quadrature with 31 from ATLASJETS36PB
      lstream >> up >> down;
      symmetriseErrors(up,down,&stmp,&dtmp);
      shift+=dtmp;
      sys88 = stmp;   //systematics are in %

      //Next 7 systematics
      for (int isystype = 6; isystype < 13; isystype++)
      {
        lstream >> up >> down;
        symmetriseErrors(up,down,&stmp,&dtmp);
        shift+=dtmp;
        fSys[idat][sysnumber[i][isystype]].mult=stmp;   //systematics are in %
        fSys[idat][sysnumber[i][isystype]].name = "CORR";
      }

      //Systematics number 88 and 89 are correlated in quadrature with 31 from ATLASJETS36PB
      lstream >> up >> down;
      symmetriseErrors(up,down,&stmp,&dtmp);
      shift+=dtmp;
      sys89 = stmp;                                            //systematics are in %
      fSys[idat][sysnumber[6][4]+1].mult=sqrt(sys88*sys88+sys89*sys89);  //in quadrature for sys 31 of this set
      fSys[idat][sysnumber[6][4]+1].name = "CORR";

      //Next 5 systematics - symmetric
      for (int isystype = 14; isystype < nsystype; isystype++)
      {
        lstream >> fSys[idat][sysnumber[i][isystype]].mult;
        fSys[idat][sysnumber[i][isystype]].name = "CORR";
      }

      //Two sources of uncorrelated systematics (%)
      lstream >> fSys[idat][2].mult;
      fSys[idat][2].name = "UNCORR";

      lstream >> fSys[idat][3].mult;
      fSys[idat][3].name = "UNCORR";

      idat++;
    }
  }

  // Final processing
  for (int i=0; i<fNData; i++)
    for (int l = 0; l < fNSys; l++)
    {
      fSys[i][l].add = fSys[i][l].mult*fData[i]*1e-2;    // Calculate absolute systematics for additive part
      fSys[i][l].type = MULT;                            // All systematics multiplicative
    }

  f1.close();
  f2.close();
}


/**
 * ATLASR06JETS2P76TEV filter
 * 1304.4739v1
 *
 */
void ATLASR06JETS2P76TEVFilter::ReadData()
{
  // Opening files
  fstream f1, f2, f3, f4, f5;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/ATLASR06JETS2P76TEV_raw.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream covfile("");
  covfile << dataPath() << "rawdata/"
  << fSetName << "/ATLASR06JETS2P76TEV_raw_SYSTYPE.data";
  f2.open(covfile.str().c_str(), ios::in);

  if (f2.fail()) {
    cerr << "Error opening data file " << covfile.str() << endl;
    exit(-1);
  }

  // Reading hebdata file
  string line;
  double s = 2760;

  // Reading systematic correlations
  const int nrapbin = 7;   // There are seven rapidity bins
  const int nsystype = 19; // There are 19 systematics per rapidity bin (ignoring luminosity for the moment)
  int sysnumber[nrapbin][nsystype];
  string tmps;

  for (int isys = 0; isys < nsystype; isys++)
  {
    getline(f2,line);
    istringstream lstream(line);

    //Need to avoid the strings at front of each line
    bool atnumbers = false;
    while(atnumbers != true) {
      lstream >> tmps;
      if (tmps[tmps.size()-1] == '\"') atnumbers = true;
    }

    for (int j = 0; j < nrapbin; j++)
    {
      lstream >> sysnumber[j][isys];      //Read in number identifing systematic
      sysnumber[j][isys] += 3;            //Save [0] for luminosity, [1] nonperturbative uncertainty, [2]-[3] for uncorrelated systematics
    }
  }

  //Reading data
  double etamin, etamax, eta, ptmax, ptmin, nonpert;
  double up,down,shift,stmp,dtmp,sys88,sys89;
  int ndatbin[nrapbin];
  int idat = 0;
  for (int i = 0; i < nrapbin; i++)
  {
    getline(f1,line);
    istringstream lstream(line);
    lstream >> etamin >> etamax;
    eta = (etamax+etamin)*0.5;

    getline(f1,line);
    istringstream lstream2(line);
    lstream2 >> ndatbin[i];
    for (int j = 0; j < ndatbin[i]; j++)
    {
      getline(f1,line);
      istringstream lstream(line);

      fKin1[idat] = eta;                   //eta

      lstream >> ptmin >> ptmax;           //ptmin and ptmax
      fKin2[idat] = (ptmin+ptmax)*0.5;     //pt
      fKin2[idat] *= fKin2[idat];          //pt2

      fKin3[idat] = s;                     //sqrt(s)

      fSys[idat][0].mult = 2.7;      //2.7% luminosity uncertainty (uncorrelated with 36PB luminosity uncertainty)
      fSys[idat][0].name = "CORR";

      shift = 0;

      lstream >> nonpert >> up >> down;   //nonperturbative correction and associated uncertainty
      up *= 100/nonpert;                  //convert to percentage
      down *= 100/nonpert;
      symmetriseErrors(up,down,&stmp,&dtmp);
      shift=dtmp;
      fSys[idat][1].mult=stmp;
      fSys[idat][1].name = "CORR";

      lstream >> fData[idat];
      fData[idat] /= nonpert;               //apply nonperturbative correction
                                            //note: multiplicative correction so additive systematics should be calcuated with corrected data

      lstream >> fStat[idat];                //statistical uncertainty in %
      fStat[idat] *= fData[idat]*0.01;  //convert to absolute

      //First 5 systematics
      for (int isystype = 0; isystype < 5; isystype++)
      {
        lstream >> up >> down;
        symmetriseErrors(up,down,&stmp,&dtmp);
        shift+=dtmp;
        fSys[idat][sysnumber[i][isystype]].mult=stmp;   //systematics are in %
        fSys[idat][sysnumber[i][isystype]].name = "CORR";
      }

      //Systematics number 88 and 89 are correlated in quadrature with 31 from ATLASJETS36PB
      lstream >> up >> down;
      symmetriseErrors(up,down,&stmp,&dtmp);
      shift+=dtmp;
      sys88 = stmp;   //systematics are in %

      //Next 7 systematics
      for (int isystype = 6; isystype < 13; isystype++)
      {
        lstream >> up >> down;
        symmetriseErrors(up,down,&stmp,&dtmp);
        shift+=dtmp;
        fSys[idat][sysnumber[i][isystype]].mult=stmp;   //systematics are in %
        fSys[idat][sysnumber[i][isystype]].name = "CORR";
      }

      //Systematics number 88 and 89 are correlated in quadrature with 31 from ATLASJETS36PB
      lstream >> up >> down;
      symmetriseErrors(up,down,&stmp,&dtmp);
      shift+=dtmp;
      sys89 = stmp;                                            //systematics are in %
      fSys[idat][sysnumber[6][4]+1].mult=sqrt(sys88*sys88+sys89*sys89);  //in quadrature for sys 31 of this set
      fSys[idat][sysnumber[6][4]+1].name = "CORR";

      //Next 5 systematics - symmetric
      for (int isystype = 14; isystype < nsystype; isystype++)
      {
        lstream >> fSys[idat][sysnumber[i][isystype]].mult;
        fSys[idat][sysnumber[i][isystype]].name = "CORR";
      }

      //Two sources of uncorrelated systematics (%)
      lstream >> fSys[idat][2].mult;
      fSys[idat][2].name = "UNCORR";

      lstream >> fSys[idat][3].mult;
      fSys[idat][3].name = "UNCORR";

      for (int i=0; i<fNData; i++)
        for (int l = 0; l < fNSys; l++)
        {
          fSys[i][l].add = fSys[i][l].mult*fData[i]*1e-2;    // Calculate absolute systematics for additive part
          fSys[i][l].type = MULT;                            // All systematics multiplicative
        }

      idat++;
    }
  }

  f1.close();
  f2.close();
}

/**
 * atlas-wz
 *
 * This data is taken from the paper arxiv:1109.5141
 *
 * It contains the combined lepton rapidity distributions
 * for W+, W- and Z production, with the full covariance matrix
 * available
 *
 * Fit separatelt W+ and W- with covariance matrix instead of
 * the asymmetry distribution
 *
 * The that the results are rather close to a direct
 * refitting with the asymmetry
 * distribution
 *
 * This data will be included using reweighting, with DYNNLO
 * for the computation of the theory predictions
 *
 * Put all data together into a single set
 * This is due to the strong correlation between the W+, W^- and
 * Z0 which makes not suitable to treat them separately
 * in the first stages of the minimization
 */
void ATLASWZRAP36PBFilter::ReadData()
{
  // Opening files
  fstream fWZ[3];

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/ATLAS-36pb-Wplrap.data";
  fWZ[0].open(datafile.str().c_str(), ios::in);

  if (fWZ[0].fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream datafile2("");
  datafile2 << dataPath() << "rawdata/"
  << fSetName << "/ATLAS-36pb-Wmlrap.data";
  fWZ[1].open(datafile2.str().c_str(), ios::in);

  if (fWZ[1].fail()) {
    cerr << "Error opening data file " << datafile2.str() << endl;
    exit(-1);
  }

  stringstream datafile3("");
  datafile3 << dataPath() << "rawdata/"
  << fSetName << "/ATLAS-36pb-Zrap.data";
  fWZ[2].open(datafile3.str().c_str(), ios::in);

  if (fWZ[2].fail()) {
    cerr << "Error opening data file " << datafile3.str() << endl;
    exit(-1);
  }

  // Starting filter
  const double lcorr = 1.0187; // correction factor due to luminosity upgrade
  const int ndataWZ[3] = {11,11,8};  //Number of data for W+, W- and Z respectively
  const double convfac = lcorr*1000.; // Must multiply from pb to fb
  const double MWZ2[3] = {pow(MW,2.0), pow(MW,2.0), pow(MZ,2.0)};   //Mass squared of W (+ and -) and Z

  string line;
  int idat = 0;
  double etamin,etamax,tmp;

  cout << "********** WARNING: Converting pb to fb to match ApplGrid output ********" << endl;

  for (int iWZ = 0; iWZ < 3; iWZ++)
  {
    // rapidity
    getline(fWZ[iWZ],line);
    istringstream lstream(line);
    for (int i = 0; i < ndataWZ[iWZ]; i++)
    {
      lstream >> etamin >> etamax;
      fKin1[idat+i] = etamin + (etamax-etamin)*0.5;
    }

    // M_W,Z
    for (int i = 0; i < ndataWZ[iWZ]; i++)
      fKin2[idat+i] = MWZ2[iWZ];

    // sqrt(s)
    for (int i = 0; i < ndataWZ[iWZ]; i++)
      fKin3[idat+i] = 7000;

    // obs
    getline(fWZ[iWZ],line);
    istringstream lstream2(line);
    for (int i = 0; i < ndataWZ[iWZ]; i++)
      {
        lstream2 >> fData[idat+i];
        fData[idat+i] *= convfac;
      }

    // stat (%, converted later)
    getline(fWZ[iWZ],line);
    istringstream lstream3(line);
    for (int i = 0; i < ndataWZ[iWZ]; i++)
      lstream3 >> fStat[idat+i];

    // uncorrelated sys
    getline(fWZ[iWZ],line);
    istringstream lstream4(line);
    for (int i = 0; i < ndataWZ[iWZ]; i++)
    {
      lstream4 >> fSys[idat+i][0].mult;
      fSys[idat+i][0].name = "UNCORR";
    }

    // total correlated sys (unused)
    getline(fWZ[iWZ],line);

    // total uncertainty (unused)
    getline(fWZ[iWZ],line);

    // correlated systematics
    for (int isys = 2; isys < fNSys; isys++)  //2 to skip uncorr and lumi
    {
      getline(fWZ[iWZ],line);
      istringstream lstream(line);
      lstream >> tmp;
      for (int i = 0; i < ndataWZ[iWZ]; i++)
      {
        lstream >> fSys[idat+i][isys].mult;
        fSys[idat+i][isys].name = "CORR";
      }
    }

    // luminosity: 3.4%
    for (int i = 0; i < ndataWZ[iWZ]; i++)
    {
      fSys[idat+i][1].mult = 3.5;
      fSys[idat+i][1].name = "ATLASLUMI10";
    }

    idat+=ndataWZ[iWZ];
  }

  // Convert additive uncertainties to absolute form
  for (int i = 0; i < fNData; i++)
  {
    fStat[i] *= fData[i]*1e-2;
    for(int l = 0; l < fNSys; l++)
    {
      fSys[i][l].type = MULT; // All systematics multiplicative
      fSys[i][l].add = fSys[i][l].mult*fData[i]*1e-2;
    }
  }


  fWZ[0].close();
  fWZ[1].close();
  fWZ[2].close();
}


/**
 * atlas-highmass-dy
 *
 * data taken from: http://hepdata.cedar.ac.uk/view/ins1234228
 * and from paper: http://arxiv.org/abs/ARXIV:1305.4192
 *
 */
void ATLASZHIGHMASS49FBFilter::ReadData()
{
  cout << "WARNING: kinematics are not implemented" << endl;
  // Opening files
  fstream f1, f2;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/ATLAS-49fb-Zhighmass.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream datafile2("");
  datafile2 << dataPath() << "rawdata/"
  << fSetName << "/ATLAS-49fb-Zhighmass.sys";
  f2.open(datafile2.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double mbin[fNData+1], stat;

  // Reading data
  string tmp;
  getline(f1,tmp);
  for (int i = 0; i < 6; i++) getline(f2,tmp);

  // Filtering data
  for (int i = 0; i < fNData; i++)
  {
    f1 >> mbin[i] >> mbin[i+1] >> fData[i] >> tmp >> tmp >> tmp;
    fData[i] *= 1E3; // converting pb to fb

    f2 >> tmp >> tmp >> stat;
    for (int j = 0; j < fNSys; j++) f2 >> fSys[i][j].mult;

    fStat[i] = stat*fData[i]*1e-2;
    fKin2[i] = pow( 0.5*(mbin[i] + mbin[i+1]) , 2.0);
    fKin1[i] = 0.0;

    fKin3[i] = 7E3;

    for (int j = 0; j < fNSys-1; j++)
    {
      fSys[i][j].add = fSys[i][j].mult*fData[i]*1e-2;
      fSys[i][j].type = MULT;
      if (j < 2)
        fSys[i][j].name = "UNCORR"; // for the uncorrelated
      else
        fSys[i][j].name = "CORR"; // for the correlated
    }

    fSys[i][fNSys-1].add = fSys[i][fNSys-1].mult*fData[i]*1e-2;
    fSys[i][fNSys-1].type = MULT;
    fSys[i][fNSys-1].name = "ATLASLUMI11";

  }

  f1.close();
  f2.close();
}

/**
 * ATLASWPT31PB data
 *
 * data taken from: http://hepdata.cedar.ac.uk/view/ins941555
 * and from paper: http://arxiv.org/abs/ARXIV:1108.6308v2
 *
 */
void ATLASWPT31PBFilter::ReadData()
{
  // Opening files
  fstream f1, f2;
  int idum,jdum;

  stringstream datafile("");
  datafile << dataPath() << "rawdata/"
  << fSetName << "/ATLASWPT31PB.data";
  f1.open(datafile.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  stringstream datafile2("");
  datafile2 << dataPath() << "rawdata/"
  << fSetName << "/ATLASWPT31PB.covmat";
  f2.open(datafile2.str().c_str(), ios::in);

  if (f1.fail()) {
    cerr << "Error opening data file " << datafile.str() << endl;
    exit(-1);
  }

  // Starting filter
  double mbin[fNData+1], systot;

  // Reading data
  string tmp;
  getline(f1,tmp);

  // Filtering data
  for (int i = 0; i < fNData; i++)
  {

    // Data are normalized to the total cross section, Units are GeV-1
    f1 >> fKin1[i] >> mbin[i] >> mbin[i+1] >> fData[i] >> fStat[i] >> systot;
    fKin2[i] = pow(MZ,2.0);
    fKin3[i] = 7E3;

    // Question: statistical uncertainties are also given in the covariance
    // matrix. Shall I ignore them and just use the total stat uncertainty given
    // in tha above file for each bin or shall I ignore the above and use the
    // one given in the cov matrix file?

  }

  // Reading covmat
  string line;
  double statmat[fNData][fNData];
  double sysmat[fNData][fNData];
  double** covmat = new double*[fNData];
  for(int i = 0; i < fNData; i++)
  {
    covmat[i] = new double[fNData];
    for(int j = i; j < fNData; j++)
    {
      getline(f2,line);
      istringstream lstream(line);
      lstream >> idum >> jdum >> idum >> jdum >> statmat[i][j] >> sysmat[i][j];
      covmat[i][j] = sysmat[i][j] + statmat[i][j]; // Here I consistently do not add the Stat uncertainty
    }                                // Unless I do not consider the one above. Not sure!
  }

  // Make it symmetric
  for(int i = 0; i < fNData; i++)
  for(int j = 0; j < i; j++)
    covmat[i][j] = covmat[j][i];

  // Generating artificial systematics
  double** syscor = new double*[fNData];
  for(int i = 0; i < fNData; i++)
    syscor[i] = new double[fNData];

  if(!genArtSys(fNData,covmat,syscor))
   {
     cerr << " in " << fSetName << endl;
     exit(-1);
   }

  for (int i = 0; i < fNData; i++)
    for (int l = 0; l < fNSys; l++)
    {
      fSys[i][l].add = syscor[i][l];
      fSys[i][l].mult = fSys[i][l].add*100/fData[i];
      fSys[i][l].type = MULT;
      fSys[i][l].name = "CORR";
    }

  f1.close();
  f2.close();
}
