// $Id: commondata.cc 639 2013-03-25 19:08:38Z s1044006 $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <sys/stat.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include "NNPDF/commondata.h"
#include "NNPDF/randomgenerator.h"

namespace NNPDF
{

  // Default verbosity
  bool CommonData::Verbose = true;

  /**
   * Generate a dataInfo struct given a target filename
   */
  dataInfo genInfoStruct(std::string const& targetfile, std::string const& sysfile)
  {
    // Open commondata file
    std::ifstream datafile;
    datafile.open(targetfile.c_str());

    if (!datafile.good())
    {
      std::cerr << "genInfoStruct Error: Cannot read commondata file from: "<<std::endl<<targetfile<<std::endl;
      exit(-1);
    }

    // Read metadata
    std::string SetName;
    int nSys, nData;

    int dummy;

    datafile >> SetName
    >> nSys
    >> nData;

    datafile.close();

    // Setup info struct
    const dataInfo infoStruct = {
      nData,
      nSys,
      SetName,
      targetfile,
      sysfile
    };

    return infoStruct;
  }

  int extractSysID(std::string const& sysfile)
  {
    const int sysid = atoi(sysfile.substr(sysfile.find_last_of("_")+1, sysfile.length()-4-sysfile.find_last_of("_")-1).c_str());
    return sysid;
  }

  /**
   * Helper function to parse strings into SYSTYPES
   */
  sysType parseSYS(std::string const& str)
  {
    if (str.compare("ADD") == 0) return ADD;
    else if (str.compare("MULT") == 0) return MULT;
    else
    {
      std::cerr << "parseSYS Error: Unrecognised systematic:" << str << std::endl;
      exit(-1);
    }

    return MULT;
  }

  /**
   * CommonData base class constructor
   */
  CommonData::CommonData(dataInfo const& info):
  fSetName(info.SetName),
  fNData(info.nData),
  fData(new double[fNData]),
  fProc(new std::string[fNData]),
  fKin1(new double[fNData]),
  fKin2(new double[fNData]),
  fKin3(new double[fNData]),
  fNSys(info.nSys),
  fSysId(extractSysID(info.systypeFile)),
  fStat(new double[fNData]),
  fSys(new sysError*[fNData])
  {
    if (Verbose)
    {
      std::cout << std::endl << "-- Reading COMMONDATA for Dataset: " << fSetName<<std::endl;

      std::cout << "nData: "<<fNData
           <<" nSys: "<<fNSys
           <<std::endl;
    }

    // Initialise arrays
    for (int i=0; i<fNData; i++)
    {
      fData[i] = 0.0;
      fKin1[i] = 0.0;
      fKin2[i] = 0.0;
      fKin3[i] = 0.0;

      fProc[i] = "";

      fStat[i] = 0.0;
      fSys[i] = new sysError[fNSys];
    }

    // Read commondata
    ReadData(info.targetFile, info.systypeFile);
  }

  /**
   * CommonData base class constructor
   */
  CommonData::CommonData(dataInfoRaw const& info):
  fSetName(info.SetName),
  fNData(info.nData),
  fData(new double[fNData]),
  fProc(new std::string[fNData]),
  fKin1(new double[fNData]),
  fKin2(new double[fNData]),
  fKin3(new double[fNData]),
  fNSys(info.nSys),
  fSysId(0),
  fStat(new double[fNData]),
  fSys(new sysError*[fNData])
  {
    std::cout << std::endl << "-- Reading COMMONDATA for Dataset: " << fSetName<<std::endl;

    std::cout << "nData: "<<fNData
              <<" nSys: "<<fNSys
              <<std::endl;

    // Initialise arrays
    for (int i=0; i<fNData; i++)
    {
      fData[i] = 0.0;
      fKin1[i] = 0.0;
      fKin2[i] = 0.0;
      fKin3[i] = 0.0;
      fProc[i] = info.ProcType;
      fStat[i] = 0.0;
      fSys[i] = new sysError[fNSys];
    }
  }

    /**
   * CommonData copy constructor
   */
  CommonData::CommonData(const CommonData& set):
  fSetName(set.fSetName),
  fNData(set.fNData),
  fData(new double[fNData]),
  fProc(new std::string[fNData]),
  fKin1(new double[fNData]),
  fKin2(new double[fNData]),
  fKin3(new double[fNData]),
  fNSys(set.fNSys),
  fSysId(set.fSysId),
  fStat(new double[fNData]),
  fSys(new sysError*[fNData])
  {
    // Setup masked data
    for (int i = 0; i < fNData; i++)
    {
      fData[i] = set.fData[i];
      fStat[i] = set.fStat[i];

      fSys[i] = new sysError[fNSys];

      fKin1[i] = set.fKin1[i];
      fKin2[i] = set.fKin2[i];
      fKin3[i] = set.fKin3[i];

      fProc[i] = set.fProc[i];

      for (int l = 0; l < fNSys; l++)
        fSys[i][l] = set.fSys[i][l];
    }

  }

  /**
   * CommonData masked copy constructor
   */
  CommonData::CommonData(const CommonData& set, std::vector<int> const& mask):
  fSetName(set.fSetName),
  fNData(mask.size()),
  fData(new double[fNData]),
  fProc(new std::string[fNData]),
  fKin1(new double[fNData]),
  fKin2(new double[fNData]),
  fKin3(new double[fNData]),
  fNSys(set.fNSys),
  fSysId(set.fSysId),
  fStat(new double[fNData]),
  fSys(new sysError*[fNData])
  {
    // Setup masked data
    for (int i = 0; i < fNData; i++)
    {
      fData[i] = set.fData[mask[i]];
      fStat[i] = set.fStat[mask[i]];

      fSys[i] = new sysError[fNSys];

      fKin1[i] = set.fKin1[mask[i]];
      fKin2[i] = set.fKin2[mask[i]];
      fKin3[i] = set.fKin3[mask[i]];

      fProc[i] = set.fProc[mask[i]];

      for (int l = 0; l < fNSys; l++)
        fSys[i][l] = set.fSys[mask[i]][l];
    }

  }


  /**
   * The destructors
   */
  CommonData::~CommonData()
  {
    for (int i=0; i<fNData; i++)
    {
      delete[] fSys[i];
    }

    delete[] fSys;

    delete[] fData;
    delete[] fStat;
    delete[] fKin1;
    delete[] fKin2;
    delete[] fKin3;
    delete[] fProc;

  }

  CommonData CommonData::ReadFile(std::string const& filename, std::string const& sysfile)
  {
    // Setup info struct
    const dataInfo infoStruct = genInfoStruct(filename, sysfile);
    return CommonData(infoStruct);
  }

  void CommonData::ReadData(std::string const & targetfile, std::string const& systypefilename)
  {
    // read datafile
    std::ifstream datafile;
    datafile.open(targetfile.c_str());

    if (!datafile.good())
    {
      std::cerr << "CommonData::ReadData Error - file "<<targetfile<<" is bad! "<<std::endl;
      exit(-1);
    }

    std::string setname, proc;
    int nsys, ndata;

    // Check metadata
    datafile >> setname
    >> nsys
    >> ndata;

    // Verification
    if (setname != fSetName)
    {
      std::cerr << "CommonData::ReadData Error: Setname Mismatch."<<std::endl;
      exit(-1);
    }

    if (nsys != fNSys)
    {
      std::cerr << "CommonData::ReadData Error: N_Uncertainty Mismatch"<<std::endl;
      exit(-1);
    }

    if (ndata != fNData)
    {
      std::cerr << "CommonData::ReadData Error: NData Mismatch"<<std::endl;
      exit(-1);
    }

    // Read data
    for (int i=0; i<fNData; i++)
    {
      int idat = 0;

      datafile >> idat
      >> fProc[i]
      >> fKin1[i]
      >> fKin2[i]
      >> fKin3[i]
      >> fData[i]
      >> fStat[i];

      if (idat != i+1)
      {
        std::cerr << "CommonData::ReadData Error: Datapoint Mismatch: "<<idat<<" vs "<< i<<std::endl;
        exit(-1);
      }

      // Read systematics
      for (int l = 0; l < fNSys; l++)
        datafile >> fSys[i][l].add >> fSys[i][l].mult;

    }

    datafile.close();

    // ************ Reading Systypes *******************

    std::ifstream h;
    h.open(systypefilename.c_str());

    if (h.fail()) {
      std::cerr << "Error opening systype file " << systypefilename << std::endl;
      exit(-1);
    }

    // Verify number of systematics in SYSTYPE file adds up
    h >> nsys;
    if (nsys != fNSys)
    {
      std::cerr << "DataSet::ReadData Error: Number of systematics for " << fSetName << " doesn't match" << std::endl;
      exit(-1);
    }

    // Read types and names for systematics
    for (int l = 0; l < fNSys; l++)
    {
      int sysnum; std::string tystr,sysname;
      h >> sysnum >> tystr >> sysname;

      const bool random = (tystr.compare("RAND") == 0);
      const sysType type = random ? MULT:parseSYS(tystr);

      for (int i = 0; i < fNData; i++)
      {
        fSys[i][l].name = sysname;
        fSys[i][l].isRAND = random;
        fSys[i][l].type = type;
      }
    }

    h.close();

    if (Verbose)
      std::cout << "-- COMMONDATA Files for "<<fSetName<<" successfully read."<<std::endl<<std::endl;

  }


  /**
   * Write data to file in CommonData format
   */
  void CommonData::Export(std::string const& targetdir, const double& minQ2) const
  {
    std::fstream g1;

    // Minimum Q^2 value for DIS data
    int nPassData = fNData;
    bool DIS = false;
    if (fProc[0].compare(0, 3, std::string("DIS")) == 0)
    {
      DIS = true;
      nPassData = 0;
      for (int i=0; i<fNData; i++)
        if (fKin2[i] > minQ2) nPassData++;

      if (nPassData != fNData)
        std::cout << "DIS Data min Q^2 > "<< minQ2
                  <<" Cut - "<<nPassData<<"/"<<fNData<<" data points remain"<< std::endl;
    }

    // output datafile
    std::string datafileout = targetdir + "/DATA_" + fSetName + ".dat";
    g1.open(datafileout.c_str(), std::ios::out);

    if (Verbose)
      std::cout << "-- Exporting "<<fSetName<<" to "<< datafileout<<std::endl;


    g1 << fSetName << "\t"
    << fNSys  << "\t"
    << nPassData <<std::endl;

    int idat = 1;
    for (int i=0; i<fNData; i++)
    {
      if (DIS && fKin2[i] <= minQ2) continue;

      g1.precision(15);
      g1 << std::scientific << std::setw(4) << idat
      << "\t" << fProc[i]
      << "\t" << fKin1[i]
      << "\t" << fKin2[i]
      << "\t" << fKin3[i]
      << "\t" << fData[i]
      << "\t" << fStat[i] << "\t";

      for (int l = 0; l < fNSys; l++)
        g1 << fSys[i][l].add << "\t" << fSys[i][l].mult << "\t";
      g1 << std::endl;

      idat++;
    }

    g1.close();

    // Export SYSTYPES
    const std::string sysNames[2] = {"ADD","MULT"};

    std::stringstream sysfileout;
    sysfileout << targetdir << "/systypes";

    // Make SYSTYPES folder
    mkdir(sysfileout.str().c_str(), 0777);
    sysfileout << "/SYSTYPE_" << fSetName <<"_" << fSysId << ".dat";

    std::ofstream h;
    h.open(sysfileout.str().c_str());

    if (h.fail()) {
      std::cerr << "Error writing systype file " << sysfileout.str() << std::endl;
      exit(-1);
    }

    if (Verbose)
      std::cout << "-- Exporting "<<fSetName<<" SYS to "<< sysfileout.str()<<std::endl;

    // Total number of systematics
    h << fNSys << std::endl;

    // Systematics
    for (int l = 0; l < fNSys; l++)
    {
      h << l+1 <<"    ";
      if (fSys[0][l].isRAND) h << "RAND";
      else h << sysNames[fSys[0][l].type];
      h << "    " << fSys[0][l].name << std::endl;
    }

    h.close();

  }
}
