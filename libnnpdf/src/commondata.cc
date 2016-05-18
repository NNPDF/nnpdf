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
#include <utility>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include "NNPDF/commondata.h"
#include "NNPDF/randomgenerator.h"

namespace NNPDF
{
  // Default verbosity
  bool CommonData::Verbose = true;

  // Kinematics type labels
  const CommonData::kinMap CommonData::kinLabel_latex = { 
    { "DIS",        {"$x$","$Q^2$","$y$"}},
    { "DYP",        {"$y$","$M^2$","$\\sqrt{s}$"}},
    { "JET",        {"$\\eta$","$p_T^2$","$\\sqrt{s}$"}},
    { "PHT",        {"$\\eta_\\gamma$","$E_{T,\\gamma}^2$","$\\sqrt{s}$"}},
    { "INC",        {"$0$","$\\mu^2$","$\\sqrt{s}$"}},
    { "EWK_RAP",    {"$\\eta/y$","$M^2$","$\\sqrt{s}$"}},
    { "EWK_PT",     {"$p_T$","$M^2$","$\\sqrt{s}$"}},
    { "EWK_PTRAP",  {"$\\eta/y$","$p_T^2$","$\\sqrt{s}$"}},
    { "EWK_MLL",    {"$M_{ll}$","$M_{ll}^2$","$\\sqrt{s}$"}},
    { "EWJ_RAP",    {"$\\eta/y$","$M^2$","$\\sqrt{s}$"}},
    { "EWJ_PT",     {"$p_T$","$M^2$","$\\sqrt{s}$"}},
    { "EWJ_PTRAP",  {"$\\eta/y$","$p_T^2$","$\\sqrt{s}$"}},
    { "EWJ_JRAP",   {"$\\eta/y$","$M^2$","$\\sqrt{s}$"}},
    { "EWJ_JPT",    {"$p_T$","$M^2$","$\\sqrt{s}$"}},
    { "EWJ_MLL",    {"$M_{ll}$","$M_{ll}^2$","$\\sqrt{s}$"}},
    { "HQP_YQQ",    {"$y^{QQ}$","$\\mu^2$","$\\sqrt{s}$"}},
    { "HQP_MQQ",    {"$M^{QQ}$","$\\mu^2$","$\\sqrt{s}$"}},
    { "HQP_PTQQ",   {"$p_T^{QQ}$","$\\mu^2$","$\\sqrt{s}$"}},
    { "HQP_YQ",     {"$y^Q$","$\\mu^2$","$\\sqrt{s}$"}},
    { "HQP_PTQ",    {"$p_T^Q$","$\\mu^2$","$\\sqrt{s}$"}},
    { "HIG_RAP",    {"$y$","$M_H^2$","$\\sqrt{s}$"}},
    { "SIA" ,       {"$z$", "$Q^2$", "$y$"}}
  };

  // Generate a dataInfo struct given a target filename
  dataInfo genInfoStruct(std::string const& targetfile, std::string const& sysfile)
  {
    // Open commondata file
    std::ifstream datafile;
    datafile.open(targetfile.c_str());

    if (!datafile.good())
      throw std::invalid_argument("genInfoStruct: Cannot read commondata file from: " + targetfile);

    // Read metadata
    std::string SetName;
    int nSys, nData;

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

  // Helper function to parse strings into SYSTYPES
  sysType parseSYS(std::string const& str)
  {
    if (str.compare("ADD") == 0) return ADD;
    else if (str.compare("MULT") == 0) return MULT;
    else
      throw std::invalid_argument("parseSYS: Unrecognised systematic: " + str);

    return MULT;
  }

  // CommonData base class constructor
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

  // CommonData base class constructor
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

    VerifyProc(info.ProcType);

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

  // CommonData copy constructor
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
      VerifyProc(fProc[i]);

      for (int l = 0; l < fNSys; l++)
        fSys[i][l] = set.fSys[i][l];
    }

  }


  void swap(CommonData & lhs, CommonData & rhs) {
    using std::swap;
    swap(lhs.fSetName, rhs.fSetName);
    swap(lhs.fNData, rhs.fNData);
    swap(lhs.fData, rhs.fData);
    swap(lhs.fProc, rhs.fProc);
    swap(lhs.fKin1, rhs.fKin1);
    swap(lhs.fKin2, rhs.fKin2);
    swap(lhs.fKin3, rhs.fKin3);
    swap(lhs.fNSys, rhs.fNSys);
    swap(lhs.fSysId, rhs.fSysId);
    swap(lhs.fStat, rhs.fStat);
    swap(lhs.fSys, rhs.fSys);
  }

  CommonData & CommonData::operator =(CommonData other){
    using std::swap;
    swap(*this, other);
    return *this;
  }


  CommonData::CommonData(CommonData&& other):
  fSetName(std::string()),
  fNData(0), //Set this to maintain the class invariant
  fData(nullptr),
  fProc(nullptr),
  fKin1(nullptr),
  fKin2(nullptr),
  fKin3(nullptr),
  fNSys(0),
  fSysId(-1),
  fStat(nullptr),
  fSys(nullptr)
  {
      using std::swap;
      swap(*this, other);

  }

  // CommonData masked copy constructor
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
      VerifyProc(fProc[i]);

      for (int l = 0; l < fNSys; l++)
        fSys[i][l] = set.fSys[mask[i]][l];
    }

  }

  // Destructor
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

  // Verify that the process type is one of allowed processes
  void CommonData::VerifyProc(std::string const& proc)
  {
    const int nProc = 17;
    const std::string validProc[nProc] = {
      "DIS",
      "DYP",
      "JET",
      "PHT",
      "INC",
      "EWK_RAP",
      "EWK_PT",
      "EWK_MLL",
      "EWJ_RAP",
      "EWJ_PT",
      "EWJ_MLL",
      "HQP_YQQ",
      "HQP_MQQ",
      "HQP_PTQQ",
      "HQP_YQ",
      "HQP_PTQ",
      "SIA"
    }; 

    bool foundString = false;
    for (int i=0; i<nProc; i++)
      foundString = foundString || (proc.find(validProc[i]) == 0);

    if (!foundString)
      throw std::invalid_argument("CommonData::VerifyProc: process " + proc + " is unsupported.");

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
      throw std::invalid_argument("CommonData::ReadData: file " + targetfile + " is bad!");

    std::string setname, proc;
    int nsys, ndata;

    // Check metadata
    datafile >> setname
    >> nsys
    >> ndata;

    // Verification
    if (setname != fSetName)
      throw std::runtime_error("CommonData::ReadData: Setname Mismatch.");

    if (nsys != fNSys)
      throw std::runtime_error("CommonData::ReadData: N_Uncertainty Mismatch.");
 
    if (ndata != fNData)
      throw std::runtime_error("CommonData::ReadData: NData Mismatch.");

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

      VerifyProc(fProc[i]);

      if (idat != i+1)
        throw std::runtime_error("CommonData::ReadData: Datapoint Mismatch.");

      // Read systematics
      for (int l = 0; l < fNSys; l++)
        datafile >> fSys[i][l].add >> fSys[i][l].mult;

    }

    datafile.close();

    // ************ Reading Systypes *******************

    std::ifstream h;
    h.open(systypefilename.c_str());

    if (h.fail())
      throw std::runtime_error("CommonData::ReadData: Error opening systype file "+systypefilename);

    // Verify number of systematics in SYSTYPE file adds up
    h >> nsys;
    if (nsys != fNSys)
      throw std::runtime_error("CommonData::ReadData: Number of systematics for " + fSetName + " doesn't match");

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


  // Write data to file in CommonData format
  void CommonData::Export(std::string const& targetdir) const
  {
    std::fstream g1;

    // output datafile
    std::string datafileout = targetdir + "/DATA_" + fSetName + ".dat";
    g1.open(datafileout.c_str(), std::ios::out);

    if (Verbose)
      std::cout << "-- Exporting "<<fSetName<<" to "<< datafileout<<std::endl;


    g1 << fSetName << "\t"
    << fNSys  << "\t"
    << fNData <<std::endl;

    int idat = 1;
    for (int i=0; i<fNData; i++)
    {
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

    if (h.fail())
      throw std::runtime_error("CommonData::ReadData: Error writing systype file "+sysfileout.str());

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
