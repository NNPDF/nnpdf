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

#include "NNPDF/common.h"
#include "NNPDF/commondata.h"
#include "NNPDF/randomgenerator.h"

namespace NNPDF
{

  // Kinematics type labels
  // Also defines all permissable process types
  const CommonData::kinMap CommonData::kinLabel_latex = {
    { "DIS",        {"$x$","$Q^2 (GeV^2)$","$y$"}},
    { "DYP",        {"$y$","$M^2 (GeV^2)$","$\\sqrt{s} (GeV)$"}},
    { "JET",        {"$\\eta$","$p_T^2 (GeV^2)$","$\\sqrt{s} (GeV)$"}},
    { "PHT",        {"$\\eta_\\gamma$","$E_{T,\\gamma}^2 (GeV^2)$","$\\sqrt{s} (GeV)$"}},
    { "INC",        {"$0$","$\\mu^2 (GeV^2)$","$\\sqrt{s} (GeV)$"}},
    { "PDF",        {"x","$Q^2 (GeV^2)$","flavour (PID)$"}},
    { "EWK_RAP",    {"$\\eta/y$","$M^2 (GeV^2)$","$\\sqrt{s} (GeV)$"}},
    { "EWK_PT",     {"$p_T$ (GeV)","$M^2 (GeV^2)$","$\\sqrt{s} (GeV)$"}},
    { "EWK_PTRAP",  {"$\\eta/y$","$p_T^2 (GeV^2)$","$\\sqrt{s} (GeV)$"}},
    { "EWK_MLL",    {"$M_{ll} (GeV)$","$M_{ll}^2 (GeV^2)$","$\\sqrt{s} (GeV)$"}},
    { "EWJ_RAP",    {"$\\eta/y$","$M^2 (GeV^2)$","$\\sqrt{s} (GeV)$"}},
    { "EWJ_PT",     {"$p_T (GeV)$","$M^2 (GeV^2)$","$\\sqrt{s} (GeV)$"}},
    { "EWJ_PTRAP",  {"$\\eta/y$","$p_T^2 (GeV^2)$","$\\sqrt{s} (GeV)$"}},
    { "EWJ_JRAP",   {"$\\eta/y$","$M^2 (GeV^2)$","$\\sqrt{s} (GeV)$"}},
    { "EWJ_JPT",    {"$p_T (GeV)$","$M^2 (GeV^2)$","$\\sqrt{s} (GeV)$"}},
    { "EWJ_MLL",    {"$M_{ll} (GeV)$","$M_{ll}^2 (GeV^2)$","$\\sqrt{s} (GeV)$"}},
    { "HQP_YQQ",    {"$y^{QQ}$","$\\mu^2 (GeV^2)$","$\\sqrt{s} (GeV)$"}},
    { "HQP_MQQ",    {"$M^{QQ} (GeV)$","$\\mu^2 (GeV^2)$","$\\sqrt{s} (GeV)$"}},
    { "HQP_PTQQ",   {"$p_T^{QQ} (GeV)$","$\\mu^2 (GeV^2)$","$\\sqrt{s} (GeV)$"}},
    { "HQP_YQ",     {"$y^Q$","$\\mu^2 (GeV^2)$","$\\sqrt{s} (GeV)$"}},
    { "HQP_PTQ",    {"$p_T^Q (GeV)$","$\\mu^2 (GeV^2)$","$\\sqrt{s} (GeV)$"}},
    { "HIG_RAP",    {"$y$","$M_H^2 (GeV^2)$","$\\sqrt{s} (GeV)$"}},
    { "SIA" ,       {"$z$", "$Q^2 (GeV^2)$", "$y$"}}
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

  std::string extractSysID(std::string const& sysfile)
  {
    const std::string sysid = sysfile.substr(sysfile.find_last_of("_")+1, sysfile.length()-4-sysfile.find_last_of("_")-1).c_str();
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
    get_logger() << std::endl << "-- Reading COMMONDATA for Dataset: " << fSetName<<std::endl;

    get_logger() << "nData: "<<fNData
           <<" nSys: "<<fNSys
           <<std::endl;

    // Initialise arrays
    for (int i=0; i<fNData; i++)
    {
      fData[i] = std::numeric_limits<double>::quiet_NaN();
      fKin1[i] = std::numeric_limits<double>::quiet_NaN();
      fKin2[i] = std::numeric_limits<double>::quiet_NaN();
      fKin3[i] = std::numeric_limits<double>::quiet_NaN();
      fStat[i] = std::numeric_limits<double>::quiet_NaN();
      fProc[i] = "";
      fSys[i] = new sysError[fNSys];
    }

    // Read commondata
    ReadData(info.targetFile, info.systypeFile);
    Verify();
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
  fSysId("DEFAULT"),
  fStat(new double[fNData]),
  fSys(new sysError*[fNData])
  {
    get_logger() << std::endl << "-- Reading COMMONDATA for Dataset: " << fSetName<<std::endl;

    get_logger() << "nData: "<<fNData
              <<" nSys: "<<fNSys
              <<std::endl;

    VerifyProc(info.ProcType);

    // Initialise arrays
    for (int i=0; i<fNData; i++)
    {
      fData[i] = std::numeric_limits<double>::quiet_NaN();
      fKin1[i] = std::numeric_limits<double>::quiet_NaN();
      fKin2[i] = std::numeric_limits<double>::quiet_NaN();
      fKin3[i] = std::numeric_limits<double>::quiet_NaN();
      fStat[i] = std::numeric_limits<double>::quiet_NaN();
      fProc[i] = info.ProcType;
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

    Verify();
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
  fSysId("-1"),
  fStat(nullptr),
  fSys(nullptr)
  {
      using std::swap;
      swap(*this, other);
      Verify();
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

    Verify();
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
    if (CommonData::kinLabel_latex.count(proc) == 0)
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

    Verify();
    get_logger() << "-- COMMONDATA Files for "<<fSetName<<" successfully read."<<std::endl<<std::endl;

  }

  // Verifies one field in CommonData
  void VerifyField( double field, std::string name, bool& verification )
  {
    if (std::isnan(field)){
        std::cerr << "CommonData::Verify: " + name + " unset"<<std::endl;
        verification = false;
    }
  }
  // Verify that all fields are set in CommonData
  void CommonData::Verify() const
  {
    bool pass_verification = true;
    for (int i=0; i<GetNData(); i++)
    {
        // These VerifyField calls set pass_verification to false in the case of failure
        const std::string dps = std::to_string(i);
        VerifyField( GetData(i), "Data point " + dps, pass_verification);
        VerifyField( GetStat(i), "Statistical Error " + dps, pass_verification);
        VerifyField( GetKinematics(i,0), "Kinematic1 for datapoint " + dps, pass_verification);
        VerifyField( GetKinematics(i,1), "Kinematic2 for datapoint " + dps, pass_verification);
        VerifyField( GetKinematics(i,2), "Kinematic3 for datapoint " + dps, pass_verification);
        for (int j=0; j<GetNSys(); j++)
        {
            const std::string sps = std::to_string(j);
            VerifyField( GetSys(i,j).add, "Additive systematic " + sps + " for datapoint "+dps, pass_verification);
            VerifyField( GetSys(i,j).mult,"Multiplicative systematic " + sps + " for datapoint "+dps, pass_verification);
            if (GetSys(i,j).type == UNSET){
                std::cerr << "CommonData::Verify: Systematic point type "+sps+" for datapoint "+dps+" unset"<<std::endl;
                pass_verification = false;
            }
        }
    }
    if (pass_verification == false)
        throw std::runtime_error("CommonData::Verify failed for set " + GetSetName());
  }

  // Write data to file in CommonData format
  void CommonData::Export(std::string const& targetdir) const
  {
    Verify();
    std::fstream g1;

    // output datafile
    std::string datafileout = targetdir + "/DATA_" + fSetName + ".dat";
    g1.open(datafileout.c_str(), std::ios::out);

    get_logger() << "-- Exporting "<<fSetName<<" to "<< datafileout<<std::endl;


    g1 << fSetName << "\t"
    << fNSys  << "\t"
    << fNData <<std::endl;

    int idat = 1;
    for (int i=0; i<fNData; i++)
    {
      g1.precision(12);
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

    get_logger() << "-- Exporting "<<fSetName<<" SYS to "<< sysfileout.str()<<std::endl;

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

  // Get uncorrelated error for datapoint i
  double CommonData::GetUncE( int i ) const
  {
    double uncorrErr = pow(GetStat(i),2);
    for (int l=0; l<GetNSys(); l++)
      if (GetSys(i,l).name == "UNCORR")
        uncorrErr += pow(GetSys(i,l).add,2);
    return sqrt(uncorrErr);
  }

  // Get total correlated error for datapoint i
  double CommonData::GetCorE( int i ) const
  {
    double sysErr = 0.0;
    for (int l=0; l<GetNSys(); l++)
      if (GetSys(i,l).name != "UNCORR")
        sysErr += pow(GetSys(i,l).add,2);
    return sqrt(sysErr);
  }
}
