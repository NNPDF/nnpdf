// $Id: commondata.h 3004 2015-06-11 16:13:23Z s1044006 $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "common.h"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <cstdlib>

namespace NNPDF
{
  // Data Info Struct
  struct dataInfo
  {
    const int nData;
    const int nSys;

    const std::string SetName;

    const std::string targetFile;
    const std::string systypeFile;
  };

  // Raw Data Info for Export
  struct dataInfoRaw
  {
    const int nData;
    const int nSys;
    const std::string SetName;
    const std::string ProcType;
  };

  // Systematics Struct
  enum sysType { ADD, MULT };
  struct sysError
  {
    double add;
    double mult;
    sysType type;
    std::string name;
    bool isRAND;
    sysError() : add(0.0), mult(0.0), type(ADD), name("CORR"), isRAND(false)
    {
    };
  };

  // Generate a dataInfo struct from file
  dataInfo genInfoStruct(std::string const& targetfile, std::string const& sysfile);
  
  // Extract systematics suffix number
  int extractSysID(std::string const& sysfile);
    
  /**
   *  \class CommonData
   *  \brief Class to handle data from experiments
   */

  class CommonData {

  private:
    CommonData();                             //disable default constructor
    CommonData& operator=(const CommonData&); //disable copy-assignment

  protected:
    CommonData(dataInfo const&);         //!< The constructor.
    CommonData(dataInfoRaw const&);

    // Set name
    const std::string fSetName;

    // Data points
    const int fNData;
    double *fData;

    // Kinematical variables
    std::string *fProc; // Process (determined kinematics)
    double *fKin1;      // x/pT
    double *fKin2;      // Scale/rapidity
    double *fKin3;      // inelasticity

    //uncertainties
    const int fNSys;
    const int fSysId;
    double *fStat;
    sysError **fSys;

    // Read data from file
    void ReadData(std::string const& targetfile, std::string const& sysfile);

  public:
    // ******************************* CommonData Public Constructors *************************************

    CommonData(const CommonData& set);  //!< copy constructor
    CommonData(const CommonData& set, std::vector<int> const& mask); //!< Masked copy constructor

    virtual ~CommonData();	                         //!< The destructor.

    // ******************************* CommonData Verbosity *************************************

    static bool Verbose;

    // ******************************* CommonData Get Methods *************************************

    std::string const& GetSetName()  const {return fSetName; };

    int       const&  GetNData()                    const { return fNData;};          //!< Returns number of data points
    int       const&  GetNSys()                     const { return fNSys;};           //!< Returns number of uncorrelated systematics

    double      const&  GetData(int i)                const { return fData[i];}          //!< Return data value
    double      const&  GetStat(int i)                const { return fStat[i];}          //!< Return statistical uncertanty

    sysError  const&  GetSys(int i, int l)          const { return fSys[i][l];}        //!< Return systematic

    double*     const   GetData()                     const { return fData; };          //!< Return data array

    // Process/Kinematics
    std::string const& GetProc(int i)     const { return fProc[i]; };
    double const& GetKinematics(int const& idat, int const& ikin) const
    {
      switch (ikin)
      {
        case 0:
          return fKin1[idat];
          break;

        case 1:
          return fKin2[idat];
          break;

        case 2:
          return fKin3[idat];
          break;

        default:
          std::cerr << "CommonData::GetKinematics Error: No such kinematical variable: "<<ikin<<std::endl;
          exit(-1);
          break;
      }
    };

    // ********************************* CommonData File IO *****************************************

    static CommonData ReadFile(std::string const& filename, std::string const& sysfile);
    void Export(std::string const& targetdir, double const& minQ2 = 0.0 ) const;
  };
}
