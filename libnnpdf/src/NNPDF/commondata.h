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
#include <stdexcept>
#include <map>


/** @defgroup commondata CommonData
 * \brief Classes and functions related to the handling of NNPDF standard CommonData files.
 */

/*! \addtogroup commondata
 *  @{
 */

namespace NNPDF
{
  /*! 
   *  \brief Consted information on a CommonData file - precursor to full CommonData class.
   *
   *  This struct is used to communicate the basic required information to build a CommonData instance.
   *  In this way, the general metadata is read from the file outside the CommonData constructor,
   *  and therefore can be consted in the CommonData instance. This should therefore be used when
   *  reading CommonData FILES.
   */
  struct dataInfo
  {
    /*! 
     *  dataInfo constructor (required by SWIG)
     */
    dataInfo(int nD, int nS, std::string sN, std::string tF, std::string sF):
    nData(nD),
    nSys(nS),
    SetName(sN),
    targetFile(tF),
    systypeFile(sF)
    {};

    const int nData;                //!< Number of datapoints in prototype CommonData
    const int nSys;                 //!< Number of systematic uncertainties in prototype CommonData

    const std::string SetName;      //!< Prototype CommonData set name

    const std::string targetFile;   //!< Target filename for prototype CommonData 
    const std::string systypeFile;  //!< Target filename for prototype SYSTYPE
  };

  /*! 
   *  \brief Minimal Consted information on CommonData - precursor to full CommonData class.
   *
   *  This struct is used to communicate the basic required information to build a CommonData instance,
   *  in the event where an input file is not used (for example, in buildmaster when constructing an empty
   *  commondata)
   */
  struct dataInfoRaw
  {
    /*! 
     *  dataInfoRaw constructor (required by SWIG)
     */
    dataInfoRaw(int nD, int nS, std::string sN, std::string pT):
    nData(nD),
    nSys(nS),
    SetName(sN),
    ProcType(pT)
    {};

    const int nData;            //!< Number of datapoints in prototype CommonData
    const int nSys;             //!< Number of systematic uncertainties in prototype CommonData
    const std::string SetName;  //!< Prototype CommonData set name
    const std::string ProcType; //!< Process type for prototype CommonData
  };

  /*! Specifies which type of systematic error */
  enum sysType {  
                  ADD, //!< Additive systematic
                  MULT //!< Multiplicative systematic
               };

   /*! 
   *  \brief Systematics struct - contains information on a single source of systematic error
   *
   *  This struct is used to detail how individual sources of systematic error should be treated,
   *  specifying whether they are additive or multiplicative, or whether their treatment should be
   *  randomised. Furthermore the actual value of the uncertainty in the two cases is also held.
   */
  struct sysError
  {
    double add;       //!< Value of the systematic error when additive
    double mult;      //!< Value of the systematic error when multiplicative
    sysType type;     //!< Type of the systematic error (ADD/MULT)
    std::string name; //!< Name of the systematic error (for correlation purposes)
    bool isRAND;      //!< Tag specifying whether the treatment should be randomised
    sysError() : add(0.0), mult(0.0), type(ADD), name("CORR"), isRAND(false)
    {
    };
  };

  /*! 
   *  \brief Generate a dataInfo struct from file
   *
   *  This function returns a new dataInfo struct from the provided filenames.
   *  The return value is used to generate a full CommonData instance.
   */
  dataInfo genInfoStruct(std::string const& targetfile, std::string const& sysfile);
  
  /*! 
   *  \brief Extract systematics suffix number
   *
   *  This function returns the systematic ID present in the filename provided.
   */
  int extractSysID(std::string const& sysfile);
    
  /**
   *  \class CommonData
   *  \brief Class to handle NNPDF CommonData file format
   */
  class CommonData {

  private:
    CommonData();                             //!< Disable default constructor
    CommonData& operator=(const CommonData&); //!< Disable copy-assignment

    static void VerifyProc(std::string const& proc); //!< Verify process types

  protected:
    CommonData(dataInfo const&);         //!< Constructor from file 
    CommonData(dataInfoRaw const&);      //!< Constructor for empty CommonData

    const std::string fSetName; //!< Dataset name
    const int fNData; //!< Number of datapoints
    double *fData;    //!< Data value array

    // Kinematical variables
    std::string *fProc; //!< Process (determines kinematics)
    double *fKin1;      //!< First kinematic variable array (x/pT)
    double *fKin2;      //!< Second kinematic variable array (Q/y)
    double *fKin3;      //!< Third kinematic variable array (inelasticity)

    // Uncertainties
    const int fNSys;    //!< Number of systematic errors
    const int fSysId;   //!< SYSTYPE ID
    double *fStat;      //!< Statistical error array
    sysError **fSys;    //!< Systematic errors

    /*! 
     *  Read data from provided CommonData filenames
     */
    void ReadData(std::string const& targetfile, std::string const& sysfile);

  public:
    // ******************************* CommonData Public Constructors *************************************

    CommonData(const CommonData& set);  //!< copy constructor
    CommonData(const CommonData& set, std::vector<int> const& mask); //!< Masked copy constructor

    virtual ~CommonData();	                         //!< The destructor.

    // ******************************* Process types ********************************************

    typedef std::map<std::string, std::vector<std::string>> kinMap;
    static const kinMap kinLabel_latex;

    // ******************************* CommonData Verbosity *************************************
    static bool Verbose;
    // ******************************* CommonData Get Methods *************************************

    std::string const& GetSetName() const {return fSetName; }; //!< Returns set name

    int const&  GetNData() const { return fNData;}; //!< Returns number of data points
    int const&  GetNSys()  const { return fNSys;};  //!< Returns number of systematic uncertainties

    double const& GetData(int i) const { return fData[i];} //!< Return data value for point i
    double const& GetStat(int i) const { return fStat[i];} //!< Return statistical uncertanty for point i

    sysError const& GetSys(int i, int l) const { return fSys[i][l];} //!< Return lth systematic for point i

    double* GetData() const { return fData; }; //!< Return whole data array

    // Process/Kinematics
    std::string const& GetProc(int i)     const { return fProc[i]; }; //!< Return the process type for point i
    double const& GetKinematics(int const& idat, int const& ikin) const //!< Return the ikinth kinematic value for point idat
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
          throw std::out_of_range("CommonData::GetKinematics: Kinematical variable out of range");
      }
    };

    // ********************************* CommonData File IO *****************************************

    static CommonData ReadFile(std::string const& filename, std::string const& sysfile); //!< Returns a new CommonData read from file
    void Export(std::string const& targetdir) const;  //!< Writes the current CommonData instance to file
  };
}
 /*! @} */