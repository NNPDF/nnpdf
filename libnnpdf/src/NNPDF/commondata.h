// $Id: commondata.h 3004 2015-06-11 16:13:23Z s1044006 $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "common.h"

#include <limits>
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
                  ADD,  //!< Additive systematic
                  MULT, //!< Multiplicative systematic
                  UNSET //!< Unset systematic
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

    // Normal constructor - sets quiet_NaNs as sentinel values
    sysError():
    add(std::numeric_limits<double>::quiet_NaN()),
    mult(std::numeric_limits<double>::quiet_NaN()),
    type(UNSET),
    name("CORR"),
    isRAND(false)
    {};

    // Copy constructor
    sysError(sysError const& o):
    add(o.add),
    mult(o.mult),
    type(o.type),
    name(o.name),
    isRAND(o.isRAND)
    {};
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
  std::string extractSysID(std::string const& sysfile);

  /**
   *  \class CommonData
   *  \brief Class to handle NNPDF CommonData file format
   */
  class CommonData {

  private:
    CommonData();                             //!< Disable default constructor

    static void VerifyProc(std::string const& proc); //!< Verify process types

  protected:
    CommonData(dataInfo const&);         //!< Constructor from file
    CommonData(dataInfoRaw const&);      //!< Constructor for empty CommonData

    std::string fSetName; //!< Dataset name
    int fNData; //!< Number of datapoints
    double *fData;    //!< Data value array

    // Kinematical variables
    std::string *fProc; //!< Process (determines kinematics)
    double *fKin1;      //!< First kinematic variable array (x/pT)
    double *fKin2;      //!< Second kinematic variable array (Q/y)
    double *fKin3;      //!< Third kinematic variable array (inelasticity)

    // Uncertainties
    int fNSys;          //!< Number of systematic errors
    std::string fSysId; //!< SYSTYPE ID
    double *fStat;      //!< Statistical error array
    sysError **fSys;    //!< Systematic errors

    /*!
     *  Read data from provided CommonData filenames
     */
    void ReadData(std::string const& targetfile, std::string const& sysfile);

  public:
    // ******************************* CommonData Public Constructors *************************************

    CommonData(const CommonData& set);  //!< copy constructor
    CommonData(CommonData&& set);  //!< Move constructor
    CommonData& operator=(CommonData); //!< Copy-assignment
    CommonData& operator=(CommonData&&); //!< Move-assignment
    CommonData(const CommonData& set, std::vector<int> const& mask); //!< Masked copy constructor
    friend void swap(CommonData&, CommonData& );

    virtual ~CommonData();	                         //!< The destructor.

    // ******************************* Process types ********************************************

    typedef std::map<std::string, std::vector<std::string> > kinMap;
    static const kinMap kinLabel_latex;

    // ******************************* CommonData Get Methods *************************************

    std::string const& GetSetName() const {return fSetName; }; //!< Returns set name

    int const&  GetNData() const { return fNData;}; //!< Returns number of data points
    int const&  GetNSys()  const { return fNSys;};  //!< Returns number of systematic uncertainties

    double const& GetData(int i) const { return fData[i];} //!< Return data value for point i
    double const& GetStat(int i) const { return fStat[i];} //!< Return statistical uncertanty for point i
    double        GetUncE(int i) const; //!< Return total uncorrelated error for point i
    double        GetCorE(int i) const; //!< Return total correlated error for point i

    sysError const& GetSys(int i, int l) const { return fSys[i][l];} //!< Return lth systematic for point i
    sysError** GetSysErrors() const { return fSys; } //!< Return full systematic matrix

    double* GetData() const { return fData; }; //!< Return whole data array

    // Process/Kinematics
    std::string const& GetProc(int i)     const { return fProc[i]; }; //!< Return the process type for point i
    double const& GetKinematics(int const& idat, int const& ikin) const //!< Return the ikinth kinematic value for point idat
    {
	  if (idat < 0 || idat >= fNData)
	  {
		throw std::out_of_range("CommonData::GetKinematics: Data point out of range");
	  }

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
    void Verify() const; //!< Verifies the current CommonData, checking that all fields are set
    void Export(std::string const& targetdir) const;  //!< Writes the current CommonData instance to file
  };

}
 /*! @} */
