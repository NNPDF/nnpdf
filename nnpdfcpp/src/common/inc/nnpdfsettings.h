// $Id: nnpdfsettings.h 2478 2015-02-03 13:23:12Z s1044006 $
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

/**
 * \class NNPDFSettings Read/write configuration files
 * \brief Reads the .ini file that contains all the configuration
 */

#pragma once

#include "common.h"
#include <string>
#include <vector>
#include <cmath>
#include <map>

#include <NNPDF/randomgenerator.h>
#include <NNPDF/fastkernel.h>
#include <NNPDF/utils.h>
#include <NNPDF/commondata.h>
#include <NNPDF/nnpdfdb.h>
#include <NNPDF/pathlib.h>
#include <APFEL/APFELdev.h>

#include <md5.h>
#include <yaml-cpp/yaml.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>

using std::fstream;
using std::ifstream;
using std::vector;
using std::map;
using std::make_pair;

using namespace NNPDF;

class PDFBasis;

/*
 *  DataSetInfo
 *  Container struct for Dataset level info - corresponds to each dataset line in config file
 */
struct DataSetInfo
{
  const string tSetName;
  const string tSysOpt;
  const real tTrainingFraction;
  const std::vector<string> tCFactors;
  const double weight;
};

/**
 * @brief Mutation property container for each flavor
 */
struct FlMutProperty
{
  vector<int> mutsize;
  vector<real> mutprob;
};

/*
 *  PosSetInfo
 *  Container struct for Positivity level info - corresponds to each positivity line in config file
 */
struct PosSetInfo
{
  const string tSetName;  //!< Set Name
  const real   tLambda; //!< Lagrange multiplier
};

class NNPDFSettings
{
private:
  string fFileName;
  string fPDFName;
  string fResultsDir;
  string fTheoryDir;

  vector<string> fExpName;                //!< Contains the experiment names
  vector<string> fPosName;                //!< Contains the positivity names
  vector<string> fSetName;                //!< Contains the dataset names
  vector< vector<string> > fExpSetName;   //!< Contains dataset names per experiment
  vector<FlMutProperty>    fFlMutProperty;//!< Contains the mutation - not really need by improves NGA performance
  vector<int>              fArch;         //!< Contains the NN architecture

  map<string,DataSetInfo>     fDataSetInfo;  //!< Contains the dataset info
  map<string, PosSetInfo>     fPosSetInfo;   //!< Map of PosSetInfo structs


  gsl_error_handler_t * fGSL_old_handler; //!< GSL error handler
  gsl_integration_workspace * fGSLWork;   //!< GSL integration workspace

  YAML::Node fConfig;   //!< main config file
  YAML::Node fPlotting; //!< plotting config file

  map<string,string> fTheory;

  bool fThUncertainties; //!< true if the fit uses an external runcard, false otherwise.
  bool fThCovSampling;   // true if the theory covariance matrix is included in the replicas generation
  bool fThCovFitting;    // true if the theroy covariance matrix is included in the fitting chi2

public:

  NNPDFSettings(const string& folder); //!< The constructor
  ~NNPDFSettings(); //!< The destructor.

  // extra set methods
  void SetPlotFile(string const&);

  // Get methods
  YAML::Node Get(const string& item) const;
  YAML::Node Get(const string& node, const string& item) const;
  YAML::Node GetPlotting(const string& item) const;
  YAML::Node GetFile() const { return fConfig; }
  bool Exists(const string& item) const;
  bool Exists(const string& node, const string& item) const;
  string const& GetTheory(const string& item) const { return fTheory.at(item); }
  string const& GetResultsDirectory() const { return fResultsDir; }
  string const& GetTheoryDirectory()  const { return fTheoryDir;  }

  int GetNExp() const { return (int) fExpName.size(); }
  int GetNSet() const { return (int) fSetName.size(); }
  int GetNPos() const { return (int) fPosName.size(); }
  int GetNFL()  const;
  string const& GetExpName(int i) const { return fExpName[i]; }
  string const& GetSetName(int i) const { return fSetName[i]; }
  string const& GetPosName(int i) const { return fPosName[i]; }
  string const& GetPDFName() const { return fPDFName; }
  vector<string> const& GetExpSets(int i) const { return fExpSetName[i]; }
  FlMutProperty  const& GetFlMutProp(int i) const { return fFlMutProperty[i]; }
  vector<int> const& GetArch() const { return fArch; }
  map<string,string> const& GetTheoryMap() const { return fTheory; }
  bool IsQED() const;
  bool IsIC()  const;
  bool IsThUncertainties() const { return fThUncertainties; }
  bool IsThCovSampling()   const { return fThCovSampling; }
  bool IsThCovFitting()    const { return fThCovFitting; }

  gsl_integration_workspace *GetGSLWorkspace() const { return fGSLWork; } //!< GSL integration workspace

  // Check methods
  bool CheckParam(string const& param, double const& p1, double const& p2) const;  //!< Check an individual parameter
  bool CheckParam(string const& param, string const& p1, string const& p2) const;  //!< Check an individual parameter
  void VerifyConfiguration() const;    //!< Checks the log hash against filter
  void VerifyFK(FKTable* const&) const;                //!< Verify FastKernel table settings

  // Print configuration
  void PrintConfiguration(const string& filename) const;
  void PrintTheory(const string& filename) const;

  vector<string>  GetDataInfo(string const& setname, filterType useFiltered) const;
  vector<int>     GetDataMask(string const& setname, filterType useFiltered) const;
  DataSetInfo     const& GetSetInfo(string const& setname) const;
  PosSetInfo      const& GetPosInfo(string const& posname) const;

  static minType        getFitMethod(string const& method);
  static paramType      getParamType(string const& method);
  static stopType       getStopType(string const& method);
  static basisType      getFitBasisType(string const& method);
private:

  void Splash() const;
  void LoadExperiments();
  void LoadPositivities();
  void CheckBasis();
  void LoadGA();
};
