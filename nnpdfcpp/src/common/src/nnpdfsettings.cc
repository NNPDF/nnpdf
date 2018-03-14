// $Id: nnpdfsettings.cc 2478 2015-02-03 13:23:12Z s1044006 $
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <sys/stat.h>

#include "nnpdfsettings.h"
#include "version.h"

// Strings for config file output
static const string minString[6]   = {"UNDEFINED", "GA", "NGA", "NGAFT","CMAES"};
static const string stopString[6]  = {"UNDEFINED", "FIXEDLENGTH", "LOOKBACK"};
static const string paramString[6] = {"UNDEFINED", "NN", "SLN", "SLNPP"};
static const string basisString[17]= {"UNDEFINED", "NN23", "NN23QED","EVOL", "EVOLQED","EVOLS",
                                      "EVOLSQED","NN30", "NN30QED","FLVR", "FLVRQED","NN30IC",
                                      "EVOLIC","NN31IC","LUX", "NN31ICQED", "NSR"};

static const vector< vector<string> > basiselem = { {},
                                     {"sng","g","v","t3","ds","sp","sm"},          //NN23
                                     {"sng","g","v","t3","ds","sp","sm","pht"},    //NN23QED
                                     {"sng","g","v","v3","v8","t3","t8"},          //EVOL
                                     {"sng","g","v","v3","v8","t3","t8","pht"},    //EVOLQED
                                     {"sng","g","v","v8","t3","t8","ds"},          //EVOLS
                                     {"sng","g","v","v8","t3","t8","ds","pht"},    //EVOLSQED
                                     {"sng","g","v","v8","t3","t8","ds"},          //NN30
                                     {"sng","g","v","v8","t3","t8","ds","pht"},    //NN30QED
                                     {"g","u","ubar","d","dbar","s","sbar"},       //FLVR
                                     {"g","u","ubar","d","dbar","s","sbar","pht"}, //FLVRQED
                                     {"sng","g","v","t3","ds","sp","sm","cp"},     //NN30IC
                                     {"sng","g","v","v3","v8","t3","t8","t15"},    //EVOLIC
                                     {"sng","g","v","v3","v8","t3","t8","cp"},     //NN31IC
                                     {"sng","g","v","v3","v8","t3","t8","cp"},     //LUX
                                     {"sng","g","v","v3","v8","t3","t8","cp","pht"}//NN31ICQED
                                     };

/* Convert string to enum */
minType NNPDFSettings::getFitMethod(string const& method)
{
  if (method.compare("GA") == 0)      return MIN_GA;
  if (method.compare("NGA") == 0)     return MIN_NGA;
  if (method.compare("NGAFT") == 0)     return MIN_NGAFT;
  if (method.compare("CMAES") == 0)     return MIN_CMAES;

  cerr << "getFitMethod Error: Invalid fit method: "<<method<<endl;
  cerr << "choices are: "<<endl;
  cerr <<" - GA (Basic GA)"<<endl;
  cerr <<" - NGA"<<endl;
  cerr <<" - NGAP"<<endl;
  cerr <<" - NGAFT"<<endl;
  exit(-1);

  return MIN_UNDEF;
}

paramType NNPDFSettings::getParamType(string const& method)
{
  if (method.compare("NN") == 0)        return PARAM_NN;
  if (method.compare("SLN") == 0)        return PARAM_SLN;
  if (method.compare("SLNPP") == 0)        return PARAM_SLNPP;

  cerr << "getParamType Error: Invalid parametrization type: "<<method;
  exit(-1);

  return PARAM_UNDEF;
}

stopType NNPDFSettings::getStopType(string const& method)
{
  if (method.compare("LOOKBACK") == 0)    return STOP_LB;
  if (method.compare("FIXEDLENGTH") == 0) return STOP_NONE;

  cerr << "getStopType Error: Invalid stopping type: "<<method;
  exit(-1);

  return STOP_UNDEF;
}

basisType NNPDFSettings::getFitBasisType(string const& method)
{
  if (method.compare("NN23") == 0)    return BASIS_NN23;
  if (method.compare("NN23QED") == 0) return BASIS_NN23QED;
  if (method.compare("EVOL") == 0)    return BASIS_EVOL;
  if (method.compare("EVOLQED") == 0) return BASIS_EVOLQED;
  if (method.compare("EVOLS") == 0)   return BASIS_EVOLS;
  if (method.compare("EVOLSQED") == 0)return BASIS_EVOLSQED;
  if (method.compare("NN30") == 0)    return BASIS_NN30;
  if (method.compare("NN30QED") == 0) return BASIS_NN30QED;
  if (method.compare("FLVR") == 0)    return BASIS_FLVR;
  if (method.compare("FLVRQED") == 0) return BASIS_FLVRQED;
  if (method.compare("NN30IC") == 0)  return BASIS_NN30IC;
  if (method.compare("EVOLIC") == 0)  return BASIS_EVOLIC;
  if (method.compare("NN31IC") == 0)  return BASIS_NN31IC;
  if (method.compare("LUX") == 0)     return BASIS_LUX;
  if (method.compare("NN31ICQED")==0) return BASIS_NN31ICQED;
  if (method.compare("NSR")==0) return BASIS_NSR;

  cerr << "getFitBasisType Error: Invalid parametrization type: "<<method;
  exit(-1);

  return BASIS_UNDEF;
}

// ******************** GSL integration functions ****************************
/**
 * @brief GSL error handler
 * @param reason
 * @param file
 * @param line
 * @param gsl_errno
 */
void nnpdf_GSLhandler (const char * reason,
                       const char * file,
                       int line,
                       int gsl_errno)
{
  cerr << "GSL Error: "<<reason<<endl;
  return;
}

/**
 * \param filename the configuration file
 */
NNPDFSettings::NNPDFSettings(const string &folder):
  fFileName(""),
  fPDFName(""),
  fResultsDir(""),
  fTheoryDir(""),
  fGSLWork(NULL)
{
  // Read current PDF grid name from file.
  Splash();

  // Understand if filename string is file or directory
  struct stat s;
  if(stat(folder.c_str(), &s) == 0)
    {
      if( s.st_mode & S_IFDIR)
        {
          // if folder contains a trailing / remove it
          fFileName = (folder.back() == '/') ? folder.substr(0, folder.length()-1) : folder;
          const int firstindex = (int) fFileName.find_last_of("/") + 1;
          fPDFName = fFileName.substr(firstindex, fFileName.length()-firstindex);
          fResultsDir = fFileName;
          fFileName += "/filter.yml";
        }
      else if (s.st_mode & S_IFREG)
        throw NNPDF::FileError("NNPDFSettings::NNPDFSettings",
                               "This program takes a configuration folder instead of a file!");
      else
        throw NNPDF::FileError("NNPDFSettings::NNPDFSettings",
                               "Configuration folder not recognized: " + folder);
    }
  else
    throw NNPDF::FileError("NNPDFSettings::NNPDFSettings",
                           "Configuration folder not found: " + folder);

  // Load yaml file
  try {
    fConfig = YAML::LoadFile(fFileName);
  } catch(YAML::BadFile &e) {
    throw FileError("NNPDFSettings::NNPDFSettings", "runcard not found: " + fFileName);
  }

  // Check for theory ID
  const int theoryID = Get("theory","theoryid").as<int>();
  if ( theoryID < 0) throw RangeError("NNPDFSettings::NNPDFSettings", "Invalid Theory ID");

  stringstream td;
  td << "theory_" << theoryID;
  fTheoryDir = td.str();

  // load theory map
  IndexDB db(get_data_path() + "/theory.db", "theoryIndex");
  db.ExtractMap(theoryID, APFEL::kValues, fTheory);

  cout << "==== Theory summary" << endl;
  for (size_t i = 0; i < APFEL::kValues.size(); i++)
    cout << "- " << APFEL::kValues[i] << " : " << fTheory.at(APFEL::kValues[i]) << endl;

  // check basis in yaml file
  CheckBasis();

  // Load GA parameters
  LoadGA();

  // Init Random Number generator
  RandomGenerator::InitRNG(Get("fitting","rngalgo").as<int>(),
                           Get("fitting","seed").as<unsigned long int>());

  // Allocate Integrator workspace
  // Error handling
  if (Get("debug").as<bool>())
    fGSL_old_handler=gsl_set_error_handler (&nnpdf_GSLhandler);
  else
    fGSL_old_handler=gsl_set_error_handler_off();
  fGSLWork = gsl_integration_workspace_alloc (10000);

  // Load experiments
  LoadExperiments();

  // Load positivity sets
  LoadPositivities();
}

/**
 * Destroys the variables
 */
NNPDFSettings::~NNPDFSettings()
{
  gsl_integration_workspace_free (fGSLWork);
  gsl_set_error_handler(fGSL_old_handler);
}

YAML::Node NNPDFSettings::Get(const string& item) const
{
  if (!fConfig[item])
    {
      cerr << Colour::FG_RED << "\nNNPDFSettings::Get error: item not available " << item << Colour::FG_DEFAULT << endl;
      exit(-1);
    }
  return fConfig[item];
}

YAML::Node NNPDFSettings::Get(const string& node, const string& item) const
{
  if (!fConfig[node][item])
    {
      cerr << Colour::FG_RED << "\nNNPDFSettings::Get error: item not available " << node << " " << item << Colour::FG_DEFAULT << endl;
      exit(-1);
    }
  return fConfig[node][item];
}

bool NNPDFSettings::Exists(const string &item) const
{
  return !fConfig[item] ? false : true;
}

bool NNPDFSettings::Exists(const string &node, const string &item) const
{
  return !fConfig[node][item] ? false : true;
}

YAML::Node NNPDFSettings::GetPlotting(const string& item) const
{
  if (!fPlotting[item])
    {
      cerr << Colour::FG_RED << "\nNNPDFSettings::Get error: item not available " << item << Colour::FG_DEFAULT << endl;
      exit(-1);
    }
  return fPlotting[item];
}

/**
 * @brief NNPDFSettings::GetdataInfo
 * @param setname
 * @return
 */
vector<string> NNPDFSettings::GetDataInfo(const string &setname, filterType useFiltered) const
{
  // Target data directory, if useFiltered is false, read from global data
  string targetpath;
  if (!useFiltered)
    targetpath = get_data_path()  + "/commondata" ;
  else
    targetpath = fResultsDir + "/filter/" + setname ;

  vector<string> basepath(2);
  basepath[0] = targetpath + "/DATA_"+ setname + ".dat";
  basepath[1] = targetpath + "/systypes/SYSTYPE_"+ setname + "_" +
                GetSetInfo(setname).tSysOpt + ".dat";

  return basepath;
}

/**
 * @brief NNPDFSettings::GetdataInfo
 * @param setname
 * @return
 */
vector<int> NNPDFSettings::GetDataMask(const string &setname, filterType useFiltered) const
{
  // Target data directory, if useFiltered is false, read from global data
  vector<int> mask;
  if (useFiltered)
    {
      string targetpath = fResultsDir + "/filter/" + setname + "/FKMASK_" + setname + ".dat";
      ifstream f(targetpath, ios::in);

      if (f.good())
        {
          cout << Colour::FG_YELLOW << "NNPDFSettings::GetDataMask: reading mask for " << setname << Colour::FG_DEFAULT << endl;
          int v;
          while(f >> v) mask.push_back(v);
          f.close();
        }
      else
        cout << Colour::FG_YELLOW << "NNPDFSettings::GetDataMask warning: no filtered points for " << setname << Colour::FG_DEFAULT << endl;
    }
  return mask;
}

/**
 * \param setname the name of the dataset under search
 */
DataSetInfo const& NNPDFSettings::GetSetInfo(string const& setname) const
{
  std::hash<std::string> str_hash;
  const size_t hashval = str_hash(setname);

  map<int,DataSetInfo>::const_iterator iMap = fDataSetInfo.find(hashval);
  if (iMap != fDataSetInfo.end())
    return (*iMap).second;
  else
  {
    cerr << Colour::FG_RED << "NNPDFSettings::GetSetInfo error: Cannot find Set info under: " << setname << Colour::FG_DEFAULT << endl;
    exit(-1);
  }

  return (*iMap).second;
}

/**
 * \param posname the name of the positivity set under search
 */
PosSetInfo const& NNPDFSettings::GetPosInfo(string const& posname) const
{
  std::hash<std::string> str_hash;
  const size_t hashval = str_hash(posname);

  map<int,PosSetInfo>::const_iterator iMap = fPosSetInfo.find(hashval);
  if (iMap != fPosSetInfo.end())
    return (*iMap).second;
  else
  {
    cerr << Colour::FG_RED << "NNPDFSettings::GetPosInfo error: Cannot find PosSet info under: " << posname << Colour::FG_DEFAULT << endl;
    exit(-1);
  }

  return (*iMap).second;
}

// Verify configuration file is unchanged w.r.t filter.log
void NNPDFSettings::VerifyConfiguration() const
{
  cout <<endl;
  cout << Colour::FG_YELLOW << " ----------------- Veriying Configuration against Filter ----------------- "<<endl << Colour::FG_DEFAULT <<endl;;
  string target = fResultsDir + "/filter.yml";
  string filter = fResultsDir + "/md5";

  ifstream targetConfig;
  targetConfig.open(target.c_str());

  ifstream filterConfig;
  filterConfig.open(filter.c_str());

  if (!targetConfig.good())
  {
    cerr << Colour::FG_RED << "NNPDFSettings::VerifyConfiguration Error - Cannot find current config file log." << Colour::FG_DEFAULT << endl;
    cerr << "Search path: "<<target<<endl;
    exit(-1);
  }

  if (!filterConfig.good())
  {
    cerr << Colour::FG_RED << "NNPDFSettings::VerifyConfiguration Error - Cannot find filter config file log." << Colour::FG_DEFAULT << endl;
    cerr << "Search path: "<<filter<<endl;
    exit(-1);
  }

  string md5;
  filterConfig >> md5;

  MD5 targetHash;
  targetHash.update(targetConfig);
  targetHash.finalize();

  cout << "  Current Log MD5: "<<targetHash.hexdigest()<<endl;
  cout << "  Filter  Log MD5: "<<md5<<endl;

  if (targetHash.hexdigest() == md5)
  {
    cout << endl<< Colour::FG_GREEN << " ----------------- Configuration Log Verification: PASSED -----------------"<< Colour::FG_DEFAULT <<endl<<endl;
  }
  else
  {
    cerr << Colour::FG_RED << endl << " ----------------- Configuration Log Verification: FAILED -----------------" << Colour::FG_DEFAULT << endl << endl;
    cerr << "-- Configuration has been modified since last filter run"<<endl;
    cerr << "-- Please rerun the filter"<<endl<<endl;
    exit(-1);
  }
}


/**
 * @brief NNPDFSettings::CheckParam
 * @param table
 */

bool NNPDFSettings::CheckParam(std::string const& param, double const& p1, double const& p2) const
{
    if ( fabs( p1 - p2 ) > 1e-8  )
    {
      cerr << Colour::FG_RED << "NNPDFSettings::VerifyFK Error: FastKernel Table "
           <<" does not satisfy global " <<param <<": " << p1<< " vs  "<< p2 << Colour::FG_DEFAULT << endl;
      return false;
    }

    return true;
}

bool NNPDFSettings::CheckParam(std::string const& param, std::string const& p1, std::string const& p2) const
{
    // Exempt comments from parameter check
    if (param == "Comments")
      return true;

    if ( p1 != p2 )
    {
      cerr << Colour::FG_RED << "NNPDFSettings::VerifyFK Error: FastKernel Table "
           <<" does not satisfy global " <<param <<": " << p1<< " vs  "<< p2 << Colour::FG_DEFAULT << endl;
      return false;
    }

    return true;
}

/**
 * @brief NNPDFSettings::VerifyFK
 * @param table
 */
void NNPDFSettings::VerifyFK(FKTable * const &table) const
{
  const NNPDF::FKHeader::section TI = NNPDF::FKHeader::THEORYINFO;

  bool pV = true;
  for (size_t i = 0; i < APFEL::kValues.size(); i++)
    pV = !( !pV || !CheckParam(APFEL::kValues[i], GetTheory(APFEL::kValues[i]), table->GetTag(TI, APFEL::kValues[i]) )  );

  if (!pV) throw RuntimeException("NNPDFSettings::VerifyFK","mismatch between db and FKTable");
}

/**
 * @brief NNPDFSettings::PrintConfiguration print the yaml file to disk
 * @param filename
 */
void NNPDFSettings::PrintConfiguration(const string& filename) const
{
  fstream i( fFileName.c_str(), ios::in | ios::binary);
  fstream f( (fResultsDir + "/" + filename).c_str(), ios::out | ios::binary);
  if (i.fail() || f.fail()) { throw FileError("NNPDFSettings::PrintConfiguration","file failed."); }
  f << i.rdbuf();
  f.close();
  i.close();
}

/**
 * @brief NNPDFSettings::PrintConfiguration
 * @param filename
 */
void NNPDFSettings::PrintTheory(const string& filename) const
{
  fstream f( (fResultsDir + "/" + filename).c_str(), ios::out | ios::binary);
  if (f.fail()) { throw FileError("NNPDFSettings::PrintConfiguration","file failed."); }
  for (size_t i = 0; i < APFEL::kValues.size(); i++)
    f << APFEL::kValues[i] << "\t: " << fTheory.at(APFEL::kValues[i]) << endl;
  f.close();
}

/**
 * @brief NNPDFSettings::Splash print the NNPDF splash
 */
void NNPDFSettings::Splash() const
{
  cout << Colour::FG_BLUE << endl;
  cout << "  ███╗   ██╗███╗   ██╗██████╗ ██████╗ ███████╗ " << endl;
  cout << "  ████╗  ██║████╗  ██║██╔══██╗██╔══██╗██╔════╝ " << endl;
  cout << "  ██╔██╗ ██║██╔██╗ ██║██████╔╝██║  ██║█████╗   " << endl;
  cout << "  ██║╚██╗██║██║╚██╗██║██╔═══╝ ██║  ██║██╔══╝ " << endl;
  cout << "  ██║ ╚████║██║ ╚████║██║     ██████╔╝██║ " << endl;
  cout << "  ╚═╝  ╚═══╝╚═╝  ╚═══╝╚═╝     ╚═════╝ ╚═╝ 2012-2015" << Colour::FG_DEFAULT <<endl;
  cout << "  ____git____ : " << SVN_REV << ", __coredevs__ : N.H., S.C.\n" << endl;
  cout << "  Convolution Alignment Target: "<< convoluteAlign << endl;
  cout << endl;
}

/**
 * @brief NNPDFSettings::LoadExperiments parser for the experiment
 * This is a possible way to preceed which simplifies the data structure
 */
void NNPDFSettings::LoadExperiments()
{
  YAML::Node exps = fConfig["experiments"];

  // loop over experiments
  for (int i = 0; i < (int) exps.size(); i++)
    {
      fExpName.push_back(exps[i]["experiment"].as<string>());
      vector<string> nsetname;
      if (exps[i]["datasets"].size() == 0) { cerr << Colour::FG_RED << "NNPDFSettings::LoadExperiments error: experiment " << exps[i]["experiment"] << " has no datasets!" << Colour::FG_DEFAULT  << endl; exit(-1); }

      // loop over datasets
      YAML::Node dsets = exps[i]["datasets"];
      for (int j = 0; j < (int) dsets.size(); j++)
        {
          const string setname = dsets[j]["dataset"].as<string>();
          const real setfrac = dsets[j]["frac"].as<real>();
          string setsys;
          try { setsys = dsets[j]["sys"].as<string>(); }
          catch (...) { setsys = "DEFAULT";}

          // Read C-factor sources
          std::vector<string> cfactors;
          YAML::Node cfac = dsets[j]["cfac"];
          for(size_t k=0; k<cfac.size(); k++)
          {
            std::stringstream cfs; cfs << cfac[k];
            cfactors.push_back(cfs.str());
          }

          // Generate hash of setname
          std::hash<std::string> str_hash;
          const size_t hashval = str_hash(setname);

          DataSetInfo info = {setname, setsys, setfrac, cfactors};
          map<int,DataSetInfo>::const_iterator iMap = fDataSetInfo.find(hashval);
          if (iMap != fDataSetInfo.end()) { cerr << Colour::FG_RED << "NNPDFSettings::LoadExperiments error: hash collision for set: " << setname << Colour::FG_DEFAULT << endl; exit(-1); }
          else { fDataSetInfo.insert(make_pair(hashval, info)); }
          nsetname.push_back(setname);
          fSetName.push_back(setname);
        }
      fExpSetName.push_back(nsetname);
    }
}

/**
 * @brief NNPDFSettings::LoadPositivities same as LoadExperiments but
 * for positivity observables
 */
void NNPDFSettings::LoadPositivities()
{
  YAML::Node pos = fConfig["positivity"]["posdatasets"];

  // loop over positivity obs
  for (int i = 0; i < (int) pos.size(); i++)
    {
      const string posname = pos[i]["dataset"].as<string>();
      const real poslambda = pos[i]["poslambda"].as<real>();

      fPosName.push_back(posname);
      PosSetInfo info = {posname, poslambda};

      // Generate hash of setname
      std::hash<std::string> str_hash;
      const size_t hashval = str_hash(posname);

      map<int,PosSetInfo>::const_iterator iMap = fPosSetInfo.find(hashval);
      if (iMap != fPosSetInfo.end()) { cerr << Colour::FG_RED << "NNPDFSettings::LoadPositivity error: hash collision for set: " << posname << Colour::FG_DEFAULT << endl; exit(-1); }
      else { fPosSetInfo.insert(make_pair(hashval, info)); }
    }
}

/**
 * @brief NNPDFSettings::LoadPlotFile
 * @param plotfile
 */
void NNPDFSettings::SetPlotFile(string const& plotfile)
{
  // read plotting options
  try {
    fPlotting = YAML::LoadFile(plotfile);
  } catch(YAML::BadFile &) {
    throw FileError("NNPDFSettings::SetPlotFile", "runcard not found: " + plotfile);
  }
}

bool NNPDFSettings::IsQED() const
{
  const basisType isqed = NNPDFSettings::getFitBasisType(Get("fitting","fitbasis").as<string>());
  if (isqed == BASIS_EVOLQED || isqed == BASIS_EVOLSQED ||
      isqed == BASIS_FLVRQED || isqed == BASIS_NN23QED ||
      isqed == BASIS_LUX || isqed == BASIS_NN31ICQED)
    return true;
  return false;
}

bool NNPDFSettings::IsIC() const
{
  const basisType isic = NNPDFSettings::getFitBasisType(Get("fitting","fitbasis").as<string>());
  if (isic == BASIS_EVOLIC || isic == BASIS_NN30IC || isic == BASIS_NN31IC ||
      isic == BASIS_LUX || isic == BASIS_NN31ICQED)
    return true;
  return false;
}

void NNPDFSettings::CheckBasis()
{
  vector<string> basis = basiselem[getFitBasisType(Get("fitting","fitbasis").as<string>())];

  if (basis.size() != Get("fitting","basis").size())
    {
      cerr << Colour::FG_RED << "NNPDFSettings::CheckBasis error, mismatch between fitbasis and basis size" << Colour::FG_DEFAULT << endl;
      exit(-1);
    }

  // check order and names
  for (int i = 0; i < (int) Get("fitting","basis").size(); i++)
    if (basis[i].compare(Get("fitting","basis")[i]["fl"].as<string>()) != 0)
      {
        cerr << Colour::FG_RED << "NNPDFSettings::CheckBasis error, mismatch between basis items, expected "
             << basis[i] << ", received " <<   Get("fitting","basis")[i]["fl"].as<string>() << Colour::FG_DEFAULT << endl;
        exit(-1);
      }
}

int NNPDFSettings::GetNFL() const
{
  return (int) Get("fitting","basis").size();
}

void NNPDFSettings::LoadGA()
{
  // for each flavor check mutation array size and probability array size
  for (int f = 0; f < GetNFL(); f++)
    {
      if (Get("fitting","basis")[f]["mutsize"].size() != Get("fitting","basis")[f]["mutprob"].size())
        {
          cerr << Colour::FG_RED << "NNPDFSettings::LoadGA error, mismatch between mutsize and mutprob for flavor "
               << Get("fitting","basis")[f]["fl"].as<string>() << Colour::FG_DEFAULT << endl;
          exit(-1);
        }
      else
        {
          FlMutProperty p = { Get("fitting","basis")[f]["mutsize"].as<vector<int> >(), Get("fitting","basis")[f]["mutprob"].as<vector<real> >()};
          fFlMutProperty.push_back(p);
        }
    }

  // Load architecture
  fArch = Get("fitting","nnodes").as<vector<int> >();
}
