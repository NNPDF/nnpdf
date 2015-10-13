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
#include "theorydb.h"
#include "svn.h"

// Strings for config file output
static const string minString[6]   = {"UNDEFINED", "GA", "NGA"};
static const string stopString[6]  = {"UNDEFINED", "FIXEDLENGTH", "TR", "GRAD", "VAR", "LOOKBACK"};
static const string paramString[4] = {"UNDEFINED", "NN", "CHEBYSHEV", "QUADNN"};
static const string basisString[13]= {"UNDEFINED", "NN23", "NN23QED","EVOL", "EVOLQED","EVOLS",
                                      "EVOLSQED","NN30", "NN30QED","FLVR", "FLVRQED","NN30IC","EVOLIC"};

static const vector< vector<string> > basiselem = { {},
                                     {"sng","g","v","t3","ds","sp","sm"},
                                     {"sng","g","v","t3","ds","sp","sm","pht"},
                                     {"sng","g","v","v3","v8","t3","t8"},
                                     {"sng","g","v","v3","v8","t3","t8","pht"},
                                     {"sng","g","v","v8","t3","t8","ds"},
                                     {"sng","g","v","v8","t3","t8","ds","pht"},
                                     {"sng","g","v","v8","t3","t8","ds"},
                                     {"sng","g","v","t3","ds","sp","sm"},
                                     {"sng","g","v","t3","ds","sp","sm","pht"},
                                     {"g","u","ubar","d","dbar","s","sbar"},
                                     {"g","u","ubar","d","dbar","s","sbar","pht"},
                                     {"sng","g","v","t3","ds","sp","sm","cp","cm"},
                                     {"sng","g","v","v3","v8","v15","t3","t8","t15"}
                                     };

/* Convert string to enum */
minType NNPDFSettings::getFitMethod(string const& method)
{
  if (method.compare("GA") == 0)      return MIN_GA;
  if (method.compare("NGA") == 0)     return MIN_NGA;
  
  cerr << "getFitMethod Error: Invalid fit method: "<<method<<endl;
  cerr << "choices are: "<<endl;
  cerr <<" - GA (Basic GA)"<<endl;
  cerr <<" - NGA"<<endl;
  exit(-1);
  
  return MIN_UNDEF;
}

paramType NNPDFSettings::getParamType(string const& method)
{
  if (method.compare("NN") == 0)        return PARAM_NN;
  if (method.compare("CHEBYSHEV") == 0) return PARAM_CHEBYSHEV;
  if (method.compare("QUADNN") == 0)    return PARAM_QUADNN;
  
  cerr << "getParamType Error: Invalid parametrization type: "<<method;
  exit(-1);
  
  return PARAM_UNDEF;
}

stopType NNPDFSettings::getStopType(string const& method)
{
  if (method.compare("TRVAL") == 0)       return STOP_TR;
  if (method.compare("GRAD") == 0)        return STOP_GRAD;
  if (method.compare("VAR") == 0)         return STOP_VAR;
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

  cerr << "getFitBasisType Error: Invalid parametrization type: "<<method;
  exit(-1);

  return BASIS_UNDEF;
}

FNS NNPDFSettings::getVFNS(string const& vfns)
{
 
  if (vfns.compare("FFNS") == 0) return FFNS;
  if (vfns.compare("ZM-VFNS") == 0) return ZMVFNS;
  if (vfns.compare("FONLL-A") == 0) return FONLLA;
  if (vfns.compare("FONLL-B") == 0) return FONLLB;
  if (vfns.compare("FONLL-C") == 0) return FONLLC;

  cerr << Colour::FG_RED << "getVFNS Error: Invalid FNS " << vfns << endl;
  exit(-1);
}

MODEV NNPDFSettings::getMODEV(const string &modev)
{
  if (modev.compare("TRN") == 0) return TRN;
  if (modev.compare("EXP") == 0) return EXP;
  if (modev.compare("EXA") == 0) return EXA;

  cerr << Colour::FG_RED << "getMODEV Error: Invalid MODEV " << modev << endl;
  exit(-1);
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
NNPDFSettings::NNPDFSettings(const string &filename, const string &plotfile):
  fFileName(filename),
  fPDFName(""),
  fResultsDir(""),
  fTheoryDir(""),
  fGSLWork(NULL)
{  
  // Read current PDF grid name from file.
  Splash();
    
  // Get raw name
  const int firstindex  = (int) fFileName.find_last_of("/") + 1;
  const int lastindex   = (int) fFileName.find_last_of(".") - firstindex;
  fPDFName = fFileName.substr(firstindex, lastindex);

  // Load yaml file
  fConfig = YAML::LoadFile(fFileName);

  // Check for theory ID
  if ( Get("theory","theoryid").as<int>() < 0) { cerr << Colour::FG_RED << "\nInvalid Theory ID" << endl; exit(-1); }
  stringstream td;
  td << "theory_" << Get("theory","theoryid").as<int>();
  fTheoryDir = td.str();

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

  // Eventually read plotting options
  if (plotfile.size() > 0)
    fPlotting = YAML::LoadFile(plotfile);

  // Results directory stuff
  fResultsDir = resultsPath();
  
  struct stat st;
  if(stat(fResultsDir.c_str(),&st) != 0)
  {
    printf("Warning: RESULTSDIR folder is NOT present!\nCreating a new folder.\n\n");
    mkdir(fResultsDir.c_str(), 0777);
  }

  // Setup results directory
  fResultsDir += fPDFName;
  mkdir(fResultsDir.c_str(), 0777);

  parseTheory(*this);
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
      cerr << Colour::FG_RED << "\nNNPDFSettings::Get error: item not available " << item << endl;
      exit(-1);
    }
  return fConfig[item];
}

YAML::Node NNPDFSettings::Get(const string& node, const string& item) const
{
  if (!fConfig[node][item])
    {
      cerr << Colour::FG_RED << "\nNNPDFSettings::Get error: item not available " << node << " " << item << endl;
      exit(-1);
    }
  return fConfig[node][item];
}

YAML::Node NNPDFSettings::GetPlotting(const string& item) const
{
  if (!fPlotting[item])
    {
      cerr << Colour::FG_RED << "\nNNPDFSettings::Get error: item not available " << item << endl;
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
    targetpath = dataPath() + "commondata" ;
  else
    targetpath = fResultsDir + "/filter/" + setname ;

  vector<string> basepath(2);
  basepath[0] = targetpath + "/DATA_"+ setname + ".dat";
  basepath[1] = targetpath + "/systypes/SYSTYPE_"+ setname + "_" +
                std::to_string(GetSetInfo(setname).tSysOpt) + ".dat";

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
  map<int,DataSetInfo>::const_iterator iMap = fDataSetInfo.find(HashNNPDF::IntHash(setname.c_str()));
  if (iMap != fDataSetInfo.end())
    return (*iMap).second;
  else
  {
    cerr << Colour::FG_RED << "NNPDFSettings::GetSetInfo error: Cannot find Set info under: "<<setname<<endl;
    exit(-1);
  }
  
  return (*iMap).second;
}

/**
 * \param posname the name of the positivity set under search
 */
PosSetInfo const& NNPDFSettings::GetPosInfo(string const& posname) const
{
  map<int,PosSetInfo>::const_iterator iMap = fPosSetInfo.find(HashNNPDF::IntHash(posname.c_str()));
  if (iMap != fPosSetInfo.end())
    return (*iMap).second;
  else
  {
    cerr << Colour::FG_RED << "NNPDFSettings::GetPosInfo error: Cannot find PosSet info under: "<<posname<<endl;
    exit(-1);
  }
  
  return (*iMap).second;
}

// Verify configuration file is unchanged w.r.t filter.log
void NNPDFSettings::VerifyConfiguration(const string &filename)
{  
  cout <<endl;
  cout << Colour::FG_YELLOW << " ----------------- Veriying Configuration against Filter ----------------- "<<endl << Colour::FG_DEFAULT <<endl;;
  string target = fResultsDir + "/"+filename;
  string filter = fResultsDir + "/filter.yml";
  
  ifstream targetConfig;
  targetConfig.open(target.c_str());
  
  ifstream filterConfig;
  filterConfig.open(filter.c_str());
  
  if (!targetConfig.good())
  {
    cerr << Colour::FG_RED << "NNPDFSettings::VerifyConfiguration Error - Cannot find current config file log."<<endl;
    cerr << "Search path: "<<target<<endl;
    exit(-1);
  }
  
  if (!filterConfig.good())
  {
    cerr << Colour::FG_RED << "NNPDFSettings::VerifyConfiguration Error - Cannot find filter config file log."<<endl;
    cerr << "Search path: "<<filter<<endl;
    exit(-1);
  }
  
  HashNNPDF::MD5 targetHash, filterHash;
  
  targetHash.update(targetConfig);
  filterHash.update(filterConfig);
  
  targetHash.finalize();
  filterHash.finalize();
  
  cout << "  Current Log MD5: "<<targetHash.hexdigest()<<endl;
  cout << "  Filter  Log MD5: "<<filterHash.hexdigest()<<endl;

  if (filterHash.hexdigest().compare(targetHash.hexdigest()) == 0)
  {
    cout << endl<< Colour::FG_GREEN << " ----------------- Configuration Log Verification: PASSED -----------------"<< Colour::FG_DEFAULT <<endl<<endl;
  }
  else
  {
    cerr << Colour::FG_RED << endl << " ----------------- Configuration Log Verification: FAILED -----------------"<<endl<<endl;
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
           <<" does not satisfy global " <<param <<": " << p1<< " vs  "<< p2 <<endl;
      return false;
    }

    return true;
}

bool NNPDFSettings::CheckParam(std::string const& param, std::string const& p1, std::string const& p2) const
{
    if ( p1 != p2 )
    {
      cerr << Colour::FG_RED << "NNPDFSettings::VerifyFK Error: FastKernel Table "
           <<" does not satisfy global " <<param <<": " << p1<< " vs  "<< p2 <<endl;
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
  pV = !( !pV || !CheckParam("ptord", Get("theory","ptord").as<int>(), table->GetTag<int>(TI, "PTO") )  );

  pV = !( !pV || !CheckParam("ptord", Get("theory","alphas").as<double>(), table->GetTag<double>(TI, "alphas") )  );
  pV = !( !pV || !CheckParam("ptord", Get("theory","q20").as<double>(), pow(table->GetTag<double>(TI, "Q0"),2.0) )  );
  pV = !( !pV || !CheckParam("ptord", Get("theory","qref").as<double>(), table->GetTag<double>(TI, "Qref") )  );

  // String->Bool
  //pV = !( !pV || !CheckParam("ptord", Get("theory","msbar").as<bool>(), table->GetTag<bool>(TI, "HQ") )  );
  pV = !( !pV || !CheckParam("ptord", Get("theory","mc").as<double>(), table->GetTag<double>(TI, "mc") )  );
  pV = !( !pV || !CheckParam("ptord", Get("theory","mb").as<double>(), table->GetTag<double>(TI, "mb") )  );
  pV = !( !pV || !CheckParam("ptord", Get("theory","mt").as<double>(), table->GetTag<double>(TI, "mt") )  );

  if (!pV)
    exit(-1);

}

/**
 * @brief NNPDFSettings::PrintConfiguration print the yaml file to disk
 * @param filename
 */
void NNPDFSettings::PrintConfiguration(const string& filename) const
{
  fstream i( fFileName.c_str(), ios::in | ios::binary);
  fstream f( (fResultsDir + "/" + filename).c_str(), ios::out | ios::binary);
  if (i.fail() || f.fail()) { cerr << "fail" << endl; }
  f << i.rdbuf();
  f.close();
  i.close();
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
  cout << "  ____svn____ : " << SVN_REV << ", __coredevs__ : N.H., S.C.\n" << endl;
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
      if (exps[i]["datasets"].size() == 0) { cerr << Colour::FG_RED << "NNPDFSettings::LoadExperiments error: experiment " << exps[i]["experiment"] << " has no datasets!" << endl; exit(-1); }

      // loop over datasets
      YAML::Node dsets = exps[i]["datasets"];
      for (int j = 0; j < (int) dsets.size(); j++)
        {
          const string setname = dsets[j]["dataset"].as<string>();
          const real setfrac = dsets[j]["frac"].as<real>();
          const int setsys   = dsets[j]["sys"].as<int>();

          // Read C-factor sources
          std::vector<string> cfactors;
          YAML::Node cfac = dsets[j]["cfac"];
          for(size_t k=0; k<cfac.size(); k++)
          {
            std::stringstream cfs; cfs << cfac[k];
            cfactors.push_back(cfs.str());
          }

          DataSetInfo info = {setname, setsys, setfrac, cfactors};

          map<int,DataSetInfo>::const_iterator iMap = fDataSetInfo.find(HashNNPDF::IntHash(setname.c_str()));
          if (iMap != fDataSetInfo.end()) { cerr << Colour::FG_RED << "NNPDFSettings::LoadExperiments error: hash collision for set: " << setname << endl; exit(-1); }
          else { fDataSetInfo.insert(make_pair(HashNNPDF::IntHash(setname.c_str()), info)); }
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

      map<int,PosSetInfo>::const_iterator iMap = fPosSetInfo.find(HashNNPDF::IntHash(posname.c_str()));
      if (iMap != fPosSetInfo.end()) { cerr << Colour::FG_RED << "NNPDFSettings::LoadPositivity error: hash collision for set: " << posname << endl; exit(-1); }
      else { fPosSetInfo.insert(make_pair(HashNNPDF::IntHash(posname.c_str()), info)); }
    }
}

bool NNPDFSettings::IsQED() const
{
  const basisType isqed = NNPDFSettings::getFitBasisType(Get("fitting","fitbasis").as<string>());
  if (isqed == BASIS_EVOLQED || isqed == BASIS_EVOLSQED || isqed == BASIS_FLVRQED || isqed == BASIS_NN23QED)
    return true;
  return false;
}

bool NNPDFSettings::IsIC() const
{
  const basisType isic = NNPDFSettings::getFitBasisType(Get("fitting","fitbasis").as<string>());
  if (isic == BASIS_EVOLIC || isic == BASIS_NN30IC)
    return true;
  return false;
}

void NNPDFSettings::CheckBasis()
{
  vector<string> basis = basiselem[getFitBasisType(Get("fitting","fitbasis").as<string>())];

  if (basis.size() != Get("fitting","basis").size())
    {
      cerr << Colour::FG_RED << "NNPDFSettings::CheckBasis error, mismatch between fitbasis and basis size" << endl;
      exit(-1);
    }

  // check order and names
  for (int i = 0; i < (int) Get("fitting","basis").size(); i++)
    if (basis[i].compare(Get("fitting","basis")[i]["fl"].as<string>()) != 0)
      {
        cerr << Colour::FG_RED << "NNPDFSettings::CheckBasis error, mismatch between basis items, expected "
             << basis[i] << ", received " <<   Get("fitting","basis")[i]["fl"].as<string>() << endl;
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
