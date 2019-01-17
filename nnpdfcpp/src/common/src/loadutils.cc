// $Id$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include "loadutils.h"
#include "nnpdfsettings.h"

/**
 * @brief LoadDataSet
 * @param settings
 * @param setname
 * @param useFilter
 * @return
 */
DataSet LoadDataSet(NNPDFSettings const &settings, std::string const &setname,
                    filterType useFilter)
{
  // get path for setname
  vector<string> fs = settings.GetDataInfo(setname, useFilter);
  vector<int> mask = settings.GetDataMask(setname, useFilter);
  auto weight = settings.GetSetInfo(setname).weight;

  // allocate commondata
  CommonData cd = CommonData::ReadFile(fs[0], fs[1]);

  // Load fkset of tables
  FKSet fk = LoadFK(settings, setname);

  // return dataset
  if (mask.size() > 0) {
    fk = FKSet(fk, mask);
  }
  return DataSet(cd, fk, weight);
}

/**
 * @brief LoadPositivitySet
 * @param settings
 * @param posname
 * @param lambda
 * @return
 */
PositivitySet LoadPositivitySet(NNPDFSettings const& settings, std::string const& posname, real const& lambda)
{
  // allocate commondata
  CommonData cd = CommonData::ReadFile(get_data_path() + "/commondata/DATA_" + posname + ".dat",
                                       get_data_path() + "/commondata/systypes/SYSTYPE_" + posname + "_DEFAULT.dat");
  // Load fkset of tables
  FKTable fk(get_data_path() + "/" + settings.GetTheoryDirectory() + "/fastkernel/FK_" + posname + ".dat");

  // return positivity set
  return PositivitySet(cd,fk,lambda);
}

/**
 * @brief LoadFK
 * @param settings
 * @param setname
 * @return
 */
FKSet LoadFK(NNPDFSettings const& settings,
             std::string const& setname)
{
  const std::string theoryDir = settings.GetTheoryDirectory();
  const std::string theoryPath = get_data_path() + "/" + theoryDir + "/";

  stringstream cfilename("");
  cfilename << theoryPath << "compound/"
  << "FK_" << setname << "-COMPOUND.dat";

  // allocating FKtables
  vector<FKTable*> nFK;

  int NSigma = 0;
  SigmaOp op = FKSet::parseOperator("NULL");

  ifstream  compound(cfilename.str().c_str());
  if (compound)
  {
    string line;
    while (getline(compound, line))
    {

      vector<string> sline = split(line);
      if (sline.size()==0)
        continue;

      if (sline[0]=="OP:")
      {
        if (sline.size() > 1)
          op = FKSet::parseOperator(sline[1]);
      }

      if (sline[0]=="FK:")
        if (sline.size()>1)
        {
          stringstream sigfilename("");
          sigfilename << theoryPath
          << "fastkernel/"
          << sline[1];

          // load cfactors
          vector<string> cfactors;
          for (int i = 0; i < (int) settings.GetSetInfo(setname).tCFactors.size(); i++)
          {
            const string cname = settings.GetSetInfo(setname).tCFactors[i];
            const string fname = sline[1].substr(3,sline[1].length());
            cfactors.push_back(get_data_path()+ "/" + theoryDir + "/cfactor/CF_"+cname+"_" + fname);
            cout << Colour::FG_BLUE << "-- Reading "+cname+" C-factors from: " << cfactors[i] << Colour::FG_DEFAULT << endl;
          }

          // Read FK Table
          FKTable *newTable = new FKTable(sigfilename.str(), cfactors);
          nFK.push_back(newTable);
          NSigma++;
      }
    }

    if (NSigma==0)
    {
      cerr << "DataSet::ReadFK Error: No FastKernel Grids loaded"<<endl;
      exit(-1);
    }

  } else {

    stringstream sigfilename("");

    sigfilename << theoryPath
    << "fastkernel/"
    << "FK_" << setname << ".dat";

    // load cfactors
    vector<string> cfactors;
    for (int i = 0; i < (int) settings.GetSetInfo(setname).tCFactors.size(); i++)
    {
      const string cname = settings.GetSetInfo(setname).tCFactors[i];
      cfactors.push_back(get_data_path() + "/" + theoryDir + "/cfactor/CF_"+cname+"_" + setname + ".dat");
      cout << Colour::FG_BLUE << "-- Reading "+cname+" C-factors from: " << cfactors[i] << Colour::FG_DEFAULT << endl;
    }

    // Read FK Table
    FKTable *newTable = new FKTable(sigfilename.str(), cfactors);

    nFK.push_back(newTable);
    NSigma++;
  }

  // Verify FastKernel Tables
  for (int i=0; i< NSigma; i++) settings.VerifyFK(nFK[i]);

  return FKSet(op, nFK);
}
