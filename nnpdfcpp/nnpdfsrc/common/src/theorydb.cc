// $Id$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include "theorydb.h"
#include <iostream>
#include <cmath>
#include <sqlite3.h>

#include "common.h"
#include "nnpdfsettings.h"
#include <NNPDF/fastkernel.h>
#include <NNPDF/exceptions.h>
using std::exit;

Database::Database(string const& filename)
{
  fdatabase = NULL;
  if(sqlite3_open(filename.c_str(), &fdatabase) != SQLITE_OK)
    {
      cout << "Error: cannot read database" << endl;
    }
}

Database::~Database()
{
  sqlite3_close(fdatabase);
}

string Database::query(string const& query) const
{
  sqlite3_stmt *statement;
  string results;

  if(sqlite3_prepare_v2(fdatabase, query.c_str(), -1, &statement, 0) == SQLITE_OK)
  {
    int result = 0;
    result = sqlite3_step(statement);
    results = (char*)sqlite3_column_text(statement, 0);
    sqlite3_finalize(statement);
  }

  string error = sqlite3_errmsg(fdatabase);
  if(error != "not an error") cout << query << " " << error << endl;

  return results;
}


bool parseTheory(const NNPDFSettings &settings)
{
  // opening the theory db
  Database *db = new Database(dataPath()+"theory.db");

  // allocating the theoryid
  const string thid = settings.Get("theory","theoryid").as<string>();

  string result = db->query("SELECT PTO FROM TheoryIndex WHERE ID=" + thid);
  if (settings.Get("theory","ptord").as<int>() != atoi(result.c_str()))    
    throw NNPDF::EvaluationError("parseTheory","mismatch theory.db and config file for PTORD");

  result = db->query("SELECT alphas FROM TheoryIndex WHERE ID=" + thid);
  if ( fabs(settings.Get("theory","alphas").as<double>() - atof(result.c_str()) ) > 1E-5 )
    throw NNPDF::EvaluationError("parseThery","mismatch theory.db and config file for alphas");

  result = db->query("SELECT Qref FROM TheoryIndex WHERE ID=" + thid);
  if (fabs(settings.Get("theory","qref").as<double>() - atof(result.c_str())) > 1e-5 )
    throw NNPDF::EvaluationError("parseThery","mismatch theory.db and config file for Qref");

  result = db->query("SELECT Q0 FROM TheoryIndex WHERE ID=" + thid);
  if (fabs(sqrt(settings.Get("theory","q20").as<double>()) - atof(result.c_str())) > 1e-3)
    throw NNPDF::EvaluationError("parseThery","mismatch theory.db and config file for Q0");

  result = db->query("SELECT FNS FROM TheoryIndex WHERE ID=" + thid);
  if (settings.Get("theory","vfns").as<string>().compare(result) != 0)
    throw NNPDF::EvaluationError("parseThery","mismatch theory.db and config file for VFNS");

  result = db->query("SELECT mc FROM TheoryIndex WHERE ID=" + thid);
  if (fabs(settings.Get("theory","mc").as<double>() - atof(result.c_str())) > 1e-3)
    throw NNPDF::EvaluationError("parseThery","mismatch theory.db and config file for Mc");

  result = db->query("SELECT mb FROM TheoryIndex WHERE ID=" + thid);
  if (fabs(settings.Get("theory","mb").as<double>() - atof(result.c_str())) > 1e-3)
    throw NNPDF::EvaluationError("parseThery","mismatch theory.db and config file for Mb");

  result = db->query("SELECT mt FROM TheoryIndex WHERE ID=" + thid);
  if (fabs(settings.Get("theory","mt").as<double>() - atof(result.c_str())) > 1e-3)
    throw NNPDF::EvaluationError("parseThery","mismatch theory.db and config file for Mt");

    result = db->query("SELECT Qmc FROM TheoryIndex WHERE ID=" + thid);
  if (fabs(settings.Get("theory","Qc").as<double>() - atof(result.c_str())) > 1e-3)
    throw NNPDF::EvaluationError("parseThery","mismatch theory.db and config file for Qc");

  result = db->query("SELECT Qmb FROM TheoryIndex WHERE ID=" + thid);
  if (fabs(settings.Get("theory","Qb").as<double>() - atof(result.c_str())) > 1e-3)
    throw NNPDF::EvaluationError("parseThery","mismatch theory.db and config file for Qb");

  result = db->query("SELECT Qmt FROM TheoryIndex WHERE ID=" + thid);
  if (fabs(settings.Get("theory","Qt").as<double>() - atof(result.c_str())) > 1e-3)
    throw NNPDF::EvaluationError("parseThery","mismatch theory.db and config file for Qt");

  result = db->query("SELECT MaxNfPDF FROM TheoryIndex WHERE ID=" + thid);
  if ( settings.Get("theory","nf_pdf").as<int>() != atoi(result.c_str()))
    throw NNPDF::EvaluationError("parseThery","mismatch theory.db and config file for NF_PDF");

  result = db->query("SELECT MaxNfAs FROM TheoryIndex WHERE ID=" + thid);
  if ( settings.Get("theory","nf_as").as<int>() != atoi(result.c_str()))
    throw NNPDF::EvaluationError("parseThery","mismatch theory.db and config file for NF_AS");

  result = db->query("SELECT HQ FROM TheoryIndex WHERE ID=" + thid);
  if ( settings.Get("theory","msbar").as<bool>() == true)
    if (result.compare("MSBAR") != 0)
      throw NNPDF::EvaluationError("parseThery","mismatch theory.db and config file for MSBAR");

  result = db->query("SELECT ModEv FROM TheoryIndex WHERE ID=" + thid);
  if (settings.Get("theory","modev").as<string>().compare(result) != 0)
    throw NNPDF::EvaluationError("parseThery","mismatch theory.db and config file for ModEv");

  result = db->query("SELECT SxRes FROM TheoryIndex WHERE ID=" + thid);
  if ( settings.Get("theory","sxres").as<bool>() != atoi(result.c_str()))
    throw NNPDF::EvaluationError("parseThery","mismatch theory.db and config file for SxRes");

  result = db->query("SELECT SxOrd FROM TheoryIndex WHERE ID=" + thid);
  if (settings.Get("theory","sxord").as<string>().compare(result) != 0)
    throw NNPDF::EvaluationError("parseThery","mismatch theory.db and config file for SxOrd");

  delete db;
  return true;
}
