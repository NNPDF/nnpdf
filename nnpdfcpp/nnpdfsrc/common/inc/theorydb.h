// $Id$
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#pragma once

#include <string>
using std::string;

class NNPDFSettings;
class sqlite3;

class Database
{
public:
  Database(const string& filename);
  ~Database();
  string query(const string& query) const;

private:
  sqlite3 *fdatabase;
};

// Parse a NNPDFSettings theory parameters into a theoryID
bool parseTheory(NNPDFSettings const& settings);
