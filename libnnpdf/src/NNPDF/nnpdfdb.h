#pragma once
/*
 *  nnpdfdb.h
 *  Database access routines for nnpdf family applications
 * *  nathan.hartland@physics.ox.ac.uk 02/15
 */

#include <sqlite3.h> 
#include <iostream>
#include <sstream>
#include <cstdlib>

#include "exceptions.h"

namespace NNPDF
{

  class IndexDB;
  template<class T>
  T dbquery(IndexDB const& db, int const& id, std::string const& field);

  // Index-based databases
  class IndexDB 
  {
  public:
    IndexDB(std::string const& _path, std::string const& _table);
    ~IndexDB();

    const int& GetNEntries() const {return fNEntries;};

  private:
    sqlite3* fDB;

    const std::string fPath;
    const std::string fTable;

    int fNEntries;

    template<class T>
    friend T dbquery(IndexDB const& db, int const& id, std::string const& field);
    static void HandleRetCode(const int retcode);

  };

    // Return a queried value
    template<class T>
    T dbquery(IndexDB const& db, int const& id, std::string const& field)
    {
      std::stringstream query;
      query << "select " << field << " from "<<db.fTable<<" where id=" <<id;

      sqlite3_stmt *statement;  
      const int retcode = sqlite3_prepare_v2(db.fDB, query.str().c_str(), -1, &statement, 0 );
          
      if ( retcode != SQLITE_OK ) 
      {
        std::stringstream err; err << "SQLITE Error: " << sqlite3_errmsg(db.fDB);
        throw RuntimeException("dbquery", err.str());
      }

      int res = sqlite3_step(statement);
      if ( res == SQLITE_ROW ) 
      {
        std::stringstream s;
        s <<(char*)sqlite3_column_text(statement, 0);
        
        // return value 
        T retval;
        s >> retval;

        return retval; 
      } else
      {
        std::stringstream err; err << "SQLITE Error: " << sqlite3_errmsg(db.fDB);
        throw RuntimeException("dbquery", err.str());
      }    
    }

    // Needs to be specialised to avoid splitting the string
    template<>
    std::string dbquery(IndexDB const& db, int const& id, std::string const& field);

}
