/*
 *  nnpdfdb.cc
 *  Database access routines for nnpdf applications
 *  nathan.hartland@physics.ox.ac.uk 02/15
 */

#include "NNPDF/nnpdfdb.h"

namespace NNPDF
{
  
  IndexDB::IndexDB(std::string const& _path, std::string const& _table):
    fDB(0),
    fPath(_path),
    fTable(_table),
    fNEntries(0)
  {
    
    // Setup db connection
    if( sqlite3_open(fPath.c_str(), &fDB) )
      {
	std::stringstream err; err << "Can't open database: " << sqlite3_errmsg(fDB);
	throw FileError("IndexDB::IndexDB", err.str());
      }
    
    // Determine number of entries
    sqlite3_stmt *statement;  
    const std::string query = "select * from "+fTable;
    const int retcode = sqlite3_prepare_v2(fDB, query.c_str(), -1, &statement, 0 );
    if ( retcode != SQLITE_OK ) 
      {
	std::stringstream err; err << "SQLITE Error: " << sqlite3_errmsg(fDB);
	throw RuntimeException("IndexDB::IndexDB", err.str());
      }
    
    while (sqlite3_step(statement) == SQLITE_ROW)
      fNEntries++;
    
  }
  
  IndexDB::~IndexDB()
  {
    sqlite3_close(fDB);
  }
  
  void IndexDB::ExtractMap( const int& id, std::vector<std::string> const& keys, std::map<std::string, std::string>& map)
  {
    for (size_t i=0; i<keys.size(); i++)
      map.emplace(keys[i], dbquery<std::string>(*this, id, keys[i]));
  }
  
  std::vector<std::string> IndexDB::ExtractString( const int& id, std::vector<std::string> const& keys)
  {
    std::map<std::string, std::string> map;
    ExtractMap(id,keys,map);
    
    std::vector<std::string> v;
    for( std::map<std::string, std::string>::iterator it = map.begin(); it != map.end(); ++it )
      v.push_back( it->second );
    
    return v;
  }

  // Needs to be specialised to avoid splitting the string and returning NULL strings
  template<>
  std::string dbquery(IndexDB const& db, int const& id, std::string const& field)
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
    if (sqlite3_column_type(statement, 0) == SQLITE_NULL) // result is NULL
    {
        return std::string("");
    }
    else if ( res == SQLITE_ROW ) 
    {
    	std::string retval = (char*)sqlite3_column_text(statement, 0);
      retval.erase(retval.find_last_not_of(" \n\r\t")+1);
    	return retval; 
    } 
    else
    {
    	std::stringstream err; err << "SQLITE Error: " << sqlite3_errmsg(db.fDB);
    	throw RuntimeException("dbquery", err.str());
    }    
  }

  std::vector<int> dbmatch(IndexDB const& db, std::string const& field, std::string const& value)
  {
    std::stringstream query;
    query << "select id from "<<db.fTable<<" where "<<field<<"='"<<value<<"' order by id";

    sqlite3_stmt *statement;  
    const int retcode = sqlite3_prepare_v2(db.fDB, query.str().c_str(), -1, &statement, 0 );
    if ( retcode != SQLITE_OK ) 
    {
      std::stringstream err; err << "SQLITE Error: " << sqlite3_errmsg(db.fDB);
      throw RuntimeException("dbquery", err.str());
    }

    std::vector<int> indices;
    while ( sqlite3_step(statement) == SQLITE_ROW ) 
    {
        std::stringstream s;
        s <<(char*)sqlite3_column_text(statement, 0);
        
        int ID; s >> ID;
        indices.push_back(ID); 
    }
    return indices;
  }
  
}
