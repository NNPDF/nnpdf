/*
 *  nnpdfdb.cc
 *  Database access routines for nnpdf applications
 * *  nathan.hartland@physics.ox.ac.uk 02/15
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
	  std::cerr << "IndexDB::IndexDB Error - Can't open database:" << sqlite3_errmsg(fDB) <<std::endl;
	  exit(-1);
	}

	// Determine number of enries
	sqlite3_stmt *statement;  
	const std::string query = "select * from "+fTable;
	const int retcode = sqlite3_prepare_v2(fDB, query.c_str(), -1, &statement, 0 );

	if ( retcode != SQLITE_OK ) 
	{
	  std::cerr << "IndexDB::IndexDB Error - " << sqlite3_errmsg(fDB)<<std::endl;
	  exit(-1);
	}

	while (sqlite3_step(statement) == SQLITE_ROW)
	  fNEntries++;

	}

	IndexDB::~IndexDB()
	{
	  sqlite3_close(fDB);
	}


	// Needs to be specialised to avoid splitting the string
	template<>
	std::string dbquery(IndexDB const& db, int const& id, std::string const& field)
	{
	  std::stringstream query;
	  query << "select " << field << " from "<<db.fTable<<" where id=" <<id;

	  sqlite3_stmt *statement;  
	  const int retcode = sqlite3_prepare_v2(db.fDB, query.str().c_str(), -1, &statement, 0 );
	      
	  if ( retcode != SQLITE_OK ) 
	  {
	    std::cerr << "IndexDB::dbQuery Error - " << sqlite3_errmsg(db.fDB)<<std::endl;
	    exit(-1);
	  }

	  int res = sqlite3_step(statement);
	  if ( res == SQLITE_ROW ) 
	  {
	    std::string retval = (char*)sqlite3_column_text(statement, 0);
	    return retval; 
	  } else
	  {
	      std::cerr << "IndexDB::dbQuery Error - " << sqlite3_errmsg(db.fDB)<<std::endl;
	      exit(1);
	  }    
	}

}
