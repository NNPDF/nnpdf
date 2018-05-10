// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#include "NNPDF/common.h"
#include "NNPDF/logger.h"
#include "NNPDF/exceptions.h"
#include "NNPDF/utils.h"

#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <cstdlib>
#include <functional>

using std::cout;
using std::cerr;
using std::endl;
using std::string;

namespace NNPDF
{

  // Symbol for logmanager singleton
  LogManager *LogManager::logInstance = 0;

  /**
   * @brief LogManager::LogManager()
   * Private constructor for the Log file manager class
   */
  LogManager::LogManager():
  fMap(),
  fBasePath("")
  {

  }

  /**
   * @brief LogManager::~LogManager()
   * LogManager singleton destructor
   */
  LogManager::~LogManager()
  {
    fMap.clear();
  }

  /**
   * @brief LogManager::AddLogger
   * Adds a new logger class to the manager.
   * Logger classes are kept in a hash map, indexed by the hash of logname
   * @param logname The identifier for the log, e.g NNFIT or GAMinimizer
   * @param filename The target filename to which the log is written out. Relative to base path in LogManager
   */
  void LogManager::AddLogger(string const& logname, string const& filename)
  {
    LogManager* LM = GetLM();

    std::hash<std::string> str_hash;
    const size_t hashval = str_hash(logname);

    // Collision
    if (LM->fMap.find(hashval)!=LM->fMap.end())
    {
      cerr << "LogManager::AddLogger Error - hash collision for new log: \'" <<logname<<"\' with hash: "<< hashval<<". ";
      throw LogError("LogManager::AddLogger","Log already exists!");
    }

    // Insert new logger
    string targetfile = LM->fBasePath + "/" + filename;
    LM->fMap.insert(std::pair<size_t,Logger>(hashval,Logger(logname, targetfile )));
  }

  /**
   * @brief LogManager::GetLogger
   * Returns a reference to specific logger, identified by name
   * @param logname The identifier for the log, e.g NNFIT or GAMinimizer
   */
  Logger& LogManager::GetLogger(string const& logname)
  {
    LogManager* LM = GetLM();
    std::hash<std::string> str_hash;
    const size_t hashval = str_hash(logname);

    // Collision
    LogMap::iterator iLog = LM->fMap.find(hashval);
    if ( iLog == LM->fMap.end())
    {
      cerr << "LogManager::GetLogger Error - log: \'" <<logname<<"\' does not exist!"<<endl;
      throw LogError("LogManager::AddLogger","Log does not exists!");
    }

    return (*iLog).second;
  }

  /**
   * @brief LogManager::AddLogEntry
   * Adds a new entry into the target Logger
   * @param logname The identifier for the target log, e.g NNFIT or GAMinimizer
   * @param entry The log entry to be inserted into the target log
   */
  void LogManager::AddLogEntry(string const& logname, string entry)
  {
    GetLogger(logname).AddLogEntry(entry);
  }

  /**
   * @brief LogManager::ExportLogs
   * Loops through all managed loggers, and calls their export method.
   */
  void LogManager::ExportLogs()
  {
    cout << "-- LogManager:: Exporting Logs"<<endl;
    LogMap::iterator iLog = GetLM()->fMap.begin();
    for (; iLog != GetLM()->fMap.end(); iLog++)
      (*iLog).second.Export();
    cout << "-- LogManager:: Log Export Completed"<<endl;

  }

  /**
   * @brief Logger::Logger
   * Constructor for an individual log class.
   * Each instance of Logger is associated to a filename
   * @param logname The identifier for the new log, e.g NNFIT or GAMinimizer
   * @param filename The filename of the new log, given as an absolute or relative path to the running directory.
   */
  Logger::Logger(string const& logname, string const& filename):
  fLogname(logname),
  fFilename(filename),
  fLogEntries()
  {
    cout << "** New Log File Generated. Log \'"<<logname<<"\' at "<<filename<<endl;
    return;
  }

  /**
   * @brief Logger::Export
   * Write the contents of the Logger class to it's target filename.
   */
  void Logger::Export()
  {
    if (fLogEntries.size() == 0)
    {
      cout << "** Log "<<fLogname<<" contains no entries. Log not exported. "<<endl;
      return;
    }

    std::stringstream datastream;
    datastream << "NNPDFCPP Log file. Identifier: "<<fLogname<<". "<<fLogEntries.size()<<" entries. "<<endl;

    for (size_t i=0; i<fLogEntries.size(); i++)
      datastream << fLogEntries[i]<<endl;

    write_to_file(fFilename, datastream.str());
    cout << "** Log "<<fLogname<<" successfully exported to "<<fFilename<<endl;
  }

}
