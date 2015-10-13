// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "common.h"

#include <string>
#include <sstream>
#include <map>
#include <vector>

namespace NNPDF
{

  class Logger;

  /**
   * \class LogManager
   * \brief Container class for individual loggers
   */
  class LogManager
  {
  public:
     static void    InitPath(std::string const& path) {GetLM()->fBasePath = path;}; //!< Initialise the Manager's base path
    
    static void    AddLogger(std::string const&, std::string const&);                   //!< Add a new log file to the manager
    static Logger& GetLogger(std::string const&);                                  //!< Get a log from the manager
    
    static void    AddLogEntry(std::string const& log, std::string ent); //!< Add entry 'ent' to log ID 'log'
    
    static void    ExportLogs();    //!< Export all active loggers
      
  private:
    LogManager();
    ~LogManager();
    
    // LogManager get/init method
    static LogManager* GetLM(){
      if (logInstance == 0)
        logInstance = new LogManager();
      
      return logInstance;
    };
    
    // Attributes
    static LogManager* logInstance; //!< Singleton log instance

    typedef std::map<size_t, Logger> LogMap;   //!< LogMap typedef
    LogMap fMap;                          //!< hash map of available loggers
    std::string fBasePath;                     //!< Base path of logger
    
  };

  /**
   * \class Logger
   * \brief Basic logging class. All methods private - can only be handled by logmanager
   */
  class Logger
  {
    Logger(std::string const&, std::string const&);
    
    void AddLogEntry(std::string log) {fLogEntries.push_back(log);};
    void Export();
    
    const std::string fLogname;
    const std::string fFilename;
    
    std::vector<std::string> fLogEntries;
    
    friend class LogManager;
  };

}