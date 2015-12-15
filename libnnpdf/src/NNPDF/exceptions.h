// $Id$
//
// NNPDF++ 2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#pragma once

#include <exception>
#include <stdexcept>

namespace NNPDF
{
  /**
   *
   * Simple and generic exception handlers
   *  - implements runtime and logic errors
   *
   */

  /// Error to be thrown when runtime error is detected
  class RuntimeException: public std::runtime_error
  {
  public:
    RuntimeException(const std::string& tag, const std::string& what) : std::runtime_error("[" + tag + "] error: " + what) {}
  };

  /// Error to be thrown when logic error is detected
  class LogicException: public std::logic_error
  {
  public:
    LogicException(const std::string& tag, const std::string& what) : std::logic_error("[" + tag + "] error: " + what) {}
  };

  //_______________________________________________________________________
  class FileError: public RuntimeException
  {
  public:
    FileError(const std::string& tag, const std::string& what) : RuntimeException(tag,what) {}
  };

  class EvaluationError: public RuntimeException
  {
  public:
    EvaluationError(const std::string& tag, const std::string& what) : RuntimeException(tag,what) {}
  };

  class InitError: public RuntimeException
  {
  public:
    InitError(const std::string& tag, const std::string& what) : RuntimeException(tag,what) {}
  };

  class RangeError: public LogicException
  {
  public:
    RangeError(const std::string& tag, const std::string& what) : LogicException(tag,what) {}
  };

  class LengthError: public LogicException
  {
  public:
    LengthError(const std::string& tag, const std::string& what) : LogicException(tag,what) {}
  };

  class LogError: public LogicException
  {
  public:
    LogError(const std::string& tag, const std::string& what) : LogicException(tag,what) {}
  };

  class UserError: public LogicException
  {
  public:
    UserError(const std::string& tag, const std::string& what) : LogicException(tag,what) {}
  };

}
