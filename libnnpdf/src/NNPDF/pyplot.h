// $Id
// NNPDF++ 2014
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk
//          Chris Deans,      C.S.Deans@sms.ed.ac.uk

#pragma once

#include "common.h"
#include <Python.h>
#include <map>

namespace NNPDF
{

  // Python plotting routines manager
  class PyPlot
  {
  public:
    
    static PyPlot* GetPlotter()
    {
      if (plotInstance == 0)
      plotInstance = new PyPlot();
      
      return plotInstance;
    }

   
    void ListPlots();
    void ExecPlot(std::string pltName, std::string srcPath, std::string trgDir);
    void RestartInt() {CloseInterpreter(); InitInterpreter();}
    
  private:
    PyPlot();
    ~PyPlot();
    
    void InitInterpreter();
    void CloseInterpreter();
    
    static PyPlot* plotInstance;
    std::map<std::string,PyObject*> fModuleMap;
  };

}