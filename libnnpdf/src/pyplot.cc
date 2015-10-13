// $Id
// NNPDF++ 2014
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk
//          Chris Deans,      C.S.Deans@sms.ed.ac.uk

#include <sstream>

#include "NNPDF/pyplot.h"

using namespace NNPDF;

// Singleton management
PyPlot *PyPlot::plotInstance = 0;

PyPlot::PyPlot()
{
  InitInterpreter();
}

PyPlot::~PyPlot()
{
  CloseInterpreter();
}

void PyPlot::InitInterpreter()
{
  if (Py_IsInitialized())
  {
    std::cerr << "PyPlot::InitInterpreter::Error - Interpreter is already initialised!"<<std::endl;
    exit(-1);
  }
  
  // Init interpreter
  Py_Initialize();
  
  // Set name
  Py_SetProgramName("NNPDF++ PyPlot");
  
  // Add scriptDir to pythonpath
  std::stringstream PathDir;
  PathDir << "sys.path.append(os.getcwd() + '/"<<scriptPath()<<"')";
  
  PyRun_SimpleString("import sys, os");
  PyRun_SimpleString(PathDir.str().c_str());
}

void PyPlot::CloseInterpreter()
{
  if (!Py_IsInitialized())
  {
    std::cerr << "PyPlot::CloseInterpreter::Error - Interpreter is not initialised!"<<std::endl;
    exit(-1);
  }
  
  Py_Finalize();
}

void PyPlot::ListPlots()
{
  PyPlot::GetPlotter();
  // Plot information
  /*
   PyObject *pPlotInfo = PyDict_GetItemString(pDict, "info");
   
   if (PyCallable_Check(pPlotInfo))
   {
   PyObject_CallObject(pPlotInfo, NULL);
   } else
   {
   PyErr_Print();
   }
   */
  
  // need to get some directory iterating structure here
}

void PyPlot::ExecPlot(std::string pltName, std::string srcPath, std::string trgPath)
{
  
  //InitInterpreter();
  
  PyObject *pPath = PyString_FromString(pltName.c_str());
  
  // Load python module and dictionary
  PyObject *pMod = PyImport_Import(pPath);
  
  if (!pMod)
  {
    std::cerr << "PyPlot::ExecPlot Fails: Cannot import plot "<<pltName<<std::endl;
    exit(-1);
  }
  
  // Fetch dictionary reference
  PyObject *pDict = PyModule_GetDict(pMod);
  
  // Plot args
  PyObject *pArgs = PyTuple_New(2);
  PyObject *pSource = PyString_FromString(srcPath.c_str());
  PyObject *pTarget = PyString_FromString(trgPath.c_str());

  if (!pSource || !pTarget ) {
    std::cerr << "PyPlot::ExecPlot::Cannot convert argument "<<std::endl;
    exit(-1);
  }
  
  // Set args
  PyTuple_SetItem(pArgs, 0, pSource);
  PyTuple_SetItem(pArgs, 1, pTarget);
  
  // Run Plot
  PyObject *pPlotRun = PyDict_GetItemString(pDict, "plot");
  
  if (PyCallable_Check(pPlotRun))
  {
    PyObject_CallObject(pPlotRun, pArgs);
  } else
  {
    PyErr_Print();
    exit(-1);
  }
  
  //CloseInterpreter();

}