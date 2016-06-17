// $Id: minimizer.h 1286 2013-10-28 11:54:20Z s0673800 $
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#pragma once

#include <cstdio>
#include <cstdlib>
#include <vector>

#include "common.h"
#include "fitpdfset.h"
using std::vector;

#include <NNPDF/experiments.h>
#include <NNPDF/positivity.h>
using NNPDF::Experiment;
using NNPDF::PositivitySet;

class NNPDFSettings;

/**
 *  \class Minimizer
 *  \brief Virtual minimisation base class
 */
class Minimizer
{
public:
  Minimizer(NNPDFSettings const&);
  virtual ~Minimizer();
  
  enum Mode {
    SetMode,
    ExpMode
  };  //!< Minimisation mode
  
  enum SortPDF {
    PDF_SORT,
    PDF_NOSORT
  };  //!< PDF sorting
  
  virtual void Init(FitPDFSet*, vector<Experiment*> const&, vector<PositivitySet> const&);
  virtual void Iterate(FitPDFSet*, vector<Experiment*> const&, vector<PositivitySet> const&) = 0;  //!< Perform an iteration of the minimisation
    
protected:
  virtual void ComputeErf(FitPDFSet*,
                  vector<Experiment*> const&,
                  vector<PositivitySet> const&,
                  Minimizer::Mode,
                  Minimizer::SortPDF); //!< Compute fChi2Mem
  
  real*  fChi2Mem;
  
  const NNPDFSettings& fSettings;
};

/**
 *  \class GAMinimizer
 *  \brief Basic Single Epoch Genetic Algorithm Minimizer
 */

class GAMinimizer : public Minimizer
{
public:
  GAMinimizer(NNPDFSettings const&);
  virtual ~GAMinimizer();
  
  virtual void Iterate(FitPDFSet*, vector<Experiment*> const&, vector<PositivitySet> const&);
  
protected:
  virtual void Mutation(FitPDFSet*, int const& nmut);
  int  Selection(FitPDFSet*);
};

/**
 *  \class NGAMinimiser
 *  \brief GA minimiser with nodal mutations
 */
 
class NGAMinimizer : public GAMinimizer
{
public:
  NGAMinimizer(NNPDFSettings const&);
  virtual ~NGAMinimizer();
  
protected:
   virtual void Mutation(FitPDFSet*, int const& nmut);
}; 

/*!
 * \brief The NGAPMinimizer class with preprocessing mutation
 */
class NGAPMinimizer : public NGAMinimizer
{
public:
  NGAPMinimizer(NNPDFSettings const&);

protected:
   void Mutation(FitPDFSet*, int const& nmut);

private:
   vector<real> falphamin;
   vector<real> falphamax;
   vector<real> fbetamin;
   vector<real> fbetamax;
};
