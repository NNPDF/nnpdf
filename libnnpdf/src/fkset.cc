// $Id$
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <sstream>

#include "NNPDF/fkset.h"
#include "NNPDF/fastkernel.h"
#include "NNPDF/utils.h"
#include "NNPDF/exceptions.h"

namespace NNPDF
{

  /**
    * DataSet Operators
    * Ratio -> for W asymmetries etc
    */
  static void OpNull(int const& Nvals, std::vector<real*> const& obs, real*out)
  {
    for (int i=0; i<Nvals; i++)
      out[i] = obs[0][i];
    
    return;
  }

  static void OpAdd(int const& Nvals, std::vector<real*> const& obs, real*out)
  {
    for (int i=0; i<Nvals; i++)
      out[i] = obs[0][i] + obs[1][i];
    
    return;
  }

  static void OpRatio(int const& Nvals, std::vector<real*> const& obs, real*out)
  {
    if (obs.size()!=2)
      throw LengthError("OpRatio","number of FK grids is incorrect");
    
    for (int i=0; i<Nvals; i++)
      out[i] = obs[0][i]/obs[1][i];
    
    return;
  }

  static void OpAsy(int const& Nvals, std::vector<real*> const& obs, real*out)
  {
    if (obs.size()!=2)
      throw LengthError("OpAsy","number of FK grids is incorrect");
    
    for (int i=0; i<Nvals; i++)
      out[i] = (obs[0][i]-obs[1][i])/(obs[0][i]+obs[1][i]);
    
    return;
  }


  // Normalised sum operation
  static void OpSmn(int const& Nvals, std::vector<real*> const& obs, real*out)
  {
    if (obs.size()!=4)
     throw LengthError("OpSmn","number of FK grids is incorrect");

    for (int i=0; i<Nvals; i++)
      out[i] = (obs[0][i]+obs[1][i])/(obs[2][i]+obs[3][i]);
    
    return;
  }


  // FKSet
  FKSet::FKSet(SigmaOp op, std::vector<FKTable*> const& fktabs):
  fOperator(op),
  fNSigma(fktabs.size()),
  fNDataFK(fktabs[0]->GetNData()),
  fHadronic(fktabs[0]->IsHadronic()),
  fDataName(fktabs[0]->GetDataName()),
  fFK(new FKTable*[fNSigma])
  {
    // Copy FKTables
    for (size_t i=0; i<fktabs.size(); i++)
      fFK[i] = fktabs[i];

    if (fNSigma == 0)
      throw UserError("FKSet::FKSet","No FK tables added to set");

    for (int i=0; i<fNSigma; i++)
    {
      if (fFK[i]->IsHadronic() != fHadronic)
        throw EvaluationError("FKSet::FKSet", "Hadronic status mismatch!");

      if (fFK[i]->GetNData() != fNDataFK)
        throw EvaluationError("FKSet::FKSet", "NData mismatch!");

      if (fFK[i]->GetDataName().compare(fDataName) != 0)
        throw EvaluationError("FKSet::FKSet", "Setname mismatch: " + fFK[i]->GetDataName() + "  " + fDataName);
    }
  };


  // FKSet copy-constructor
  FKSet::FKSet(FKSet const& set):
  fOperator(set.fOperator),
  fNSigma(set.fNSigma),
  fNDataFK(set.fNDataFK),
  fHadronic(set.fHadronic),
  fDataName(set.fDataName),
  fFK(new FKTable*[fNSigma])
  {
    // Copy FKTables
    for (int i=0; i<set.fNSigma; i++)
      fFK[i] = new FKTable(*set.fFK[i]);

    // Verify masking (uneccesary after one go)
    for (int i=0; i<fNSigma; i++)
      if (fFK[i]->GetNData() != fNDataFK)
        throw RangeError("FKSet::FKSet", "NData mismatch!");
  };

  void swap(FKSet & lhs, FKSet & rhs)
  {
    using std::swap;
    swap(lhs.fOperator, rhs.fOperator);
    swap(lhs.fNSigma, rhs.fNSigma);
    swap(lhs.fNDataFK, rhs.fNDataFK);
    swap(lhs.fHadronic, rhs.fHadronic);
    swap(lhs.fDataName, rhs.fDataName);
    swap(lhs.fFK, rhs.fFK);
  }

  FKSet & FKSet::operator =(FKSet other){
    using std::swap;
    swap(*this, other);
    return *this;
  }


  FKSet::FKSet(FKSet && other):
  fOperator(OpNull),
  fNSigma(0),
  fNDataFK(0),
  fHadronic(false),
  fDataName(std::string()),
  fFK(nullptr)
  {
    using std::swap;
    swap(*this, other);
  }


  // FKSet masked copy-constructor
  FKSet::FKSet(FKSet const& set, std::vector<int> const& mask):
  fOperator(set.fOperator),
  fNSigma(set.fNSigma),
  fNDataFK( mask.size() ),
  fHadronic(set.fHadronic),
  fDataName(set.fDataName),
  fFK(new FKTable*[fNSigma])
  {
    // Copy FKTables
    for (int i=0; i<set.fNSigma; i++)
      fFK[i] = new FKTable(*set.fFK[i], mask);

    // Verify masking (uneccesary after one go)
    for (int i=0; i<fNSigma; i++)
      if (fFK[i]->GetNData() != fNDataFK)
        throw RangeError("FKSet::FKSet","NData mismatch!");
  };

  FKSet::~FKSet()
  {
    for (int i=0; i<fNSigma; i++)
      delete fFK[i];
    delete[] fFK;
  }

    // Parse dataset operators
  SigmaOp FKSet::parseOperator(std::string const& op)
  {

    if (op.compare("RATIO") == 0)
      return OpRatio;

    if (op.compare("ASY") == 0)
      return OpAsy;

    if (op.compare("ADD") == 0)
      return OpAdd;

    if (op.compare("SMN") == 0)
      return OpSmn;

    if (op.compare("NULL") == 0)
      return OpNull;

      // Add other operations here if required

    throw EvaluationError("FKSet::parseOperator","Unknown operator!");

    return OpNull;
  }

}

