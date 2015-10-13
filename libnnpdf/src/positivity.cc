// $Id: positivity.cc 3187 2015-08-23 11:08:42Z stefano.carrazza@mi.infn.it $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#include "NNPDF/positivity.h"

#include "NNPDF/utils.h"
#include "NNPDF/fastkernel.h"
#include "NNPDF/pdfset.h"
#include "NNPDF/thpredictions.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>

namespace NNPDF
{

/**
 * @brief Positivity Constraint Set
 * @param The base commondata structure
 * @param The FKTable corresponding to the commondata
 * @param The Lagrange multiplier
 */
PositivitySet::PositivitySet(CommonData const& data, FKTable const& table, real const& lambda):
CommonData(data),
FKTable(table),
fLambda(lambda)
{
  std::cout << "Positivity Lagrange Multiplier: "<<fLambda<<std::endl;
}

/**
 * @brief Positivity Constraint Set
 * @param The base commondata structure
 * @param The FKTable corresponding to the commondata
 * @param The Lagrange multiplier
 */
PositivitySet::PositivitySet(PositivitySet const& set):
CommonData(set),
FKTable(set),
fLambda(set.fLambda),
fBounds(set.fBounds)
{
}

/**
 * @brief The positivity set destructor
 */
PositivitySet::~PositivitySet()
{
  fBounds.clear();
}

/**
  *
  */
void PositivitySet::SetBounds(const NNPDF::PDFSet* pdf)
{
  if (pdf->GetMembers() > 1) { std::cerr << "Positivity::SetBounds error, bound PDF contains more than 1 replica" << std::endl; exit(-1); }

  const int Ndat = CommonData::GetNData();
  real *tmp = new real[Ndat];
  ComputePoints(pdf,tmp);

  fBounds.clear();
  for (int i = 0; i < Ndat; i++)
    fBounds.push_back(-0.25*fabs(tmp[i]));
  
  delete[] tmp;
}

/**
 * @brief Compute contribution to error function from Positivity violation
 * @param pdf The input PDF set which will be used to compute the error function
 * @param res The error function for each pdf member (res[i])
 * if the theoretical prediction for such observable is < 0 we penalize the error function
 */
void PositivitySet::ComputeErf(const NNPDF::PDFSet* pdf, real* res) const
{  
  const int Npdf = pdf->GetMembers();
  const int Ndat = CommonData::GetNData();
  
  // Compute positivity points
  real* tmp = new real[Npdf*Ndat];
  ComputePoints(pdf,tmp);
  
  // Contribution to Error Function
  for (int i = 0; i< Ndat; i++)
    for (int j = 0; j< Npdf; j++)
      if (tmp[i*Npdf+j] < 0)
        res[j] -= fLambda*tmp[i*Npdf+j];

  delete[] tmp;
}

/**
 * @brief Computes the theoretical predictions for Positivity observables
 * @param pdf The PDF set used in the convolution
 * @param res The observable, in vector of members dimention res[i]
 */
void PositivitySet::ComputePoints(const PDFSet* pdf, real* res) const
{
  for (int i=0; i<CommonData::GetNData()*pdf->GetMembers(); i++)
    res[i] = 0;
  
  ThPredictions::Convolute(pdf,this,res);
}


/**
 * @brief Computes the total number of violated data points (negative observables)
 * @param pdf The PDF set used to compute the theoretical predictions
 * @param res The total number of points which violates the positivity observable per PDF member
 */
void PositivitySet::ComputeNViolated( const PDFSet* pdf, int* res) const
{  
  const int Npdf = pdf->GetMembers();
  const int Ndat = CommonData::GetNData();
  
  // Compute positivity points
  real* tmp = new real[Npdf*Ndat];
  ComputePoints(pdf,tmp);

  // Contribution to Error Function
  for (int j = 0; j< pdf->GetMembers(); j++)
  {
    res[j] = 0;
    for (int i = 0; i< Ndat; i++)
        if (tmp[i*Npdf+j] < 0)
          res[j] += 1;
  }

  delete[] tmp;
}


/**
 * @brief Computes the total number of unacceptable data points (observables less than bounds). Ignores first and last points.
 * @param pdf The PDF set used to compute the theoretical predictions
 * @param pdf The PDF set used to compute the theoretical predictions
 * @param res The total number of points which are unacceptable per PDF member
 */
void PositivitySet::ComputeNUnacceptable( const PDFSet* pdf, int* res) const
{  
  const int Npdf = pdf->GetMembers();
  const int Ndat = CommonData::GetNData();
  
  // Compute positivity points
  real* tmp = new real[Npdf*Ndat];
  ComputePoints(pdf,tmp);

  // Contribution to Error Function
  for (int j = 0; j< pdf->GetMembers(); j++)
  {
    res[j] = 0;
    for (int i = 0; i< Ndat; i++)
        if (tmp[i*Npdf+j] < fBounds[i])
          res[j] += 1;
  }

  delete[] tmp;
}

}
