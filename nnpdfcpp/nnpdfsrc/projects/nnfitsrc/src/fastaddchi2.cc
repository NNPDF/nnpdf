// $Id$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include "fastaddchi2.h"
#include "datautils.h"
#include <NNPDF/thpredictions.h>
using NNPDF::ThPredictions;

void Convolute(const PDFSet* pdf, const Experiment* exp, real * theory)
{
  int index = 0;
  for (int i = 0; i < exp->GetNSet(); i++)
    {
      ThPredictions::Convolute(pdf,&exp->GetSet(i),theory+index);
      index += pdf->GetMembers()*exp->GetSet(i).GetNData();
    }
}

void FastAddChi2(const PDFSet* pdf, const DataSet* set, real* chi2)
{
  // Set up theory array
  const int nMem = pdf->GetMembers();
  real* theory = new real[set->GetNData()*nMem];

  // Perform convolution and chi^2 calculation
  ThPredictions::Convolute(pdf,set,theory);

  // Compute chi2
  ComputeChi2(set,nMem,theory,chi2);

  delete[] theory;
}

void FastAddChi2(const PDFSet* pdf, const Experiment* exp, real* chi2)
{
  // Set up theory array
  const int nMem = pdf->GetMembers();
  real *theory = new real[exp->GetNData()*nMem];

  // Perform convolution and chi^2 calculation
  Convolute(pdf,exp,theory);
  ComputeChi2(exp,nMem,theory,chi2);

  delete[] theory;
}
