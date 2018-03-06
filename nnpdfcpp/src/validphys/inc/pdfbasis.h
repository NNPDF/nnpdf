// $Id$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#pragma once

#include "common.h"
#include <vector>
using std::vector;

class PDFSet;
class NNPDFSettings;

/**
 *  \class PDFBasis
 *  \brief Base class for the storage of preprocessing factors
 */
class PreprocParam
{
public:
  PreprocParam(const int npdf):
  fNPDF(npdf),
  fPDFNorm(new real[fNPDF]),
  fPDFAux(new real[fNPDF])
  {
    for (int i=0; i<fNPDF; i++)
    {
      fPDFNorm[i] = 1.0;
      fPDFAux[i] = 0.0;
    }
  }

  PreprocParam( PreprocParam const& copy):
  fNPDF(copy.fNPDF),
  fPDFNorm(new real[fNPDF]),
  fPDFAux(new real[fNPDF])
  {
    for (int i=0; i<fNPDF; i++)
    {
      fPDFNorm[i] = copy.fPDFNorm[i];
      fPDFAux[i] = copy.fPDFAux[i];
    }

  }

  ~PreprocParam()
  {
    delete[] fPDFNorm;
    delete[] fPDFAux;
  }

  const int  fNPDF;     //!< Number of members
  real* const fPDFNorm; //!< PDF normalisations
  real* const fPDFAux;  //!< Normalisation of auxilliary terms
};

/**
 *  \class PDFBasis
 *  \brief Base class for all PDF basis definitions
 */
class PDFBasis
{
public:
  PDFBasis(string name, int nPDF):
  fBasisName(name),
  fPDFNames(new string[nPDF]),
  fNPDF(nPDF)
  {
    cout << "PDFBasis:: initialised basis: "<<fBasisName<<endl;
    return;
  }

  virtual ~PDFBasis()
  {
    delete[] fPDFNames;
  }

  const int& GetNPDF() const {return fNPDF;} //!< Return number of PDFs in basis
  const string& GetPDFName(int const& pdf){ return fPDFNames[pdf]; }  //!< Return name of PDF in basis

  // Common transformation methods for all bases
  void LHA2EVLN(const real *LHA,  real *EVLN) const;
  void EVLN2LHA(const real *EVLN, real *LHA ) const;

  // Utility method for basis->LHA
  void BASIS2LHA(real const* basis, real* lha) const;
  void LHA2BASIS(real const* lha, real* basis) const;

  // Basis to EVLN and back
  virtual void BASIS2EVLN(real const* basis, real* evln) const = 0;
  virtual void EVLN2BASIS(real const* evln, real* basis) const = 0;

  // Compute associated sum rules
  virtual real ComputeSumRules(sumRule, int mem, PDFSet*, bool&) const = 0;

protected:
  const string fBasisName;    //!< Name of the basis
  string* const fPDFNames;    //!< Names of the PDFs in the basis (for plotting)

  const int  fNPDF;           //!< Number of PDFs in the basis
};

/**
 *  \class LHABasis
 *  \brief nFL flavour basis
 */
class LHABasis: public PDFBasis
{
public:
  LHABasis(): PDFBasis("LHABasis", 13) {}

  void BASIS2EVLN(real const* basis, real* evln) const;
  void EVLN2BASIS(real const* evln, real* basis) const;

  real ComputeSumRules(sumRule,int mem, PDFSet*, bool&) const;
};


