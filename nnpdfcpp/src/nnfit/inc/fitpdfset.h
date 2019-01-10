// $Id: fitpdfset.h 1506 2014-01-15 11:26:09Z s0673800 $
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#pragma once

#include <cstdio>
#include <cstdlib>

#include "common.h"
#include "fitbases.h"
#include <NNPDF/pdfset.h>
#include <NNPDF/parametrisation.h>
#include <nnpdfsettings.h>
#include <gsl/gsl_integration.h>

using std::vector;
using std::pair;
using NNPDF::PDFSet;
using NNPDF::Parametrisation;
class PreprocParam;

/**
 *  \class NNPDFSet
 *  \brief Neural Network PDFSet to be minimized
 */

class FitPDFSet : public PDFSet
{
public:
  ~FitPDFSet();

  template <class T, class A> static FitPDFSet* Generate(NNPDFSettings const& settings, FitBasis* basis)
  {
    const int nfl = basis->GetNPDF();
    FitPDFSet* ns = new FitPDFSet(settings, basis);

    Parametrisation** newpdf = new Parametrisation*[nfl];
    for (int i=0; i<nfl; i++)
    {
      T* np = new T(settings.GetArch());
      newpdf[i] = static_cast<Parametrisation*>(np);
    }

    ns->fBestFit=newpdf;
    return ns;
  }

  void InitPDFSet() const;
  void ExportMeta(int const& rep, real const& erf_val, real const& erf_trn, real const& chi2, bool posVeto);

  bool ComputeIntegrals( int const& mem ); //!< Compute all associated integrals and sum rules
  void ComputeSumRules();                  //!< Compute preprocessing sum rule constraints over all members
  void ValidateStartingPDFs();             //!< Validate initial PDFs

  void SetNMembers(int const& mem) {fMembers = mem; ExpandMembers();}

  vector<Parametrisation**>& GetPDFs() {return fPDFs;}
  Parametrisation**          GetBestFit() {return fBestFit;}
  void SortMembers(real*);

  void SetBestFit(int const&);
  real GetEbf()         { return fEbf; }
  void SetEbf( real const& e) {fEbf = e;} //!< Sets the new best fit error function

  int  GetNIte()        const { return fNIte; }
  void SetNIte( int const& newIte) { fNIte = newIte; }
  void Iterate()        { fNIte++; }

  void GetPDF (real const& x, real const& Q2, int const& n, real* pdf) const; //!< Get evolution basis PDF
  real GetPDF  (real const& x, real const& Q2, int const& n, int const& fl) const; //!< Get preprocessed Fit basis PDF

  real CalculateArcLength(int const& mem, int const& fl, real const& dampfact, real xmin = 1e-15, real xmax = 1.0) const;

  real GetQ20() const { return fQ20; }

private:
  FitPDFSet(NNPDFSettings const&, FitBasis*);

  void ExpandMembers();                      //!< Expand internal vectors to fMembers
  void DisableMember(int mem);               //!< Disable a member PDF by moving it to the end of the vector and decrementing fMembers

  NNPDFSettings const& fSettings;

  FitBasis* fFitBasis;                      //!< Fitting basis for PDF

  const int fNfl;
  const real fQ20;
  vector<PreprocParam*>     fPreprocParam; //!< PDF preprocessing parameters by member

  Parametrisation**         fBestFit;     //!< Best fit PDF
  vector<Parametrisation**> fPDFs;        //!< Vector of PDF members
  real  fEbf;                             //!< Figure of merit for best fit PDF
  int fNIte;                              //!< Counts the number of fit iterations
  basisType fbtype;                       //!< store the basis type

  friend class Minimizer;
};
