// $Id$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#pragma once

#include "common.h"
#include "pdfbasis.h"

namespace NNPDF {
  class LHAPDFSet;
}

/**
 * Return an initialised fitting basis
 */
class FitBasis;
FitBasis* getFitBasis(NNPDFSettings const& settings, basisType btype, int const& rep = 0);

/**
 *  \class FitBasis
 *  \brief Base class for all Fitting PDF basis definitions
 */
class FitBasis: public PDFBasis
{
public:
  FitBasis(NNPDFSettings const&, string name, int nPDF); //!< FitBasis constructor

  virtual void BASIS2EVLN(real const* basis, real* evln) const = 0;
  virtual void EVLN2BASIS(real const* evln, real* basis) const = 0;

  // Compute associated sum rules
  virtual real ComputeSumRules(sumRule, int mem, PDFSet*, bool&) const = 0;

  // Preprocessing
  virtual void ComputeParam(PDFSet*, int mem, PreprocParam&, bool&) const = 0;
  virtual void Preprocess(real const& x, int const& fl, real& pdf, PreprocParam const&);
  virtual void Preprocess(real const& x, real* pdf, PreprocParam const&);
  virtual void NetTransform(int const& fl, int const& nfl, int*);

  real const& GetAlpha(int const& fl) const { return fAlpha[fl]; }
  real const& GetBeta(int const& fl) const { return fBeta[fl]; }
  void SetAlpha(int const& fl, real const& v) { fAlpha[fl] = v; }
  void SetBeta(int const& fl, real const& v)  { fBeta[fl] = v; }

  bool const& GetPDFSqrPos(int const& fl) const { return fPDFSqrPos[fl]; }

  double* fArcDampFactor;
protected:
  // Basic Preprocessing Constants
  bool* const fPDFSqrPos; //!< Is the PDF to be squared?
  real* const fAlpha;       //!< Low-x Preprocessing exponents
  real* const fBeta;        //!< High-x Preprocessing exponents
  real fQ2;                 //!< Fit initial scale for integration
  gsl_integration_workspace* fGSLWork; //!< workspace for integration
};


/**
 *  \class NN23FitBasis
 *  \brief Fit basis used in NNPDF releases 2.3 and below
 */
class NN23FitBasis: public FitBasis
{
public:
  NN23FitBasis(NNPDFSettings const&);

  enum fitBasis {FIT_SNG, FIT_GLU, FIT_VAL, FIT_T3, FIT_DS, FIT_SP, FIT_SM, FIT_GAM };

  void BASIS2EVLN(real const* basis, real* evln) const;
  void EVLN2BASIS(real const* evln, real* basis) const;

  real ComputeSumRules(sumRule, int mem, PDFSet*, bool&) const;

  // Preprocessing
  void ComputeParam(PDFSet*, int mem, PreprocParam&, bool&) const;
  void Preprocess(real const& x, int const& fl, real& pdf, PreprocParam const&);

private:
  const bool fQED;

  // Strange Auxilliary terms
  real fSauxAlpha;
  real fSauxBeta;
  real fSauxGamma;
};

/**
 *  \class EvolFitBasis
 *  \brief Evol basis for fitting
 */
class EvolFitBasis: public FitBasis
{
public:
  EvolFitBasis(NNPDFSettings const&);

  // Σ, g, V, V3, V8, V15, V24, V35, T3, T8, T15, T24, T35, γ

  enum fitBasis {FIT_SNG, FIT_GLU, FIT_VAL, FIT_V3, FIT_V8, FIT_T3, FIT_T8, FIT_GAM };

  void BASIS2EVLN(real const* basis, real* evln) const;
  void EVLN2BASIS(real const* evln, real* basis) const;

  real ComputeSumRules(sumRule, int mem, PDFSet*, bool&) const;

  // Preprocessing
  void ComputeParam(PDFSet*, int mem, PreprocParam&, bool&) const;

protected:
  bool fQED;
};

/**
 * @brief The LuxBasis class
 */
class LuxBasis: public FitBasis
{
public:
  LuxBasis(NNPDFSettings const&set, int const& replica);
  ~LuxBasis();

  enum fitBasis { FIT_SNG, FIT_GLU, FIT_VAL, FIT_V3, FIT_V8, FIT_T3, FIT_T8, FIT_CP, FIT_GAM };

  void BASIS2EVLN(real const* basis, real* evln) const;
  void EVLN2BASIS(real const* evln, real* basis) const;

  real ComputeSumRules(sumRule, int mem, PDFSet*, bool&) const;

  // Preprocessing
  void ComputeParam(PDFSet*, int mem, PreprocParam&, bool&) const;

  void Preprocess(real const& x, int const& fl, real& pdf, PreprocParam const&);
  void Preprocess(real const& x, real* pdf, PreprocParam const& par);

private:
  double fQ0;
  NNPDF::LHAPDFSet* fPhotonSet;
};


/**
 *  \class EvolSFitBasis
 *  \brief Evol Fit basis for strangeness, otherwise NN23
 */
class EvolSFitBasis: public FitBasis
{
public:
  EvolSFitBasis(NNPDFSettings const&);

  enum fitBasis {FIT_SNG, FIT_GLU, FIT_VAL, FIT_V8, FIT_T3, FIT_T8, FIT_DS, FIT_GAM };

  void BASIS2EVLN(real const* basis, real* evln) const;
  void EVLN2BASIS(real const* evln, real* basis) const;

  real ComputeSumRules(sumRule, int mem, PDFSet*, bool&) const;

  // Preprocessing
  void ComputeParam(PDFSet*, int mem, PreprocParam&, bool&) const;

private:
  const bool fQED;

};

class NN30FitBasis: public EvolFitBasis
{
public:
   NN30FitBasis(NNPDFSettings const&);

   enum netBasis {NET_SNG, NET_GLU, NET_VAL, NET_T3, NET_DS, NET_SP, NET_SM, NET_GAM };

   void NetTransform(int const& fl, int const& nfl, int*);
   void Preprocess(real const& x, real* pdf, PreprocParam const&);
};

class FLVRFitBasis: public EvolFitBasis
{
public:
   FLVRFitBasis(NNPDFSettings const&);

   enum netBasis {NET_GLU, NET_U, NET_UBAR, NET_D, NET_DBAR, NET_S, NET_SBAR, NET_GAM };

   void NetTransform(int const& fl, int const& nfl, int*);
   void Preprocess(real const& x, real* pdf, PreprocParam const&);
};

/**
 *  \class EvolICFitBasis
 *  \brief Evol basis for fitting intrinsic charm
 */
class EvolICFitBasis: public FitBasis
{
public:
  EvolICFitBasis(NNPDFSettings const&);

  // γ, Σ, g, V, V3, V8, V15, V24, V35, T3, T8, T15, T24, T35

  enum fitBasis {FIT_SNG, FIT_GLU, FIT_VAL, FIT_V3, FIT_V8, FIT_T3, FIT_T8, FIT_T15, FIT_GAM };

  void BASIS2EVLN(real const* basis, real* evln) const;
  void EVLN2BASIS(real const* evln, real* basis) const;

  real ComputeSumRules(sumRule, int mem, PDFSet*, bool&) const;

  // Preprocessing
  void ComputeParam(PDFSet*, int mem, PreprocParam&, bool&) const;

protected:
  const bool fQED;

};

/**
 * @brief The NN30ICFitBasis class
 */
class NN30ICFitBasis: public EvolICFitBasis
{
public:
   NN30ICFitBasis(NNPDFSettings const&);

   //enum netBasis {NET_SNG, NET_GLU, NET_VAL, NET_T3, NET_DS, NET_SP, NET_SM, NET_CP, NET_CM, NET_GAM };
   enum netBasis {NET_SNG, NET_GLU, NET_VAL, NET_T3, NET_DS, NET_SP, NET_SM, NET_CP, NET_GAM }; // c- = 0

   void NetTransform(int const& fl, int const& nfl, int*);
   void Preprocess(real const& x, real* pdf, PreprocParam const&);
};

/**
 *  \class NN31ICFitBasis
 *  \brief Evol basis for fitting intrinsic charm
 */
class NN31ICFitBasis: public FitBasis
{
public:
  NN31ICFitBasis(NNPDFSettings const&);

  // Σ, g, V, V3, V8, T3, T8, c+, (γ)

  enum fitBasis {FIT_SNG, FIT_GLU, FIT_VAL, FIT_V3, FIT_V8, FIT_T3, FIT_T8, FIT_CP, FIT_GAM };

  void BASIS2EVLN(real const* basis, real* evln) const;
  void EVLN2BASIS(real const* evln, real* basis) const;

  real ComputeSumRules(sumRule, int mem, PDFSet*, bool&) const;

  // Preprocessing
  void ComputeParam(PDFSet*, int mem, PreprocParam&, bool&) const;

protected:
  const bool fQED;

};

/**
 *  \class NoSumRuleBasis
 *  \brief Test basis, identical to NN31ICFitBasis but with no sum rules applied
 */
class NoSumRuleBasis: public FitBasis
{
public:
  NoSumRuleBasis(NNPDFSettings const&);

  // Σ, g, V, V3, V8, T3, T8, c+, (γ)

  enum fitBasis {FIT_SNG, FIT_GLU, FIT_VAL, FIT_V3, FIT_V8, FIT_T3, FIT_T8, FIT_CP, FIT_GAM };

  void BASIS2EVLN(real const* basis, real* evln) const;
  void EVLN2BASIS(real const* evln, real* basis) const;

  real ComputeSumRules(sumRule, int mem, PDFSet*, bool&) const;

  // Preprocessing
  void ComputeParam(PDFSet*, int mem, PreprocParam&, bool&) const;

protected:
  const bool fQED;

};


