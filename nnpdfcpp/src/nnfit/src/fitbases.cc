// $Id$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include "fitbases.h"
#include "nnpdfsettings.h"

#include <NNPDF/pdfset.h>
#include <NNPDF/lhapdfset.h>
#include <gsl/gsl_sf_gamma.h>

// Initialise a fit basis
FitBasis* getFitBasis(NNPDFSettings const& settings, basisType btype, const int &rep)
{
  // Fit Basis
  FitBasis* fitbasis = NULL;
  switch ( btype ) {

    case BASIS_NN23:
    case BASIS_NN23QED:
    {
      fitbasis = new NN23FitBasis(settings);
      cout << Colour::FG_BLUE << "Selecting FitBasis: NN23" << Colour::FG_DEFAULT << endl;
      break;
    }

    case BASIS_EVOL:
    case BASIS_EVOLQED:
    {
      fitbasis = new EvolFitBasis(settings);
      cout << Colour::FG_BLUE << "Selecting FitBasis: EVOL" << Colour::FG_DEFAULT << endl;
      break;
    }
    case BASIS_DISEVOL:
    {
      fitbasis = new DISEvolFitBasis(settings);
      cout << Colour::FG_BLUE << "Selecting FitBasis: DISEVOL" << Colour::FG_DEFAULT << endl;
      break;
    }
    case BASIS_LUX:
    {
      fitbasis = new LuxBasis(settings, rep);
      cout << Colour::FG_BLUE << "Selecting FitBasis: LUX" << Colour::FG_DEFAULT << endl;
      break;
    }

    case BASIS_EVOLS:
    case BASIS_EVOLSQED:
    {
      fitbasis = new EvolSFitBasis(settings);
      cout << Colour::FG_BLUE <<"Selecting FitBasis: EVOLS" << Colour::FG_DEFAULT << endl;
      break;
    }

    case BASIS_NN30:
    case BASIS_NN30QED:
    {
      fitbasis = new NN30FitBasis(settings);
      cout << Colour::FG_BLUE << "Selecting FitBasis: NN30" << Colour::FG_DEFAULT << endl;
      break;
    }

    case BASIS_FLVR:
    case BASIS_FLVRQED:
    {
      fitbasis = new FLVRFitBasis(settings);
      cout << Colour::FG_BLUE << "Selecting FitBasis: FLVR" << Colour::FG_DEFAULT << endl;
      break;
    }

    case BASIS_NN30IC:
    {
      fitbasis = new NN30ICFitBasis(settings);
      cout << Colour::FG_BLUE << "Selecting FitBasis: NN30IC" << Colour::FG_DEFAULT << endl;
      break;
    }

    case BASIS_EVOLIC:
    {
      fitbasis = new EvolICFitBasis(settings);
      cout << Colour::FG_BLUE << "Selecting FitBasis: EVOLIC" << Colour::FG_DEFAULT << endl;
      break;
    }

    case BASIS_NN31IC:
    case BASIS_NN31ICQED:
    {
      fitbasis = new NN31ICFitBasis(settings);
      cout << Colour::FG_BLUE << "Selecting FitBasis: NN31IC" << Colour::FG_DEFAULT << endl;
      break;
    }

    case BASIS_NSR:
    {
      fitbasis = new NoSumRuleBasis(settings);
      cout << Colour::FG_BLUE << "Selecting FitBasis: NSR" << Colour::FG_DEFAULT << endl;
      break;
    }

    default:
      cerr << Colour::FG_RED << "[getFitBasis] error: Invalid Fitting Basis" << Colour::FG_DEFAULT << endl;
      exit(-1);
      break;
  }

  return fitbasis;

}

/**
 *  Basic FitBasis constructor
 *  Common attributes for all fit bases
 **/
FitBasis::FitBasis(NNPDFSettings const& nnset, string name, int nPDF):
PDFBasis(name, nPDF),
fArcDampFactor(new double[fNPDF]),
fPDFSqrPos(new bool[fNPDF]),
fAlpha(new real[fNPDF]),
fBeta(new real[fNPDF]),
fQ2( static_cast<real>(pow(stod(nnset.GetTheory(APFEL::kQ0)),2))),
fGSLWork(nnset.GetGSLWorkspace())
{
  // Squared Positivity
  for (int i = 0; i < fNPDF; i++)
    fPDFSqrPos[i] = nnset.Get("fitting","basis")[i]["pos"].as<bool>();

  // Preprocessing constants
  RandomGenerator* rg = RandomGenerator::GetRNG();

    for (int i = 0; i < fNPDF; i++)
      {
        // Small x exponents
        fAlpha[i] = -rg->GetRandomUniform(nnset.Get("fitting","basis")[i]["smallx"][0].as<real>(),
                                          nnset.Get("fitting","basis")[i]["smallx"][1].as<real>());

        // Large x exponents
        fBeta[i]  =  rg->GetRandomUniform(nnset.Get("fitting","basis")[i]["largex"][0].as<real>(),
                                          nnset.Get("fitting","basis")[i]["largex"][1].as<real>());
      }

  return;
}

/**
 * @brief The FitBasis destructor
 */
FitBasis::~FitBasis()
{
  delete[] fArcDampFactor;
  delete[] fPDFSqrPos;
  delete[] fAlpha;
  delete[] fBeta;
}

// Preprocess a supplied PDF
void FitBasis::Preprocess(real const& x, int const& fl, real& pdf, PreprocParam const& par)
{
  //Preprocessing and normalisation
  pdf *= par.fPDFNorm[fl]*pow(1-x,fBeta[fl])*pow(x,fAlpha[fl]+1);

  return;
}

void FitBasis::Preprocess(real const& x, real* pdf, PreprocParam const& par)
{
  for (int i = 0; i < fNPDF; i++)
  {
    // Positive definite PDFs
    if (fPDFSqrPos[i])
      pdf[i] *= pdf[i];

    Preprocess(x,i,pdf[i],par);
  }

  return;
}

void FitBasis::NetTransform(int const& fl, int const& nfl, int* transform)
{
  for (int i = 0; i < nfl; i++)
    transform[i] = 0;
  transform[fl] = 1;
}

/**
 *  NNPDF2.3 Fit Basis
 **/

NN23FitBasis::NN23FitBasis(NNPDFSettings const& nnset):
FitBasis(nnset, "NN23FitBasis", 7+nnset.IsQED()),
fQED(nnset.IsQED())
{
  // PDF Names for plotting
  fPDFNames[FIT_SNG] = "Singlet";
  fPDFNames[FIT_GLU] = "Gluon";
  fPDFNames[FIT_VAL] = "Valence";
  fPDFNames[FIT_T3] = "Triplet";
  fPDFNames[FIT_DS] = "Sea Asymmetry";
  fPDFNames[FIT_SP] = "Strange Sea";
  fPDFNames[FIT_SM] = "Strange Valence";
  if (fQED)
    fPDFNames[FIT_GAM] = "Photon";

  // Damping factor for arclengths
  fArcDampFactor[FIT_SNG] = 1;
  fArcDampFactor[FIT_GLU] = 1;
  fArcDampFactor[FIT_VAL] = 0;
  fArcDampFactor[FIT_T3] = 1;
  fArcDampFactor[FIT_DS] = 1;
  fArcDampFactor[FIT_SP] = 1;
  fArcDampFactor[FIT_SM] = 0;
  if (fQED)
    fArcDampFactor[FIT_GAM] = 1;

  RandomGenerator* rg = RandomGenerator::GetRNG();
  fSauxBeta = 3.5 + rg->GetRandomUniform(0.0,1.0);
  fSauxAlpha = fSauxBeta/2;
  fSauxGamma = gsl_sf_gamma(fSauxBeta + fSauxAlpha + 2)/(gsl_sf_gamma(fSauxAlpha+1)*gsl_sf_gamma(fSauxBeta+1));
}

void NN23FitBasis::ComputeParam(PDFSet* pdf, int mem, PreprocParam& param, bool &status) const
{
  // status
  status = false;

  // Clear old normalisations
  for (int i=0; i<fNPDF; i++)
  {
    param.fPDFNorm[i] = 1.0;
    param.fPDFAux[i] = 0.0;
  }

  // Normalisations pointer
  real* norm = param.fPDFNorm;

  // Quantum number sum rules
  norm[FIT_VAL] = 3.0f/pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork); // Total valence
  norm[ FIT_DS] = (1-pdf->IntegratePDF(mem,FIT_T3,fQ2, PDFSet::FX,status,fGSLWork))/
  (2*pdf->IntegratePDF(mem,FIT_DS,fQ2, PDFSet::FX,status,fGSLWork)); // D_S - 1-t3/2d_s

  // Strange valence sum rule
  param.fPDFAux[FIT_SM] = pdf->IntegratePDF(mem,FIT_SM,fQ2,PDFSet::FX,status,fGSLWork)*fSauxGamma;

  // ************ QED dependent normalisations **************
  if (fQED)
  {
    norm[FIT_GLU] = (1-pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork) - pdf->IntegratePDF(mem,FIT_GAM,fQ2,PDFSet::XFX,status,fGSLWork))
    / pdf->IntegratePDF(mem,FIT_GLU,fQ2,PDFSet::XFX,status,fGSLWork);
  }
  else
  {
    norm[FIT_GLU] = (1-pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork))/
    pdf->IntegratePDF(mem,FIT_GLU,fQ2,PDFSet::XFX,status,fGSLWork);
  }

  return;

}

// Preprocess a supplied PDF
void NN23FitBasis::Preprocess(real const& x, int const& fl, real& pdf, PreprocParam const& par)
{
  // Basic Preprocessing
  FitBasis::Preprocess(x,fl,pdf,par);

  // Strange auxilliary term
  if (fl == FIT_SM)
    pdf -= par.fPDFAux[FIT_SM]*pow(1-x,fSauxBeta)*pow(x,fSauxAlpha+1);

  return;
}

/**
 * @brief NN23FitBasis::BASIS2EVLN
 * @param FIT
 * @param EVLN
 */
void NN23FitBasis::BASIS2EVLN(const real *FIT, real *EVLN) const
{
  if ( fQED )
    EVLN[EVLN_GAM] = FIT[FIT_GAM];
  else
    EVLN[EVLN_GAM] = 0;

  EVLN[EVLN_SNG] = FIT[FIT_SNG]; //Singlet
  EVLN[EVLN_GLU] = FIT[FIT_GLU]; //Gluon
  EVLN[EVLN_VAL] = FIT[FIT_VAL]; //Valence
  EVLN[EVLN_V3]  = FIT[FIT_T3] + 2*FIT[FIT_DS]; //V3 = T3+2*Ds
  EVLN[EVLN_V8]  = FIT[FIT_VAL] -3*FIT[FIT_SM]; //V8 = V - 3sm
  EVLN[EVLN_V15] = FIT[FIT_VAL]; // V15
  EVLN[EVLN_V24] = FIT[FIT_VAL]; // V24
  EVLN[EVLN_V35] = FIT[FIT_VAL]; // V35
  EVLN[EVLN_T3]  = FIT[FIT_T3];  //T3
  EVLN[EVLN_T8]  = FIT[FIT_SNG] - 3*FIT[FIT_SP]; //T8 = S - 3sp
  EVLN[EVLN_T15] = FIT[FIT_SNG]; //T15
  EVLN[EVLN_T24] = FIT[FIT_SNG]; //T24
  EVLN[EVLN_T35] = FIT[FIT_SNG]; //T35

  return;
}

void NN23FitBasis::EVLN2BASIS(const real *EVLN, real *FIT) const
{
  // Order in fitting basis
  // S g V T3 Ds sp sm gam

  // Order in Evln bassi
  // γ, Σ, g, V, V3, V8, V15, V24, V35, T3, T8, T15, T24, T35

  //V3 = T3+2*Ds
  //V8 = V - 3sm
  //T8 = S - 3sp

  FIT[FIT_SNG] = EVLN[EVLN_SNG]; //Singlet
  FIT[FIT_GLU] = EVLN[EVLN_GLU]; //gluon
  FIT[FIT_VAL] = EVLN[EVLN_VAL]; //valence
  FIT[FIT_T3]  = EVLN[EVLN_T3];   // T3
  FIT[FIT_DS]  = 0.5*(EVLN[EVLN_V3] - EVLN[EVLN_T3]);   // Ds = (V3-T3)/2

  FIT[FIT_SP]  = (EVLN[EVLN_SNG] - EVLN[EVLN_T8])/3.0;  //sp = (S-T8)/3
  FIT[FIT_SM]  = (EVLN[EVLN_VAL] - EVLN[EVLN_V8])/3.0;  //sm = (V - V8)/3

  if (fQED)
    FIT[FIT_GAM] =  EVLN[EVLN_GAM];  // photon

  return;
}

real NN23FitBasis::ComputeSumRules(sumRule rule, int mem, PDFSet *pdf, bool &status) const
{
  // status
  status = false;

  // sum rule calculations
  switch (rule) {
     // total momentum
     case SUM_MSR:
       {
         real xsng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real xglu = pdf->IntegratePDF(mem,FIT_GLU,fQ2,PDFSet::XFX,status,fGSLWork);
         real msr = xsng+xglu;
         if (fQED)
           {
             real xgam = pdf->IntegratePDF(mem,FIT_GAM,fQ2,PDFSet::XFX,status,fGSLWork);
             msr += xgam;
           }
         return msr;
       }
       break;
     // up valence
     case SUM_UVL:
       {
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real t3 = pdf->IntegratePDF(mem,FIT_T3,fQ2,PDFSet::FX,status,fGSLWork);
         real delta = pdf->IntegratePDF(mem,FIT_DS,fQ2,PDFSet::FX,status,fGSLWork);
         real sm = pdf->IntegratePDF(mem,FIT_SM,fQ2,PDFSet::FX,status,fGSLWork);
         return 0.5*( val + t3 + 2.0*delta - sm );
       }
       break;
     // down valence
     case SUM_DVL:
       {
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real t3 = pdf->IntegratePDF(mem,FIT_T3,fQ2,PDFSet::FX,status,fGSLWork);
         real delta = pdf->IntegratePDF(mem,FIT_DS,fQ2,PDFSet::FX,status,fGSLWork);
         real sm = pdf->IntegratePDF(mem,FIT_SM,fQ2,PDFSet::FX,status,fGSLWork);
         return 0.5*( val - t3 - 2.0*delta - sm );
       }
       break;
     case SUM_SVL:
     // strange valence
        return pdf->IntegratePDF(mem,FIT_SM,fQ2,PDFSet::FX,status,fGSLWork);
        break;
     case SUM_USM:
     // up momentum fraction
       {
         real sng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real t3 = pdf->IntegratePDF(mem,FIT_T3,fQ2,PDFSet::XFX,status,fGSLWork);
         real sp = pdf->IntegratePDF(mem,FIT_SP,fQ2,PDFSet::XFX,status,fGSLWork);
         return 0.5*( sng + t3 - sp );
       }
       break;
     case SUM_DSM:
     // down momentum fraction
       {
         real sng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real t3 = pdf->IntegratePDF(mem,FIT_T3,fQ2,PDFSet::XFX,status,fGSLWork);
         real sp = pdf->IntegratePDF(mem,FIT_SP,fQ2,PDFSet::XFX,status,fGSLWork);
         return 0.5*( sng - t3 - sp );
       }
       break;
     case SUM_SSM:
     // strange momentum fraction
       return pdf->IntegratePDF(mem,FIT_SP,fQ2,PDFSet::XFX,status,fGSLWork);
       break;
     default:
       cerr << "NN23FitBasis::ComputeSumRules error: unknown sum rule"<<endl;
       exit(-1);
       break;
  }
}



// ****************** EVOLUTION BASIS *************************


/**
 *  Evolution Fit Basis
 **/
EvolFitBasis::EvolFitBasis(NNPDFSettings const& nnset):
FitBasis(nnset, "EvolFitBasis", 7+nnset.IsQED()),
fQED(nnset.IsQED())
{
  // PDF Names for plotting
  fPDFNames[FIT_SNG] = "Singlet";
  fPDFNames[FIT_GLU] = "Gluon";
  fPDFNames[FIT_VAL] = "V";
  fPDFNames[FIT_V3] = "V3";
  fPDFNames[FIT_V8] = "V8";
  fPDFNames[FIT_T3] = "T3";
  fPDFNames[FIT_T8] = "T8";
  if (fQED)
    fPDFNames[FIT_GAM] = "Photon";

  // Damping factor for arclengths
  fArcDampFactor[FIT_SNG] = 1;
  fArcDampFactor[FIT_GLU] = 1;
  fArcDampFactor[FIT_VAL] = 0;
  fArcDampFactor[FIT_V3] = 0;
  fArcDampFactor[FIT_V8] = 0;
  fArcDampFactor[FIT_T3] = 1;
  fArcDampFactor[FIT_T8] = 1;
  if (fQED)
    fArcDampFactor[FIT_GAM] = 1;

}

void EvolFitBasis::ComputeParam(PDFSet* pdf, int mem, PreprocParam& param, bool &status) const
{
  // status
  status = false;

  // Clear old normalisations
  for (int i=0; i<fNPDF; i++)
  {
    param.fPDFNorm[i] = 1.0;
    param.fPDFAux[i] = 0.0;
  }

  // Normalisations pointer
  real* norm = param.fPDFNorm;

  // Quantum number sum rules
  norm[FIT_VAL] = 3.0f/pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork); // Total valence
  norm[FIT_V3] = 1.0f/pdf->IntegratePDF(mem,FIT_V3,fQ2,PDFSet::FX,status,fGSLWork); // V3
  norm[FIT_V8] = 3.0f/pdf->IntegratePDF(mem,FIT_V8,fQ2,PDFSet::FX,status,fGSLWork); // V8

  // ************ QED dependent normalisations **************
  if (fQED)
  {
    norm[FIT_GLU] = (1-pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork) - pdf->IntegratePDF(mem,FIT_GAM,fQ2,PDFSet::XFX,status,fGSLWork))
    / pdf->IntegratePDF(mem,FIT_GLU,fQ2,PDFSet::XFX,status,fGSLWork);
  }
  else
  {
    norm[FIT_GLU] = (1-pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork))/
    pdf->IntegratePDF(mem,FIT_GLU,fQ2,PDFSet::XFX,status,fGSLWork);
  }

  return;

}

/**
 * @brief EvolFitBasis::BASIS2EVLN
 * @param FIT
 * @param EVLN
 */
void EvolFitBasis::BASIS2EVLN(const real *FIT, real *EVLN) const
{
  if ( fQED )
    EVLN[EVLN_GAM] = FIT[FIT_GAM];
  else
    EVLN[EVLN_GAM] = 0;

  EVLN[EVLN_SNG]  = FIT[FIT_SNG]; //Singlet
  EVLN[EVLN_GLU]  = FIT[FIT_GLU]; //Gluon
  EVLN[EVLN_VAL]  = FIT[FIT_VAL]; //Valence
  EVLN[EVLN_V3]   = FIT[FIT_V3];  //V3
  EVLN[EVLN_V8]   = FIT[FIT_V8];  //V8
  EVLN[EVLN_V15]  = FIT[FIT_VAL]; //V15 = V
  EVLN[EVLN_V24]  = FIT[FIT_VAL]; //V24 = V
  EVLN[EVLN_V35]  = FIT[FIT_VAL]; //V35 = V
  EVLN[EVLN_T3]   = FIT[FIT_T3];  //T3
  EVLN[EVLN_T8]   = FIT[FIT_T8];  //T8
  EVLN[EVLN_T15]  = FIT[FIT_SNG]; //T15 = S
  EVLN[EVLN_T24]  = FIT[FIT_SNG]; //T24 = S
  EVLN[EVLN_T35]  = FIT[FIT_SNG]; //T35 = S

  return;
}

void EvolFitBasis::EVLN2BASIS(const real *EVLN, real *FIT) const
{
  // Order in fitting basis
  // S g V V3 V8 T3 T8 gam

  // Order in Evln bassi
  // γ, Σ, g, V, V3, V8, V15, V24, V35, T3, T8, T15, T24, T35

  FIT[FIT_SNG]  = EVLN[EVLN_SNG]; //Singlet
  FIT[FIT_GLU]  = EVLN[EVLN_GLU]; //gluon
  FIT[FIT_VAL]  = EVLN[EVLN_VAL]; //valence
  FIT[FIT_V3]   = EVLN[EVLN_V3]; // V3
  FIT[FIT_V8]   = EVLN[EVLN_V8]; // V8
  FIT[FIT_T3]   = EVLN[EVLN_T3]; // T3
  FIT[FIT_T8]   = EVLN[EVLN_T8]; // T8

  if (fQED)
    FIT[FIT_GAM] =  EVLN[0];  // photon

  return;
}

real EvolFitBasis::ComputeSumRules(sumRule rule, int mem, PDFSet *pdf, bool &status) const
{
  // status
  status = false;

  // sum rule calculations
  switch (rule) {
     case SUM_MSR:
       {
         real xsng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real xglu = pdf->IntegratePDF(mem,FIT_GLU,fQ2,PDFSet::XFX,status,fGSLWork);
         real msr = xsng+xglu;
         if (fQED)
           {
             real xgam = pdf->IntegratePDF(mem,FIT_GAM,fQ2,PDFSet::XFX,status,fGSLWork);
             msr += xgam;
           }
         return msr;
       }
       break;
     case SUM_UVL:
       {
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real v3 = pdf->IntegratePDF(mem,FIT_V3,fQ2,PDFSet::FX,status,fGSLWork);
         real v8 = pdf->IntegratePDF(mem,FIT_V8,fQ2,PDFSet::FX,status,fGSLWork);
         return ( 2.0*val + 3.0*v3 + v8 )/6.0;
       }
       break;
     case SUM_DVL:
       {
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real v3 = pdf->IntegratePDF(mem,FIT_V3,fQ2,PDFSet::FX,status,fGSLWork);
         real v8 = pdf->IntegratePDF(mem,FIT_V8,fQ2,PDFSet::FX,status,fGSLWork);
         return ( 2.0*val - 3.0*v3 + v8 )/6.0;
       }
       break;
     case SUM_SVL:
       {
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real v8 = pdf->IntegratePDF(mem,FIT_V8,fQ2,PDFSet::FX,status,fGSLWork);
         return ( val - v8 )/3.0;
       }
        break;
     case SUM_USM:
       {
         real sng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real t3 = pdf->IntegratePDF(mem,FIT_T3,fQ2,PDFSet::XFX,status,fGSLWork);
         real t8 = pdf->IntegratePDF(mem,FIT_T8,fQ2,PDFSet::XFX,status,fGSLWork);
         return ( 2.0*sng + 3.0*t3 + t8 )/6.0;
       }
       break;
     case SUM_DSM:
       {
         real sng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real t3 = pdf->IntegratePDF(mem,FIT_T3,fQ2,PDFSet::XFX,status,fGSLWork);
         real t8 = pdf->IntegratePDF(mem,FIT_T8,fQ2,PDFSet::XFX,status,fGSLWork);
         return ( 2.0*sng - 3.0*t3 + t8 )/6.0;
       }
       break;
     case SUM_SSM:
       {
         real sng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real t8 = pdf->IntegratePDF(mem,FIT_T8,fQ2,PDFSet::XFX,status,fGSLWork);
         return ( sng - t8 )/3.0;
       }
       break;
     default:
       cerr << "EvolFitBasis::ComputeSumRules error: unknown sum rule"<<endl;
       exit(-1);
       break;
  }
}

// ******************** DIS EVOLUTION BASIS (ISOSCALAR CASE)*******


/**
 *  DIS Evolution Fit Basis
 **/
DISEvolFitBasis::DISEvolFitBasis(NNPDFSettings const& nnset):
FitBasis(nnset, "EvolFitBasis", 3)
{
  //PDF names for plotting
  fPDFNames[FIT_SNG] = "Singlet";
  fPDFNames[FIT_GLU] = "Gluon";
  //fPDFNames[FIT_T3] = "T3";
  fPDFNames[FIT_T8] = "T8";

  // Damping factor for arclengths
  fArcDampFactor[FIT_SNG] = 1;
  fArcDampFactor[FIT_GLU] = 1;
  //fArcDampFactor[FIT_T3] = 1;
  fArcDampFactor[FIT_T8] = 1;
}

void DISEvolFitBasis::ComputeParam(PDFSet* pdf, int mem, PreprocParam& param, bool &status) const
{
  // status
  status = false;

  // Clear old normalisations
  for (int i=0; i<fNPDF; i++)
  {
    param.fPDFNorm[i] = 1.0;
    param.fPDFAux[i] = 0.0;
  }

  // Normalisations pointer
  real* norm = param.fPDFNorm;

  // ************ QED dependent normalisations **************
    norm[FIT_GLU] = (1-pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork))/
    pdf->IntegratePDF(mem,FIT_GLU,fQ2,PDFSet::XFX,status,fGSLWork);

  return;
}

/**
 * @brief DISEvolFitBasis::BASIS2EVLN
 * @param FIT
 * @param EVLN
 */
void DISEvolFitBasis::BASIS2EVLN(const real *FIT, real *EVLN) const
{

  EVLN[EVLN_GAM] = 0;
  EVLN[EVLN_SNG]  = FIT[FIT_SNG]; //Singlet
  EVLN[EVLN_GLU]  = FIT[FIT_GLU]; //Gluon
  EVLN[EVLN_VAL]  = 0; //Valence
  EVLN[EVLN_V3]   = 0; //V3 = V
  EVLN[EVLN_V8]   = 0; //V8 = V
  EVLN[EVLN_V15]  = 0; //V15 = V
  EVLN[EVLN_V24]  = 0; //V24 = V
  EVLN[EVLN_V35]  = 0; //V35 = V
  EVLN[EVLN_T3]   = 0;  //T3
  EVLN[EVLN_T8]   = FIT[FIT_T8];  //T8
  EVLN[EVLN_T15]  = 0; //T15 = S
  EVLN[EVLN_T24]  = 0; //T24 = S
  EVLN[EVLN_T35]  = 0; //T35 = S

  return;
}

void DISEvolFitBasis::EVLN2BASIS(const real *EVLN, real *FIT) const
{
  // Order in fitting basis
  // S g V V3 V8 T3 T8 gam

  // Order in Evln bassi
  // γ, Σ, g, V, V3, V8, V15, V24, V35, T3, T8, T15, T24, T35

  FIT[FIT_SNG]  = EVLN[EVLN_SNG]; //Singlet
  FIT[FIT_GLU]  = EVLN[EVLN_GLU]; //gluon
  FIT[FIT_T8]   = EVLN[EVLN_T8]; // T8

  return;
}

real DISEvolFitBasis::ComputeSumRules(sumRule rule, int mem, PDFSet *pdf, bool &status) const
{
// status
  status = false;

  // sum rule calculations
  switch (rule) {
     case SUM_MSR:
       {
         real xsng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real xglu = pdf->IntegratePDF(mem,FIT_GLU,fQ2,PDFSet::XFX,status,fGSLWork);
         real msr = xsng+xglu;
         return msr;
       }
       break;
     case SUM_UVL:
       {
        /*
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real v3 = pdf->IntegratePDF(mem,FIT_V3,fQ2,PDFSet::FX,status,fGSLWork);
         real v8 = pdf->IntegratePDF(mem,FIT_V8,fQ2,PDFSet::FX,status,fGSLWork);
         return ( 2.0*val + 3.0*v3 + v8 )/6.0;
         */
        return 0;
       }
       break;
     case SUM_DVL:
       {
        /*
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real v3 = pdf->IntegratePDF(mem,FIT_V3,fQ2,PDFSet::FX,status,fGSLWork);
         real v8 = pdf->IntegratePDF(mem,FIT_V8,fQ2,PDFSet::FX,status,fGSLWork);
         return ( 2.0*val - 3.0*v3 + v8 )/6.0;
         */
        return 0;
       }
       break;
     case SUM_SVL:
       {
        /*
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real v8 = pdf->IntegratePDF(mem,FIT_V8,fQ2,PDFSet::FX,status,fGSLWork);
         return ( val - v8 )/3.0;
         */
        return 0;
       }
        break;
     case SUM_USM:
       {
         real sng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real t3 = 0; //pdf->IntegratePDF(mem,FIT_T3,fQ2,PDFSet::XFX,status,fGSLWork);
         real t8 = pdf->IntegratePDF(mem,FIT_T8,fQ2,PDFSet::XFX,status,fGSLWork);
         return ( 2.0*sng + 3.0*t3 + t8 )/6.0;
       }
       break;
     case SUM_DSM:
       {
         real sng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real t3 = 0; //pdf->IntegratePDF(mem,FIT_T3,fQ2,PDFSet::XFX,status,fGSLWork);
         real t8 = pdf->IntegratePDF(mem,FIT_T8,fQ2,PDFSet::XFX,status,fGSLWork);
         return ( 2.0*sng - 3.0*t3 + t8 )/6.0;
       }
       break;
     case SUM_SSM:
       {
         real sng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real t8 = pdf->IntegratePDF(mem,FIT_T8,fQ2,PDFSet::XFX,status,fGSLWork);
         return ( sng - t8 )/3.0;
       }
       break;
     default:
       cerr << "EvolFitBasis::ComputeSumRules error: unknown sum rule"<<endl;
       exit(-1);
       break;
  }
}
// ******************** EVOLUTION with LUXQED *********************

/**
 * @brief LuxBasis::LuxBasis
 */
LuxBasis::LuxBasis(NNPDFSettings const& nnset, int const& replica):
  FitBasis(nnset, "LuxBasis", 8),
  fQ0(stod(nnset.GetTheory(APFEL::kQ0)))
{
  // photon from t0
  fPhotonSet = new LHAPDFSet(nnset.Get("datacuts","t0pdfset").as<string>(), replica);

  // PDF Names for plotting
  fPDFNames[FIT_SNG] = "Singlet";
  fPDFNames[FIT_GLU] = "Gluon";
  fPDFNames[FIT_VAL] = "V";
  fPDFNames[FIT_V3] = "V3";
  fPDFNames[FIT_V8] = "V8";
  fPDFNames[FIT_T3] = "T3";
  fPDFNames[FIT_T8] = "T8";
  fPDFNames[FIT_CP] = "c+";

  // Damping factor for arclengths
  fArcDampFactor[FIT_SNG] = 1;
  fArcDampFactor[FIT_GLU] = 1;
  fArcDampFactor[FIT_VAL] = 0;
  fArcDampFactor[FIT_V3] = 0;
  fArcDampFactor[FIT_V8] = 0;
  fArcDampFactor[FIT_T3] = 1;
  fArcDampFactor[FIT_T8] = 1;
  fArcDampFactor[FIT_CP] = 1;
}

/**
 * @brief LuxBasis::~LuxBasis
 */
LuxBasis::~LuxBasis()
{
  delete fPhotonSet;
}

/**
 * @brief LuxBasis::BASIS2EVLN
 */
void LuxBasis::BASIS2EVLN(const real *FIT, real *EVLN) const
{
  EVLN[EVLN_GAM]  = FIT[FIT_GAM]; //Photon
  EVLN[EVLN_SNG]  = FIT[FIT_SNG]; //Singlet
  EVLN[EVLN_GLU]  = FIT[FIT_GLU]; //Gluon
  EVLN[EVLN_VAL]  = FIT[FIT_VAL]; //Valence
  EVLN[EVLN_V3]   = FIT[FIT_V3];  //V3
  EVLN[EVLN_V8]   = FIT[FIT_V8];  //V8
  EVLN[EVLN_V15]  = FIT[FIT_VAL]; //V15 = V (c- = 0)
  EVLN[EVLN_V24]  = FIT[FIT_VAL]; //V24 = V
  EVLN[EVLN_V35]  = FIT[FIT_VAL]; //V35 = V
  EVLN[EVLN_T3]   = FIT[FIT_T3];  //T3
  EVLN[EVLN_T8]   = FIT[FIT_T8];  //T8
  EVLN[EVLN_T15]  = FIT[FIT_SNG] - 4*FIT[FIT_CP]; //T15
  EVLN[EVLN_T24]  = FIT[FIT_SNG]; //T24 = S
  EVLN[EVLN_T35]  = FIT[FIT_SNG]; //T35 = S
}

/**
 * @brief LuxBasis::EVLN2BASIS
 */
void LuxBasis::EVLN2BASIS(const real *EVLN, real *FIT) const
{
  // Order in fitting basis
  // S g V V3 V8 T3 T8 c+ gam

  // Order in Evln bassi
  // γ, Σ, g, V, V3, V8, V15, V24, V35, T3, T8, T15, T24, T35

  FIT[FIT_SNG]  = EVLN[EVLN_SNG]; //Singlet
  FIT[FIT_GLU]  = EVLN[EVLN_GLU]; //gluon
  FIT[FIT_VAL]  = EVLN[EVLN_VAL]; //valence
  FIT[FIT_V3]   = EVLN[EVLN_V3];  //V3
  FIT[FIT_V8]   = EVLN[EVLN_V8];  //V8
  FIT[FIT_T3]   = EVLN[EVLN_T3];  //T3
  FIT[FIT_T8]   = EVLN[EVLN_T8];  //T8
  FIT[FIT_CP]   = (EVLN[EVLN_SNG]-EVLN[EVLN_T15])/4;  // T15
  FIT[FIT_GAM]  = EVLN[EVLN_GAM]; //photon
}

/**
 * @brief LuxBasis::ComputeParam
 */
void LuxBasis::ComputeParam(PDFSet* pdf, int mem, PreprocParam& param, bool &status) const
{
  // status
  status = false;

  // Clear old normalisations
  for (int i=0; i<fNPDF; i++)
  {
    param.fPDFNorm[i] = 1.0;
    param.fPDFAux[i] = 0.0;
  }

  // Normalisations pointer
  real* norm = param.fPDFNorm;

  // Quantum number sum rules
  norm[FIT_VAL] = 3.0f/pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork); // Total valence
  norm[FIT_V3] = 1.0f/pdf->IntegratePDF(mem,FIT_V3,fQ2,PDFSet::FX,status,fGSLWork); // V3
  norm[FIT_V8] = 3.0f/pdf->IntegratePDF(mem,FIT_V8,fQ2,PDFSet::FX,status,fGSLWork); // V8
  norm[FIT_GLU] = (1-pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork) - pdf->IntegratePDF(mem,FIT_GAM,fQ2,PDFSet::XFX,status,fGSLWork))
  / pdf->IntegratePDF(mem,FIT_GLU,fQ2,PDFSet::XFX,status,fGSLWork);

}

/**
 * @brief LuxBasis::ComputeSumRules
 */
real LuxBasis::ComputeSumRules(sumRule rule, int mem, PDFSet *pdf, bool &status) const
{
  // status
  status = false;

  // sum rule calculations
  switch (rule) {
     case SUM_MSR:
       {
         real xsng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real xglu = pdf->IntegratePDF(mem,FIT_GLU,fQ2,PDFSet::XFX,status,fGSLWork);
         real xgam = pdf->IntegratePDF(mem,FIT_GAM,fQ2,PDFSet::XFX,status,fGSLWork);
         real msr = xsng+xglu+xgam;
         return msr;
       }
       break;
     case SUM_UVL:
       {
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real v3 = pdf->IntegratePDF(mem,FIT_V3,fQ2,PDFSet::FX,status,fGSLWork);
         real v8 = pdf->IntegratePDF(mem,FIT_V8,fQ2,PDFSet::FX,status,fGSLWork);
         real v15 = val;
         return ( 3.0*val + 6.0*v3 + 2.0*v8 + v15 )/12.0;
       }
       break;
     case SUM_DVL:
       {
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real v3 = pdf->IntegratePDF(mem,FIT_V3,fQ2,PDFSet::FX,status,fGSLWork);
         real v8 = pdf->IntegratePDF(mem,FIT_V8,fQ2,PDFSet::FX,status,fGSLWork);
         real v15 = val;
         return ( 3.0*val - 6.0*v3 + 2.0*v8 + v15 )/12.0;
       }
       break;
     case SUM_SVL:
       {
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real v8 = pdf->IntegratePDF(mem,FIT_V8,fQ2,PDFSet::FX,status,fGSLWork);
         real v15 = val;
         return ( 3.0*val - 4.0*v8 + v15 )/12.0;
       }
       break;
      case SUM_CVL:
       {
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real v15 = val;
         return ( val - v15 )/4.0;
       }
       break;
     case SUM_USM:
       {
         real sng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real t3 = pdf->IntegratePDF(mem,FIT_T3,fQ2,PDFSet::XFX,status,fGSLWork);
         real t8 = pdf->IntegratePDF(mem,FIT_T8,fQ2,PDFSet::XFX,status,fGSLWork);
         real cp = pdf->IntegratePDF(mem,FIT_CP,fQ2,PDFSet::XFX,status,fGSLWork);
         real t15 = sng-4*cp;
         return ( 3.0*sng + 6.0*t3 + 2.0*t8 + t15 )/12.0;
      }
       break;
     case SUM_DSM:
       {
         real sng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real t3 = pdf->IntegratePDF(mem,FIT_T3,fQ2,PDFSet::XFX,status,fGSLWork);
         real t8 = pdf->IntegratePDF(mem,FIT_T8,fQ2,PDFSet::XFX,status,fGSLWork);
         real cp = pdf->IntegratePDF(mem,FIT_CP,fQ2,PDFSet::XFX,status,fGSLWork);
         real t15 = sng-4*cp;
         return ( 3.0*sng - 6.0*t3 + 2.0*t8 + t15)/12.0;
       }
       break;
     case SUM_SSM:
       {
         real sng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real t8 = pdf->IntegratePDF(mem,FIT_T8,fQ2,PDFSet::XFX,status,fGSLWork);
         real cp = pdf->IntegratePDF(mem,FIT_CP,fQ2,PDFSet::XFX,status,fGSLWork);
         real t15 = sng-4*cp;
         return ( 3.0*sng - 4.0*t8 + t15)/12.0;
       }
       break;
     case SUM_CSM:
       {
        real cp = pdf->IntegratePDF(mem,FIT_CP,fQ2,PDFSet::XFX,status,fGSLWork);
        return cp;
       }
      break;
     default:
       cerr << "LuxBasis::ComputeSumRules error: unknown sum rule"<<endl;
       exit(-1);
       break;
  }
}

/**
 * @brief LuxBasis::Preprocess
 */
void LuxBasis::Preprocess(const real &x, const int &fl, real &pdf, const PreprocParam &par)
{
  FitBasis::Preprocess(x,fl,pdf,par);
  if (fl == FIT_GAM)
    pdf = fPhotonSet->xfxQ(x, fQ0, 0, 22);
}

/**
 * @brief FitBasis::Preprocess
 */
void LuxBasis::Preprocess(real const& x, real* pdf, PreprocParam const& par)
{
  FitBasis::Preprocess(x, pdf, par);

  // Force photon preprocess
  Preprocess(x,fNPDF,pdf[fNPDF],par);
}


// ******************** EVOLUTION ONLY STRANGENESS *********************

/**
 *  Evolution only for strangeness
 **/

EvolSFitBasis::EvolSFitBasis(NNPDFSettings const& nnset):
FitBasis(nnset, "EvolSFitBasis", 7+nnset.IsQED()),
fQED(nnset.IsQED())
{
  // PDF Names for plotting
  fPDFNames[FIT_SNG] = "Singlet";
  fPDFNames[FIT_GLU] = "Gluon";
  fPDFNames[FIT_VAL] = "Valence";
  fPDFNames[FIT_V8] = "V8";
  fPDFNames[FIT_T3] = "Triplet";
  fPDFNames[FIT_T8] = "T8";
  fPDFNames[FIT_DS] = "Sea Asymmetry";
  if (fQED)
    fPDFNames[FIT_GAM] = "Photon";

  // Damping factor for arclengths
  fArcDampFactor[FIT_SNG] = 1;
  fArcDampFactor[FIT_GLU] = 1;
  fArcDampFactor[FIT_VAL] = 0;
  fArcDampFactor[FIT_V8] = 0;
  fArcDampFactor[FIT_T3] = 1;
  fArcDampFactor[FIT_T8] = 1;
  fArcDampFactor[FIT_DS] = 1;
  if (fQED)
    fArcDampFactor[FIT_GAM] = 1;
}

void EvolSFitBasis::ComputeParam(PDFSet* pdf, int mem, PreprocParam& param, bool &status) const
{
  // status
  status = false;

  // Clear old normalisations
  for (int i=0; i<fNPDF; i++)
  {
    param.fPDFNorm[i] = 1.0;
    param.fPDFAux[i] = 0.0;
  }

  // Normalisations pointer
  real* norm = param.fPDFNorm;

  // Quantum number sum rules
  norm[FIT_VAL] = 3.0f/pdf->IntegratePDF(mem,FIT_VAL, fQ2,PDFSet::FX,status,fGSLWork); // Total valence
  norm[FIT_V8] = 3.0f/pdf->IntegratePDF(mem,FIT_V8, fQ2,PDFSet::FX,status,fGSLWork); // V8
  norm[ FIT_DS] = (1-pdf->IntegratePDF(mem,FIT_T3, fQ2,PDFSet::FX,status,fGSLWork))/
  (2*pdf->IntegratePDF(mem,FIT_DS, fQ2,PDFSet::FX,status,fGSLWork)); // D_S - 1-t3/2d_s

  // ************ QED dependent normalisations **************
  if (fQED)
  {
    norm[FIT_GLU] = (1-pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork) - pdf->IntegratePDF(mem,FIT_GAM,fQ2,PDFSet::XFX,status,fGSLWork))
    / pdf->IntegratePDF(mem,FIT_GLU, fQ2,PDFSet::XFX,status,fGSLWork);
  }
  else
  {
    norm[FIT_GLU] = (1-pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork))/
    pdf->IntegratePDF(mem,FIT_GLU,fQ2,PDFSet::XFX,status,fGSLWork);
  }

  return;

}

/**
 * @brief NN23FitBasis::BASIS2EVLN
 * @param FIT
 * @param EVLN
 */
void EvolSFitBasis::BASIS2EVLN(const real *FIT, real *EVLN) const
{
  if ( fQED )
    EVLN[EVLN_GAM] = FIT[FIT_GAM];
  else
    EVLN[EVLN_GAM] = 0;


  enum fitBasis {FIT_SNG, FIT_GLU, FIT_VAL, FIT_V8, FIT_T3, FIT_T8, FIT_DS, FIT_GAM };

  EVLN[EVLN_SNG] = FIT[FIT_SNG]; //Singlet
  EVLN[EVLN_GLU] = FIT[FIT_GLU]; //Gluon
  EVLN[EVLN_VAL] = FIT[FIT_VAL]; //Valence
  EVLN[EVLN_V3]  = FIT[FIT_T3]+2*FIT[FIT_DS]; //V3 = T3+2*Ds
  EVLN[EVLN_V8]  = FIT[FIT_V8];  //V8
  EVLN[EVLN_V15] = FIT[FIT_VAL]; //V15 = V
  EVLN[EVLN_V24] = FIT[FIT_VAL]; //V24 = V
  EVLN[EVLN_V35] = FIT[FIT_VAL]; //V35 = V
  EVLN[EVLN_T3]  = FIT[FIT_T3];  //T3
  EVLN[EVLN_T8]  = FIT[FIT_T8];  //T8
  EVLN[EVLN_T15] = FIT[FIT_SNG]; //T15 = S
  EVLN[EVLN_T24] = FIT[FIT_SNG]; //T24 = S
  EVLN[EVLN_T35] = FIT[FIT_SNG]; //T35 = S

  return;
}

void EvolSFitBasis::EVLN2BASIS(const real *EVLN, real *FIT) const
{
  // Order in fitting basis
  // S g V T3 Ds sp sm gam

  // Order in Evln bassi
  // γ, Σ, g, V, V3, V8, V15, V24, V35, T3, T8, T15, T24, T35

  //V3 = T3+2*Ds
  //V8 = V - 3sm
  //T8 = S - 3sp

  FIT[FIT_SNG]  = EVLN[EVLN_SNG]; //Singlet
  FIT[FIT_GLU]  = EVLN[EVLN_GLU]; //gluon
  FIT[FIT_VAL]  = EVLN[EVLN_VAL]; //valence
  FIT[FIT_V8]  = EVLN[EVLN_V8]; // v8
  FIT[FIT_T3]  = EVLN[EVLN_T3]; // t3
  FIT[FIT_T8]  = EVLN[EVLN_T8]; // t8
  FIT[FIT_DS]  = 0.5*(EVLN[EVLN_V3] - EVLN[EVLN_T3]); // Ds = (V3-T3)/2

  if (fQED)
    FIT[FIT_GAM] =  EVLN[EVLN_GAM];  // photon

  return;
}

real EvolSFitBasis::ComputeSumRules(sumRule rule, int mem, PDFSet *pdf, bool &status) const
{
  // status
  status = false;

  // sum rule calculations
  switch (rule) {
     case SUM_MSR:
       {
         real xsng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real xglu = pdf->IntegratePDF(mem,FIT_GLU,fQ2,PDFSet::XFX,status,fGSLWork);
         real msr = xsng+xglu;
         if (fQED)
           {
             real xgam = pdf->IntegratePDF(mem,FIT_GAM,fQ2,PDFSet::XFX,status,fGSLWork);
             msr += xgam;
           }
         return msr;
       }
       break;
     case SUM_UVL:
       {
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real t3 = pdf->IntegratePDF(mem,FIT_T3,fQ2,PDFSet::FX,status,fGSLWork);
         real ds = pdf->IntegratePDF(mem,FIT_DS,fQ2,PDFSet::FX,status,fGSLWork);
         real v8 = pdf->IntegratePDF(mem,FIT_V8,fQ2,PDFSet::FX,status,fGSLWork);
         return ( 2.0*val + 3.0*t3 + 6.0*ds + v8 )/6.0;
       }
       break;
     case SUM_DVL:
       {
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real t3 = pdf->IntegratePDF(mem,FIT_T3,fQ2,PDFSet::FX,status,fGSLWork);
         real ds = pdf->IntegratePDF(mem,FIT_DS,fQ2,PDFSet::FX,status,fGSLWork);
         real v8 = pdf->IntegratePDF(mem,FIT_V8,fQ2,PDFSet::FX,status,fGSLWork);
         return ( 2.0*val - 3.0*t3 - 6.0*ds + v8 )/6.0;
       }
       break;
     case SUM_SVL:
       {
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real v8 = pdf->IntegratePDF(mem,FIT_V8,fQ2,PDFSet::FX,status,fGSLWork);
         return ( val - v8 )/3.0;
       }
        break;
     case SUM_USM:
       {
         real sng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real t3 = pdf->IntegratePDF(mem,FIT_T3,fQ2,PDFSet::XFX,status,fGSLWork);
         real t8 = pdf->IntegratePDF(mem,FIT_T8,fQ2,PDFSet::XFX,status,fGSLWork);
         return ( 2.0*sng + 3.0*t3 + t8 )/6.0;
       }
       break;
     case SUM_DSM:
       {
         real sng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real t3 = pdf->IntegratePDF(mem,FIT_T3,fQ2,PDFSet::XFX,status,fGSLWork);
         real t8 = pdf->IntegratePDF(mem,FIT_T8,fQ2,PDFSet::XFX,status,fGSLWork);
         return ( 2.0*sng - 3.0*t3 + t8 )/6.0;
       }
       break;
     case SUM_SSM:
       {
         real sng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real t8 = pdf->IntegratePDF(mem,FIT_T8,fQ2,PDFSet::XFX,status,fGSLWork);
         return ( sng - t8 )/3.0;
       }
       break;
     default:
       cerr << "EvolSFitBasis::ComputeSumRules error: unknown sum rule"<<endl;
       exit(-1);
       break;
  }
}

NN30FitBasis::NN30FitBasis(NNPDFSettings const& nnset):
EvolFitBasis(nnset)
{
  // PDF Names for plotting
  fPDFNames[NET_SNG] = "Singlet";
  fPDFNames[NET_GLU] = "Gluon";
  fPDFNames[NET_VAL] = "Valence";
  fPDFNames[NET_T3] = "Triplet";
  fPDFNames[NET_DS] = "Sea Asymmetry";
  fPDFNames[NET_SP] = "Strange Sea";
  fPDFNames[NET_SM] = "Strange Valence";
  if (fQED)
    fPDFNames[NET_GAM] = "Photon";
}

void NN30FitBasis::Preprocess(real const& x, real* pdf, PreprocParam const& par)
{
   // Positive definite PDFs
   for (int i = 0; i < fNPDF; i++)
     if (fPDFSqrPos[i])
       pdf[i] *= pdf[i];

   // Transform to evol basis
   real tmppdf[15];
   tmppdf[FIT_SNG] = pdf[NET_SNG];
   tmppdf[FIT_GLU] = pdf[NET_GLU];
   tmppdf[FIT_VAL] = pdf[NET_VAL];
   tmppdf[FIT_V3] = pdf[NET_T3]+2*pdf[NET_DS];
   tmppdf[FIT_V8] = pdf[NET_VAL]-3*pdf[NET_SM];
   tmppdf[FIT_T3] = pdf[NET_T3];
   tmppdf[FIT_T8] = pdf[NET_SNG]-3*pdf[NET_SP];
   if (fQED)
     tmppdf[FIT_GAM] = pdf[NET_GAM];

   for (int i = 0; i < fNPDF; i++)
     pdf[i] = tmppdf[i];

   // Preproccess in evol basis
   for (int i = 0; i < fNPDF; i++)
     FitBasis::Preprocess(x,i,pdf[i],par);

   return;
}

void NN30FitBasis::NetTransform(int const& fl, int const& nfl, int* transform)
{
  for (int i = 0; i < nfl; i++)
    transform[i] = 0;

  switch (fl)
  {
    case FIT_SNG:
    {
      transform[NET_SNG] = 1;
      break;
    }
    case FIT_GLU:
    {
      transform[NET_GLU] = 1;
      break;
    }
    case FIT_VAL:
    {
      transform[NET_VAL] = 1;
      break;
    }
    case FIT_V3:
    {
      transform[NET_T3] = 1;
      transform[NET_DS] = 2;
      break;
    }
    case FIT_V8:
    {
      transform[NET_VAL] = 1;
      transform[NET_SM] = -3;
      break;
    }
    case FIT_T3:
    {
      transform[NET_T3] = 1;
      break;
    }
    case FIT_T8:
    {
      transform[NET_SNG] = 1;
      transform[NET_SP] = -3;
      break;
    }
    case FIT_GAM:
    {
      transform[NET_GAM] = 1;
      break;
    }
  }

  return;
}

// Flavour basis
FLVRFitBasis::FLVRFitBasis(NNPDFSettings const& nnset):
EvolFitBasis(nnset)
{
  // PDF Names for plotting
  fPDFNames[NET_GLU] = "Gluon";
  fPDFNames[NET_U] = "Up";
  fPDFNames[NET_UBAR] = "AntiUp";
  fPDFNames[NET_D] = "Down";
  fPDFNames[NET_DBAR] = "AntiDown";
  fPDFNames[NET_S] = "Strange";
  fPDFNames[NET_SBAR] = "AntiStrange";
  if (fQED)
    fPDFNames[NET_GAM] = "Photon";
}

void FLVRFitBasis::Preprocess(real const& x, real* pdf, PreprocParam const& par)
{
   // Positive definite PDFs
   for (int i = 0; i < fNPDF; i++)
     if (fPDFSqrPos[i])
       pdf[i] *= pdf[i];

   // Transform to evol basis
   real tmppdf[15];
   tmppdf[FIT_SNG] = pdf[NET_U]+pdf[NET_UBAR]+pdf[NET_D]+pdf[NET_DBAR]+pdf[NET_S]+pdf[NET_SBAR];
   tmppdf[FIT_GLU] = pdf[NET_GLU];
   tmppdf[FIT_VAL] = pdf[NET_U]-pdf[NET_UBAR]+pdf[NET_D]-pdf[NET_DBAR]+pdf[NET_S]-pdf[NET_SBAR];
   tmppdf[FIT_V3] = pdf[NET_U]-pdf[NET_UBAR]-pdf[NET_D]+pdf[NET_DBAR];
   tmppdf[FIT_V8] = pdf[NET_U]-pdf[NET_UBAR]+pdf[NET_D]-pdf[NET_DBAR]-2*pdf[NET_S]+2*pdf[NET_SBAR];
   tmppdf[FIT_T3] = pdf[NET_U]+pdf[NET_UBAR]-pdf[NET_D]-pdf[NET_DBAR];
   tmppdf[FIT_T8] = pdf[NET_U]+pdf[NET_UBAR]+pdf[NET_D]+pdf[NET_DBAR]-2*pdf[NET_S]-2*pdf[NET_SBAR];
   if (fQED)
     tmppdf[FIT_GAM] = pdf[NET_GAM];

   for (int i = 0; i < fNPDF; i++)
     pdf[i] = tmppdf[i];

   // Preprocess in evol basis
   for (int i = 0; i < fNPDF; i++)
     FitBasis::Preprocess(x,i,pdf[i],par);

   return;
}

void FLVRFitBasis::NetTransform(int const& fl, int const& nfl, int* transform)
{
  for (int i = 0; i < nfl; i++)
    transform[i] = 0;

  switch (fl)
  {
    case FIT_SNG:
    {
      transform[NET_U] = 1;
      transform[NET_UBAR] = 1;
      transform[NET_D] = 1;
      transform[NET_DBAR] = 1;
      transform[NET_S] = 1;
      transform[NET_SBAR] = 1;
      break;
    }
    case FIT_GLU:
    {
      transform[NET_GLU] = 1;
      break;
    }
    case FIT_VAL:
    {
      transform[NET_U] = 1;
      transform[NET_UBAR] = -1;
      transform[NET_D] = 1;
      transform[NET_DBAR] = -1;
      transform[NET_S] = 1;
      transform[NET_SBAR] = -1;
      break;
    }
    case FIT_V3:
    {
      transform[NET_U] = 1;
      transform[NET_UBAR] = -1;
      transform[NET_D] = -1;
      transform[NET_DBAR] = 1;
      break;
    }
    case FIT_V8:
    {
      transform[NET_U] = 1;
      transform[NET_UBAR] = -1;
      transform[NET_D] = 1;
      transform[NET_DBAR] = -1;
      transform[NET_S] = -2;
      transform[NET_SBAR] = 2;
      break;
    }
    case FIT_T3:
    {
      transform[NET_U] = 1;
      transform[NET_UBAR] = 1;
      transform[NET_D] = -1;
      transform[NET_DBAR] = -1;
      break;
    }
    case FIT_T8:
    {
      transform[NET_U] = 1;
      transform[NET_UBAR] = 1;
      transform[NET_D] = 1;
      transform[NET_DBAR] = 1;
      transform[NET_S] = -2;
      transform[NET_SBAR] = -2;
      break;
    }
    case FIT_GAM:
    {
      transform[NET_GAM] = 1;
      break;
    }
  }

  return;
}

/************************ Intrinsic Charm ****************************/
/**
 *  Evolution Fit Basis
 **/
EvolICFitBasis::EvolICFitBasis(NNPDFSettings const& nnset):
FitBasis(nnset, "EvolICFitBasis", 8+nnset.IsQED()),
fQED(nnset.IsQED())
{
  // PDF Names for plotting
  fPDFNames[FIT_SNG] = "Singlet";
  fPDFNames[FIT_GLU] = "Gluon";
  fPDFNames[FIT_VAL] = "V";
  fPDFNames[FIT_V3] = "V3";
  fPDFNames[FIT_V8] = "V8";
  fPDFNames[FIT_T3] = "T3";
  fPDFNames[FIT_T8] = "T8";
  fPDFNames[FIT_T15]= "T15";
  if (fQED)
    fPDFNames[FIT_GAM] = "Photon";

  // Damping factor for arclengths
  fArcDampFactor[FIT_SNG] = 1;
  fArcDampFactor[FIT_GLU] = 1;
  fArcDampFactor[FIT_VAL] = 0;
  fArcDampFactor[FIT_V3] = 0;
  fArcDampFactor[FIT_V8] = 0;
  fArcDampFactor[FIT_T3] = 1;
  fArcDampFactor[FIT_T8] = 1;
  fArcDampFactor[FIT_T15] = 1;
  if (fQED)
    fArcDampFactor[FIT_GAM] = 1;
}

void EvolICFitBasis::ComputeParam(PDFSet* pdf, int mem, PreprocParam& param, bool &status) const
{
  // status
  status = false;

  // Clear old normalisations
  for (int i=0; i<fNPDF; i++)
  {
    param.fPDFNorm[i] = 1.0;
    param.fPDFAux[i] = 0.0;
  }

  // Normalisations pointer
  real* norm = param.fPDFNorm;

  // Quantum number sum rules
  norm[FIT_VAL] = 3.0f/pdf->IntegratePDF(mem,FIT_VAL, fQ2,PDFSet::FX,status,fGSLWork); // Total valence
  norm[FIT_V3] = 1.0f/pdf->IntegratePDF(mem,FIT_V3, fQ2,PDFSet::FX,status,fGSLWork); // V3
  norm[FIT_V8] = 3.0f/pdf->IntegratePDF(mem,FIT_V8, fQ2,PDFSet::FX,status,fGSLWork); // V8

  // ************ QED dependent normalisations **************
  if (fQED)
  {
    norm[FIT_GLU] = (1-pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork) - pdf->IntegratePDF(mem,FIT_GAM,fQ2,PDFSet::XFX,status,fGSLWork))
    / pdf->IntegratePDF(mem,FIT_GLU, fQ2,PDFSet::XFX,status,fGSLWork);
  }
  else
  {
    norm[FIT_GLU] = (1-pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork))/
    pdf->IntegratePDF(mem,FIT_GLU,fQ2,PDFSet::XFX,status,fGSLWork);
  }

  return;

}

/**
 * @brief EvolFitBasis::BASIS2EVLN
 * @param FIT
 * @param EVLN
 */
void EvolICFitBasis::BASIS2EVLN(const real *FIT, real *EVLN) const
{
  if ( fQED )
    EVLN[EVLN_GAM] = FIT[FIT_GAM];
  else
    EVLN[EVLN_GAM] = 0;

  EVLN[EVLN_SNG]  = FIT[FIT_SNG]; //Singlet
  EVLN[EVLN_GLU]  = FIT[FIT_GLU]; //Gluon
  EVLN[EVLN_VAL]  = FIT[FIT_VAL]; //Valence
  EVLN[EVLN_V3]   = FIT[FIT_V3];  //V3
  EVLN[EVLN_V8]   = FIT[FIT_V8];  //V8
  EVLN[EVLN_V15]  = FIT[FIT_VAL]; //V15 = V (c- = 0)
  EVLN[EVLN_V24]  = FIT[FIT_VAL]; //V24 = V
  EVLN[EVLN_V35]  = FIT[FIT_VAL]; //V35 = V
  EVLN[EVLN_T3]   = FIT[FIT_T3];  //T3
  EVLN[EVLN_T8]   = FIT[FIT_T8];  //T8
  EVLN[EVLN_T15]  = FIT[FIT_T15]; //T15
  EVLN[EVLN_T24]  = FIT[FIT_SNG]; //T24 = S
  EVLN[EVLN_T35]  = FIT[FIT_SNG]; //T35 = S

  return;
}

void EvolICFitBasis::EVLN2BASIS(const real *EVLN, real *FIT) const
{
  // Order in fitting basis
  // S g V V3 V8 T3 T8 gam

  // Order in Evln bassi
  // γ, Σ, g, V, V3, V8, V15, V24, V35, T3, T8, T15, T24, T35

  FIT[FIT_SNG]  = EVLN[EVLN_SNG]; //Singlet
  FIT[FIT_GLU]  = EVLN[EVLN_GLU]; //gluon
  FIT[FIT_VAL]  = EVLN[EVLN_VAL]; //valence
  FIT[FIT_V3]   = EVLN[EVLN_V3];  // V3
  FIT[FIT_V8]   = EVLN[EVLN_V8];  // V8
  FIT[FIT_T3]   = EVLN[EVLN_T3];  // T3
  FIT[FIT_T8]   = EVLN[EVLN_T8];  // T8
  FIT[FIT_T15]  = EVLN[EVLN_T15];  // T15

  if (fQED)
    FIT[FIT_GAM] =  EVLN[0];  // photon

  return;
}

real EvolICFitBasis::ComputeSumRules(sumRule rule, int mem, PDFSet *pdf, bool &status) const
{
  // status
  status = false;

  // sum rule calculations
  switch (rule) {
     case SUM_MSR:
       {
         real xsng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real xglu = pdf->IntegratePDF(mem,FIT_GLU,fQ2,PDFSet::XFX,status,fGSLWork);
         real msr = xsng+xglu;
         if (fQED)
           {
             real xgam = pdf->IntegratePDF(mem,FIT_GAM,fQ2,PDFSet::XFX,status,fGSLWork);
             msr += xgam;
           }
         return msr;
       }
       break;
     case SUM_UVL:
       {
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real v3 = pdf->IntegratePDF(mem,FIT_V3,fQ2,PDFSet::FX,status,fGSLWork);
         real v8 = pdf->IntegratePDF(mem,FIT_V8,fQ2,PDFSet::FX,status,fGSLWork);
         real v15 = val;
         return ( 3.0*val + 6.0*v3 + 2.0*v8 + v15 )/12.0;
       }
       break;
     case SUM_DVL:
       {
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real v3 = pdf->IntegratePDF(mem,FIT_V3,fQ2,PDFSet::FX,status,fGSLWork);
         real v8 = pdf->IntegratePDF(mem,FIT_V8,fQ2,PDFSet::FX,status,fGSLWork);
         real v15 = val;
         return ( 3.0*val - 6.0*v3 + 2.0*v8 + v15 )/12.0;
       }
       break;
     case SUM_SVL:
       {
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real v8 = pdf->IntegratePDF(mem,FIT_V8,fQ2,PDFSet::FX,status,fGSLWork);
         real v15 = val;
         return ( 3.0*val - 4.0*v8 + v15 )/12.0;
       }
        break;
     case SUM_CVL:
       {
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real v15 = val;
         return ( val - v15 )/4.0;
       }
        break;

     case SUM_USM:
       {
         real sng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real t3 = pdf->IntegratePDF(mem,FIT_T3,fQ2,PDFSet::XFX,status,fGSLWork);
         real t8 = pdf->IntegratePDF(mem,FIT_T8,fQ2,PDFSet::XFX,status,fGSLWork);
         real t15 = pdf->IntegratePDF(mem,FIT_T15,fQ2,PDFSet::XFX,status,fGSLWork);
         return ( 3.0*sng + 6.0*t3 + 2.0*t8 + t15 )/12.0;
       }
       break;
     case SUM_DSM:
       {
         real sng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real t3 = pdf->IntegratePDF(mem,FIT_T3,fQ2,PDFSet::XFX,status,fGSLWork);
         real t8 = pdf->IntegratePDF(mem,FIT_T8,fQ2,PDFSet::XFX,status,fGSLWork);
         real t15 = pdf->IntegratePDF(mem,FIT_T15,fQ2,PDFSet::XFX,status,fGSLWork);
         return ( 3.0*sng - 6.0*t3 + 2.0*t8 + t15)/12.0;
       }
       break;
     case SUM_SSM:
       {
         real sng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real t8 = pdf->IntegratePDF(mem,FIT_T8,fQ2,PDFSet::XFX,status,fGSLWork);
         real t15 = pdf->IntegratePDF(mem,FIT_T15,fQ2,PDFSet::XFX,status,fGSLWork);
         return ( 3.0*sng - 4.0*t8 + t15)/12.0;
       }
       break;
    case SUM_CSM:
      {
        real sng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
        real t15 = pdf->IntegratePDF(mem,FIT_T15,fQ2,PDFSet::XFX,status,fGSLWork);
        return ( sng - t15)/4.0;
      }
     default:
       cerr << "EvolFitBasis::ComputeSumRules error: unknown sum rule"<<endl;
       exit(-1);
       break;
  }
}


NN30ICFitBasis::NN30ICFitBasis(NNPDFSettings const& nnset):
EvolICFitBasis(nnset)
{
  // PDF Names for plotting
  fPDFNames[NET_SNG] = "Singlet";
  fPDFNames[NET_GLU] = "Gluon";
  fPDFNames[NET_VAL] = "Valence";
  fPDFNames[NET_T3] = "Triplet";
  fPDFNames[NET_DS] = "Sea Asymmetry";
  fPDFNames[NET_SP] = "Strange Sea";
  fPDFNames[NET_SM] = "Strange Valence";
  fPDFNames[NET_CP] = "Charm Sea";
  //fPDFNames[NET_CM] = "Charm Valence";
  if (fQED)
    fPDFNames[NET_GAM] = "Photon";
}

void NN30ICFitBasis::Preprocess(real const& x, real* pdf, PreprocParam const& par)
{
   // Positive definite PDFs
   for (int i = 0; i < fNPDF; i++)
     if (fPDFSqrPos[i])
       pdf[i] *= pdf[i];

   // Transform to evol basis
   real tmppdf[15];
   tmppdf[FIT_SNG] = pdf[NET_SNG];
   tmppdf[FIT_GLU] = pdf[NET_GLU];
   tmppdf[FIT_VAL] = pdf[NET_VAL];
   tmppdf[FIT_V3] = pdf[NET_T3]+2*pdf[NET_DS];
   tmppdf[FIT_V8] = pdf[NET_VAL]-3*pdf[NET_SM]; //-pdf[NET_CM];
   tmppdf[FIT_T3] = pdf[NET_T3];
   tmppdf[FIT_T8] = pdf[NET_SNG]-3*pdf[NET_SP]-pdf[NET_CP];
   tmppdf[FIT_T15] = pdf[NET_SNG]-4*pdf[NET_CP];
   if (fQED)
     tmppdf[FIT_GAM] = pdf[NET_GAM];

   for (int i = 0; i < fNPDF; i++)
     pdf[i] = tmppdf[i];

   // Preprocess in evol basis
   for (int i = 0; i < fNPDF; i++)
     FitBasis::Preprocess(x,i,pdf[i],par);

   return;
}

void NN30ICFitBasis::NetTransform(int const& fl, int const& nfl, int* transform)
{
  for (int i = 0; i < nfl; i++)
    transform[i] = 0;

  switch (fl)
  {
    case FIT_SNG:
    {
      transform[NET_SNG] = 1;
      break;
    }
    case FIT_GLU:
    {
      transform[NET_GLU] = 1;
      break;
    }
    case FIT_VAL:
    {
      transform[NET_VAL] = 1;
      break;
    }
    case FIT_V3:
    {
      transform[NET_T3] = 1;
      transform[NET_DS] = 2;
      break;
    }
    case FIT_V8:
    {
      transform[NET_VAL] = 1;
      transform[NET_SM] = -3;
      //transform[NET_CM] = -1;
      break;
    }
    case FIT_T3:
    {
      transform[NET_T3] = 1;
      break;
    }
    case FIT_T8:
    {
      transform[NET_SNG] = 1;
      transform[NET_SP] = -3;
      transform[NET_CP] = -1;
      break;
    }
    case FIT_T15:
    {
      transform[NET_SNG] = 1;
      transform[NET_CP] = -4;
      break;
    }
    case FIT_GAM:
    {
      transform[NET_GAM] = 1;
      break;
    }
  }
}

/**
 *  Evolution Fit Basis with c+
 **/
NN31ICFitBasis::NN31ICFitBasis(NNPDFSettings const& nnset):
FitBasis(nnset, "NN31ICFitBasis", 8+nnset.IsQED()),
fQED(nnset.IsQED())
{
  // PDF Names for plotting
  fPDFNames[FIT_SNG] = "Singlet";
  fPDFNames[FIT_GLU] = "Gluon";
  fPDFNames[FIT_VAL] = "V";
  fPDFNames[FIT_V3] = "V3";
  fPDFNames[FIT_V8] = "V8";
  fPDFNames[FIT_T3] = "T3";
  fPDFNames[FIT_T8] = "T8";
  fPDFNames[FIT_CP] = "c+";
  if (fQED)
    fPDFNames[FIT_GAM] = "Photon";

  // Damping factor for arclengths
  fArcDampFactor[FIT_SNG] = 1;
  fArcDampFactor[FIT_GLU] = 1;
  fArcDampFactor[FIT_VAL] = 0;
  fArcDampFactor[FIT_V3] = 0;
  fArcDampFactor[FIT_V8] = 0;
  fArcDampFactor[FIT_T3] = 1;
  fArcDampFactor[FIT_T8] = 1;
  fArcDampFactor[FIT_CP] = 1;
  if (fQED)
    fArcDampFactor[FIT_GAM] = 1;
}

/*!
 * \brief NN31ICFitBasis::ComputeParam
 * \param pdf
 * \param mem
 * \param param
 * \param status
 */
void NN31ICFitBasis::ComputeParam(PDFSet* pdf, int mem, PreprocParam& param, bool &status) const
{
  // status
  status = false;

  // Clear old normalisations
  for (int i=0; i<fNPDF; i++)
  {
    param.fPDFNorm[i] = 1.0;
    param.fPDFAux[i] = 0.0;
  }

  // Normalisations pointer
  real* norm = param.fPDFNorm;

  // Quantum number sum rules
  norm[FIT_VAL] = 3.0f/pdf->IntegratePDF(mem,FIT_VAL, fQ2,PDFSet::FX,status,fGSLWork); // Total valence
  norm[FIT_V3] = 1.0f/pdf->IntegratePDF(mem,FIT_V3, fQ2,PDFSet::FX,status,fGSLWork); // V3
  norm[FIT_V8] = 3.0f/pdf->IntegratePDF(mem,FIT_V8, fQ2,PDFSet::FX,status,fGSLWork); // V8

  // ************ QED dependent normalisations **************
  if (fQED)
  {
    norm[FIT_GLU] = (1-pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork) - pdf->IntegratePDF(mem,FIT_GAM,fQ2,PDFSet::XFX,status,fGSLWork))
    / pdf->IntegratePDF(mem,FIT_GLU, fQ2,PDFSet::XFX,status,fGSLWork);
  }
  else
  {
    norm[FIT_GLU] = (1-pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork))/
    pdf->IntegratePDF(mem,FIT_GLU,fQ2,PDFSet::XFX,status,fGSLWork);
  }

  return;

}

/**
 * @brief EvolFitBasis::BASIS2EVLN
 * @param FIT
 * @param EVLN
 */
void NN31ICFitBasis::BASIS2EVLN(const real *FIT, real *EVLN) const
{
  if ( fQED )
    EVLN[EVLN_GAM] = FIT[FIT_GAM];
  else
    EVLN[EVLN_GAM] = 0;

  EVLN[EVLN_SNG]  = FIT[FIT_SNG]; //Singlet
  EVLN[EVLN_GLU]  = FIT[FIT_GLU]; //Gluon
  EVLN[EVLN_VAL]  = FIT[FIT_VAL]; //Valence
  EVLN[EVLN_V3]   = FIT[FIT_V3];  //V3
  EVLN[EVLN_V8]   = FIT[FIT_V8];  //V8
  EVLN[EVLN_V15]  = FIT[FIT_VAL]; //V15 = V (c- = 0)
  EVLN[EVLN_V24]  = FIT[FIT_VAL]; //V24 = V
  EVLN[EVLN_V35]  = FIT[FIT_VAL]; //V35 = V
  EVLN[EVLN_T3]   = FIT[FIT_T3];  //T3
  EVLN[EVLN_T8]   = FIT[FIT_T8];  //T8
  EVLN[EVLN_T15]  = FIT[FIT_SNG] - 4*FIT[FIT_CP]; //T15
  EVLN[EVLN_T24]  = FIT[FIT_SNG]; //T24 = S
  EVLN[EVLN_T35]  = FIT[FIT_SNG]; //T35 = S

  return;
}

/*!
 * \brief NN31ICFitBasis::EVLN2BASIS
 * \param EVLN
 * \param FIT
 */
void NN31ICFitBasis::EVLN2BASIS(const real *EVLN, real *FIT) const
{
  // Order in fitting basis
  // S g V V3 V8 T3 T8 c+ gam

  // Order in Evln bassi
  // γ, Σ, g, V, V3, V8, V15, V24, V35, T3, T8, T15, T24, T35

  FIT[FIT_SNG]  = EVLN[EVLN_SNG]; //Singlet
  FIT[FIT_GLU]  = EVLN[EVLN_GLU]; //gluon
  FIT[FIT_VAL]  = EVLN[EVLN_VAL]; //valence
  FIT[FIT_V3]   = EVLN[EVLN_V3];  // V3
  FIT[FIT_V8]   = EVLN[EVLN_V8];  // V8
  FIT[FIT_T3]   = EVLN[EVLN_T3];  // T3
  FIT[FIT_T8]   = EVLN[EVLN_T8];  // T8
  FIT[FIT_CP]   = (EVLN[EVLN_SNG]-EVLN[EVLN_T15])/4;  // T15

  if (fQED)
    FIT[FIT_GAM] =  EVLN[0];  // photon

  return;
}

/*!
 * \brief NN31ICFitBasis::ComputeSumRules
 * \param rule
 * \param mem
 * \param pdf
 * \param status
 * \return
 */
real NN31ICFitBasis::ComputeSumRules(sumRule rule, int mem, PDFSet *pdf, bool &status) const
{
  // status
  status = false;

  // sum rule calculations
  switch (rule) {
     case SUM_MSR:
       {
         real xsng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real xglu = pdf->IntegratePDF(mem,FIT_GLU,fQ2,PDFSet::XFX,status,fGSLWork);
         real msr = xsng+xglu;
         if (fQED)
           {
             real xgam = pdf->IntegratePDF(mem,FIT_GAM,fQ2,PDFSet::XFX,status,fGSLWork);
             msr += xgam;
           }
         return msr;
       }
       break;
     case SUM_UVL:
       {
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real v3 = pdf->IntegratePDF(mem,FIT_V3,fQ2,PDFSet::FX,status,fGSLWork);
         real v8 = pdf->IntegratePDF(mem,FIT_V8,fQ2,PDFSet::FX,status,fGSLWork);
         real v15 = val;
         return ( 3.0*val + 6.0*v3 + 2.0*v8 + v15 )/12.0;
       }
       break;
     case SUM_DVL:
       {
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real v3 = pdf->IntegratePDF(mem,FIT_V3,fQ2,PDFSet::FX,status,fGSLWork);
         real v8 = pdf->IntegratePDF(mem,FIT_V8,fQ2,PDFSet::FX,status,fGSLWork);
         real v15 = val;
         return ( 3.0*val - 6.0*v3 + 2.0*v8 + v15 )/12.0;
       }
       break;
     case SUM_SVL:
       {
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real v8 = pdf->IntegratePDF(mem,FIT_V8,fQ2,PDFSet::FX,status,fGSLWork);
         real v15 = val;
         return ( 3.0*val - 4.0*v8 + v15 )/12.0;
       }
        break;
     case SUM_CVL:
       {
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real v15 = val;
         return ( val - v15 )/4.0;
       }
        break;

     case SUM_USM:
       {
         real sng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real t3 = pdf->IntegratePDF(mem,FIT_T3,fQ2,PDFSet::XFX,status,fGSLWork);
         real t8 = pdf->IntegratePDF(mem,FIT_T8,fQ2,PDFSet::XFX,status,fGSLWork);
         real cp = pdf->IntegratePDF(mem,FIT_CP,fQ2,PDFSet::XFX,status,fGSLWork);
         real t15 = sng-4*cp;
         return ( 3.0*sng + 6.0*t3 + 2.0*t8 + t15 )/12.0;
       }
       break;
     case SUM_DSM:
       {
         real sng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real t3 = pdf->IntegratePDF(mem,FIT_T3,fQ2,PDFSet::XFX,status,fGSLWork);
         real t8 = pdf->IntegratePDF(mem,FIT_T8,fQ2,PDFSet::XFX,status,fGSLWork);
         real cp = pdf->IntegratePDF(mem,FIT_CP,fQ2,PDFSet::XFX,status,fGSLWork);
         real t15 = sng-4*cp;
         return ( 3.0*sng - 6.0*t3 + 2.0*t8 + t15)/12.0;
       }
       break;
     case SUM_SSM:
       {
         real sng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real t8 = pdf->IntegratePDF(mem,FIT_T8,fQ2,PDFSet::XFX,status,fGSLWork);
         real cp = pdf->IntegratePDF(mem,FIT_CP,fQ2,PDFSet::XFX,status,fGSLWork);
         real t15 = sng-4*cp;
         return ( 3.0*sng - 4.0*t8 + t15)/12.0;
       }
       break;
    case SUM_CSM:
      {
        real cp = pdf->IntegratePDF(mem,FIT_CP,fQ2,PDFSet::XFX,status,fGSLWork);
        return cp;
      }
      break;
     default:
       cerr << "NN31ICFitBasis::ComputeSumRules error: unknown sum rule"<<endl;
       exit(-1);
       break;
  }
}

// ***********************************************************************************

/**
 *  Evolution Fit Basis with c+
 **/
NoSumRuleBasis::NoSumRuleBasis(NNPDFSettings const& nnset):
FitBasis(nnset, "NoSumRuleBasis", 8)
{
  // PDF Names for plotting
  fPDFNames[FIT_SNG] = "Singlet";
  fPDFNames[FIT_GLU] = "Gluon";
  fPDFNames[FIT_VAL] = "V";
  fPDFNames[FIT_V3] = "V3";
  fPDFNames[FIT_V8] = "V8";
  fPDFNames[FIT_T3] = "T3";
  fPDFNames[FIT_T8] = "T8";
  fPDFNames[FIT_CP] = "c+";

  // Damping factor for arclengths
  fArcDampFactor[FIT_SNG] = 1;
  fArcDampFactor[FIT_GLU] = 1;
  fArcDampFactor[FIT_VAL] = 0;
  fArcDampFactor[FIT_V3] = 0;
  fArcDampFactor[FIT_V8] = 0;
  fArcDampFactor[FIT_T3] = 1;
  fArcDampFactor[FIT_T8] = 1;
  fArcDampFactor[FIT_CP] = 1;
}

/*!
 * \brief NoSumRuleBasis::ComputeParam
 * \param pdf
 * \param mem
 * \param param
 * \param status
 */
void NoSumRuleBasis::ComputeParam(PDFSet*, int, PreprocParam& param, bool &) const
{
  // Clear old normalisations
  for (int i=0; i<fNPDF; i++)
  {
    param.fPDFNorm[i] = 1.0;
    param.fPDFAux[i] = 0.0;
  }
  return;
}

/**
 * @brief EvolFitBasis::BASIS2EVLN
 * @param FIT
 * @param EVLN
 */
void NoSumRuleBasis::BASIS2EVLN(const real *FIT, real *EVLN) const
{
  EVLN[EVLN_GAM] = 0;
  EVLN[EVLN_SNG]  = FIT[FIT_SNG]; //Singlet
  EVLN[EVLN_GLU]  = FIT[FIT_GLU]; //Gluon
  EVLN[EVLN_VAL]  = FIT[FIT_VAL]; //Valence
  EVLN[EVLN_V3]   = FIT[FIT_V3];  //V3
  EVLN[EVLN_V8]   = FIT[FIT_V8];  //V8
  EVLN[EVLN_V15]  = FIT[FIT_VAL]; //V15 = V (c- = 0)
  EVLN[EVLN_V24]  = FIT[FIT_VAL]; //V24 = V
  EVLN[EVLN_V35]  = FIT[FIT_VAL]; //V35 = V
  EVLN[EVLN_T3]   = FIT[FIT_T3];  //T3
  EVLN[EVLN_T8]   = FIT[FIT_T8];  //T8
  EVLN[EVLN_T15]  = FIT[FIT_SNG] - 4*FIT[FIT_CP]; //T15
  EVLN[EVLN_T24]  = FIT[FIT_SNG]; //T24 = S
  EVLN[EVLN_T35]  = FIT[FIT_SNG]; //T35 = S

  return;
}

/*!
 * \brief NoSumRuleBasis::EVLN2BASIS
 * \param EVLN
 * \param FIT
 */
void NoSumRuleBasis::EVLN2BASIS(const real *EVLN, real *FIT) const
{
  // Order in fitting basis
  // S g V V3 V8 T3 T8 c+ gam

  // Order in Evln bassi
  // γ, Σ, g, V, V3, V8, V15, V24, V35, T3, T8, T15, T24, T35

  FIT[FIT_SNG]  = EVLN[EVLN_SNG]; //Singlet
  FIT[FIT_GLU]  = EVLN[EVLN_GLU]; //gluon
  FIT[FIT_VAL]  = EVLN[EVLN_VAL]; //valence
  FIT[FIT_V3]   = EVLN[EVLN_V3];  // V3
  FIT[FIT_V8]   = EVLN[EVLN_V8];  // V8
  FIT[FIT_T3]   = EVLN[EVLN_T3];  // T3
  FIT[FIT_T8]   = EVLN[EVLN_T8];  // T8
  FIT[FIT_CP]   = (EVLN[EVLN_SNG]-EVLN[EVLN_T15])/4;  // T15

  return;
}

/*!
 * \brief NoSumRuleBasis::ComputeSumRules
 * \param rule
 * \param mem
 * \param pdf
 * \param status
 * \return
 */
real NoSumRuleBasis::ComputeSumRules(sumRule rule, int mem, PDFSet *pdf, bool &status) const
{
  // status
  status = false;

  // sum rule calculations
  switch (rule) {
     case SUM_MSR:
       {
         real xsng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real xglu = pdf->IntegratePDF(mem,FIT_GLU,fQ2,PDFSet::XFX,status,fGSLWork);
         real msr = xsng+xglu;
         return msr;
       }
       break;
     case SUM_UVL:
       {
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real v3 = pdf->IntegratePDF(mem,FIT_V3,fQ2,PDFSet::FX,status,fGSLWork);
         real v8 = pdf->IntegratePDF(mem,FIT_V8,fQ2,PDFSet::FX,status,fGSLWork);
         real v15 = val;
         return ( 3.0*val + 6.0*v3 + 2.0*v8 + v15 )/12.0;
       }
       break;
     case SUM_DVL:
       {
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real v3 = pdf->IntegratePDF(mem,FIT_V3,fQ2,PDFSet::FX,status,fGSLWork);
         real v8 = pdf->IntegratePDF(mem,FIT_V8,fQ2,PDFSet::FX,status,fGSLWork);
         real v15 = val;
         return ( 3.0*val - 6.0*v3 + 2.0*v8 + v15 )/12.0;
       }
       break;
     case SUM_SVL:
       {
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real v8 = pdf->IntegratePDF(mem,FIT_V8,fQ2,PDFSet::FX,status,fGSLWork);
         real v15 = val;
         return ( 3.0*val - 4.0*v8 + v15 )/12.0;
       }
        break;
     case SUM_CVL:
       {
         real val = pdf->IntegratePDF(mem,FIT_VAL,fQ2,PDFSet::FX,status,fGSLWork);
         real v15 = val;
         return ( val - v15 )/4.0;
       }
        break;

     case SUM_USM:
       {
         real sng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real t3 = pdf->IntegratePDF(mem,FIT_T3,fQ2,PDFSet::XFX,status,fGSLWork);
         real t8 = pdf->IntegratePDF(mem,FIT_T8,fQ2,PDFSet::XFX,status,fGSLWork);
         real cp = pdf->IntegratePDF(mem,FIT_CP,fQ2,PDFSet::XFX,status,fGSLWork);
         real t15 = sng-4*cp;
         return ( 3.0*sng + 6.0*t3 + 2.0*t8 + t15 )/12.0;
       }
       break;
     case SUM_DSM:
       {
         real sng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real t3 = pdf->IntegratePDF(mem,FIT_T3,fQ2,PDFSet::XFX,status,fGSLWork);
         real t8 = pdf->IntegratePDF(mem,FIT_T8,fQ2,PDFSet::XFX,status,fGSLWork);
         real cp = pdf->IntegratePDF(mem,FIT_CP,fQ2,PDFSet::XFX,status,fGSLWork);
         real t15 = sng-4*cp;
         return ( 3.0*sng - 6.0*t3 + 2.0*t8 + t15)/12.0;
       }
       break;
     case SUM_SSM:
       {
         real sng = pdf->IntegratePDF(mem,FIT_SNG,fQ2,PDFSet::XFX,status,fGSLWork);
         real t8 = pdf->IntegratePDF(mem,FIT_T8,fQ2,PDFSet::XFX,status,fGSLWork);
         real cp = pdf->IntegratePDF(mem,FIT_CP,fQ2,PDFSet::XFX,status,fGSLWork);
         real t15 = sng-4*cp;
         return ( 3.0*sng - 4.0*t8 + t15)/12.0;
       }
       break;
    case SUM_CSM:
      {
        real cp = pdf->IntegratePDF(mem,FIT_CP,fQ2,PDFSet::XFX,status,fGSLWork);
        return cp;
      }
      break;
     default:
       cerr << "NoSumRuleBasis::ComputeSumRules error: unknown sum rule"<<endl;
       exit(-1);
       break;
  }
}
