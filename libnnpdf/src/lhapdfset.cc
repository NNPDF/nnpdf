// $Id: lhapdfset.cc 3177 2015-08-18 14:43:31Z stefano.carrazza@mi.infn.it $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#include <cmath>
#include <cstdlib>
#include <algorithm>

#include "LHAPDF/LHAPDF.h"
#include "NNPDF/lhapdfset.h"
#include "NNPDF/pdfset.h"
#include "NNPDF/utils.h"
#include "NNPDF/exceptions.h"

using namespace NNPDF;

/**
 * The constructor
 */
LHAPDFSet::LHAPDFSet(std::string const& pdfname, erType etype):
LHAPDF::PDFSet(pdfname),
NNPDF::PDFSet(pdfname, size(), etype),
fLHA(new real[14])
{
  // Verify erType
  const std::string LHError = errorType();
  switch (etype)
  {
     case erType::ER_NONE:
     case erType::ER_MCT0:
      fMembers = 1;
      break;

     case erType::ER_MC:
     case erType::ER_MC68:
      if (LHError.compare("replicas") != 0 )
        throw RuntimeException("LHAPDFSet::LHAPDFSet", "Error: ErrorSet Types do not match: ER_MC(68) and " + LHError);
      fMembers-=1; // Remove replica zero
      break;

     case erType::ER_EIG:
     case erType::ER_EIG90:
      if (LHError.compare("hessian") != 0 )
        throw RuntimeException("LHAPDFSet::LHAPDFSet", "Error: ErrorSet Types do not match: ER_EIG(90) and " + LHError);
      break;

     case erType::ER_SYMEIG:
      if (LHError.compare("symmhessian") != 0 )
        throw RuntimeException("LHAPDFSet::LHAPDFSet", "Error: ErrorSet Types do not match: ER_SYMEIG and " + LHError);
      break;
  }

  // Load member PDFs
  if (fMembers == 1)
    fMemberPDFs.push_back(mkPDF(0));
  else
    mkPDFs( fMemberPDFs );

    get_logger() << pdfname<< " Initialised with " << fMembers<<" members and errorType "<<LHError<<std::endl;
}

LHAPDFSet::LHAPDFSet(const std::string & pdfname, const int &replica):
  LHAPDF::PDFSet(pdfname),
  NNPDF::PDFSet(pdfname, 1, NNPDF::PDFSet::erType::ER_NONE),
  fLHA(new real[14])
{
  fMemberPDFs.push_back(mkPDF(replica));
  get_logger() << pdfname << " Initialised with " << fMembers <<" members and no errorType." << std::endl;
}

LHAPDFSet::~LHAPDFSet()
{
  delete[] fLHA;

  for (size_t i=0; i<fMemberPDFs.size(); i++)
    delete fMemberPDFs[i];
  fMemberPDFs.clear();
};

bool LHAPDFSet::hasFlavor(int pdgid) const
{
  return fMemberPDFs[0]->hasFlavor(pdgid);
}

void LHAPDFSet::GetPDF(real const& x, real const& Q2, int const& n, real* pdf)      const
{
  // Skip replica 0 for MC like ensembles
  const int iMem = (fEtype == erType::ER_MC || fEtype == erType::ER_MC68) ? n+1:n;

  if (!fMemberPDFs[iMem]->inPhysicalRangeXQ2(x, Q2))
    throw RangeError("LHAPDFSet::GetPDF","kinematics out of set range: x="+ std::to_string(x) +" Q2="+ std::to_string(Q2));

  // Clear LHA (probably unneeded)
  for (int i=0; i<14; i++)
    fLHA[i] = 0.0;

   //** PDG CODES **//
  // Let's be super explicit here

  fLHA[TBAR] = fMemberPDFs[iMem]->xfxQ2(-6, x, Q2);
  fLHA[BBAR] = fMemberPDFs[iMem]->xfxQ2(-5, x, Q2);
  fLHA[CBAR] = fMemberPDFs[iMem]->xfxQ2(-4, x, Q2);
  fLHA[SBAR] = fMemberPDFs[iMem]->xfxQ2(-3, x, Q2);
  fLHA[UBAR] = fMemberPDFs[iMem]->xfxQ2(-2, x, Q2);
  fLHA[DBAR] = fMemberPDFs[iMem]->xfxQ2(-1, x, Q2);

  fLHA[T] = fMemberPDFs[iMem]->xfxQ2(6, x, Q2);
  fLHA[B] = fMemberPDFs[iMem]->xfxQ2(5, x, Q2);
  fLHA[C] = fMemberPDFs[iMem]->xfxQ2(4, x, Q2);
  fLHA[S] = fMemberPDFs[iMem]->xfxQ2(3, x, Q2);
  fLHA[U] = fMemberPDFs[iMem]->xfxQ2(2, x, Q2);
  fLHA[D] = fMemberPDFs[iMem]->xfxQ2(1, x, Q2);

  fLHA[GLUON] = fMemberPDFs[iMem]->xfxQ2(21, x, Q2);
  fLHA[PHT]   = fMemberPDFs[iMem]->xfxQ2(22, x, Q2);

  // Transform
  LHA2EVLN(fLHA, pdf);
}

real LHAPDFSet::xfxQ(real const& x, real const& Q, int const& n, int const& fl) const
{
  // Skip replica 0 for MC like ensembles
  const int iMem = (fEtype == erType::ER_MC || fEtype == erType::ER_MC68) ? n+1:n;

  if (!fMemberPDFs[iMem]->inPhysicalRangeXQ(x, Q))
    throw RangeError("LHAPDFSet::GetPDF", "kinematics out of set range: x=" + std::to_string(x) + " Q=" + std::to_string(Q));

  return fMemberPDFs[iMem]->xfxQ(fl, x, Q);
}




