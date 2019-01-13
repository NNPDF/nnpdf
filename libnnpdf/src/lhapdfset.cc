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
NNPDF::PDFSet(pdfname, size(), etype)
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
  NNPDF::PDFSet(pdfname, 1, NNPDF::PDFSet::erType::ER_NONE)
{
  fMemberPDFs.push_back(mkPDF(replica));
  get_logger() << pdfname << " Initialised with " << fMembers <<" members and no errorType." << std::endl;
}

LHAPDFSet::~LHAPDFSet()
{
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

   //** PDG CODES **//
  // Let's be super explicit here
  std::array<real, 14> LHA;

  LHA[TBAR] = fMemberPDFs[iMem]->xfxQ2(-6, x, Q2);
  LHA[BBAR] = fMemberPDFs[iMem]->xfxQ2(-5, x, Q2);
  LHA[CBAR] = fMemberPDFs[iMem]->xfxQ2(-4, x, Q2);
  LHA[SBAR] = fMemberPDFs[iMem]->xfxQ2(-3, x, Q2);
  LHA[UBAR] = fMemberPDFs[iMem]->xfxQ2(-2, x, Q2);
  LHA[DBAR] = fMemberPDFs[iMem]->xfxQ2(-1, x, Q2);

  LHA[T] = fMemberPDFs[iMem]->xfxQ2(6, x, Q2);
  LHA[B] = fMemberPDFs[iMem]->xfxQ2(5, x, Q2);
  LHA[C] = fMemberPDFs[iMem]->xfxQ2(4, x, Q2);
  LHA[S] = fMemberPDFs[iMem]->xfxQ2(3, x, Q2);
  LHA[U] = fMemberPDFs[iMem]->xfxQ2(2, x, Q2);
  LHA[D] = fMemberPDFs[iMem]->xfxQ2(1, x, Q2);

  LHA[GLUON] = fMemberPDFs[iMem]->xfxQ2(21, x, Q2);
  LHA[PHT]   = fMemberPDFs[iMem]->xfxQ2(22, x, Q2);

  // Transform
  LHA2EVLN(LHA.data(), pdf);
}

real LHAPDFSet::xfxQ(real const& x, real const& Q, int const& n, int const& fl) const
{
  // Skip replica 0 for MC like ensembles
  const int iMem = (fEtype == erType::ER_MC || fEtype == erType::ER_MC68) ? n+1:n;

  if (!fMemberPDFs[iMem]->inPhysicalRangeXQ(x, Q))
    throw RangeError("LHAPDFSet::GetPDF", "kinematics out of set range: x=" + std::to_string(x) + " Q=" + std::to_string(Q));

  return fMemberPDFs[iMem]->xfxQ(fl, x, Q);
}




