// $Id: lhapdfset.h 3177 2015-08-18 14:43:31Z stefano.carrazza@mi.infn.it $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once
#include "pdfset.h"
#include "LHAPDF/LHAPDF.h"

#include <vector>

namespace NNPDF
{

  /**
   *  \class LHAPDFSet
   *  \brief Class for handle LHAPDF methods
   */
  class LHAPDFSet : private LHAPDF::PDFSet, public NNPDF::PDFSet
  {
  public:
    //!< Constructor
    LHAPDFSet(std::string const&, erType);
    virtual ~LHAPDFSet();   //!< Destructor

    //!< Get PDF array in evolution basis
    virtual void GetPDF(real const& x, real const& Q2, int const& n, real* pdf) const;
    
    //!< Check if flavor exists
    bool hasFlavor(int pdgid) const;
    
    //!< Return single flavor
    real xfxQ(real const& x, real const& Q, int const& n, int const& fl) const;

  protected:
    // Member PDFs
    std::vector<LHAPDF::PDF*> fMemberPDFs;

    // Internal transformation array
    real* fLHA;
  };
}
