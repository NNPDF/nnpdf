// $Id: pdfset.h 3177 2015-08-18 14:43:31Z stefano.carrazza@mi.infn.it $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include <string>
#include "common.h"

#include <gsl/gsl_integration.h>

namespace NNPDF
{
    /**
     * \class PDFSet
     * \brief An abstract class that provides the basic interface to PDF sets
     */

    class PDFSet
    {
      public:
           //! Error types (Etype)
        enum erType 
        { ER_NONE,  //!< PDF set is not an error set
          ER_MC,    //!< 1sigma error for NNPDF
          ER_MC68,  //!< 68cl error for NNPDF
          ER_MCT0,  //!< NNPDF set use only replica 0 for T0
          ER_EIG,   //!< 1sigma error for CTEQ & MSTW
          ER_EIG90, //!< 90cl error for CTEQ & MSTW
          ER_SYMEIG //!< Symmetric eigenvectors for PDF4LHC
        };

        // LHA-style flavour basis
        enum {TBAR,BBAR,CBAR,SBAR,UBAR,DBAR,GLUON,D,U,S,C,B,T,PHT};
        // NNPDF-style EVLN basis
        enum evlnBasis {  EVLN_GAM, EVLN_SNG, EVLN_GLU, EVLN_VAL, EVLN_V3,
                          EVLN_V8, EVLN_V15, EVLN_V24, EVLN_V35,  EVLN_T3, 
                          EVLN_T8, EVLN_T15, EVLN_T24, EVLN_T35 };

      protected:
        PDFSet(std::string const&, int const& fMembers, erType const& etype ); //!< Constructor
      
        const std::string fSetName;   //!< Name of the PDF Set

        int fMembers;           //!< Number of members in the PDF set
        const erType fEtype;    //!< Type of error handling for the set

      public:
        virtual ~PDFSet();

        // Verbosity control
        static bool Verbose;

        //! Configuration of PDFs
        virtual real GetPDF(real const& x, real const& Q2, int const& n, int const& fl) const;   //!< Get single PDF
        virtual void GetPDF(real const& x, real const& Q2, int const& n, real* pdf)     const = 0;   //!< Get PDF array

        // Get methods
        const std::string& GetSetName() const { return fSetName; }  //!< Return the set name
        const erType& GetEtype()        const { return fEtype;   }  //!< Fetch the error type of the PDF set
        const int& GetMembers()         const { return fMembers; }  //!< Return number of memeber PDFs

        enum intType {  XFX = true,  FX = false  }; // Integration type helper
        // Integrate the PDF
        real IntegratePDF(  int const& mem, int const& fl, real const& Q2, 
                            PDFSet::intType xfx, bool& gslerror, 
                            gsl_integration_workspace *,
                            real xmin = 0.0,real xmax = 1.0) const; 

        static std::string errString(erType const& type)
        {
            switch (type)
            {
                case ER_NONE:
                    return "No Errors";
                case ER_MC:
                    return "Monte Carlo";
                case ER_MC68:
                    return "Monte Carlo 68pc";
                case ER_MCT0:
                    return "Monte Carlo t0";
                case ER_EIG:
                    return "Eigenvector 68pc";
                case ER_EIG90:
                    return "Eigenvector 90pc";
                case ER_SYMEIG:
                    return "Symmetric eigenvector";
            }

            return "Unrecognised Error Type";
        }

	// Transformations
	static void LHA2EVLN(const real *LHA, real *EVLN);
	static void EVLN2LHA(const real* EVL, real* LHA);
    };

}
