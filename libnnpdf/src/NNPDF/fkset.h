// $Id$
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

// libnnpdf 2014 nph

#pragma once

#include <string>
#include <vector>

#include "common.h"
#include "fastkernel.h"

namespace NNPDF
{
    class FKTable;

    // Operator for compound FK tables
    typedef void (*SigmaOp) (int const&, std::vector<real*> const&, real* );

    /**
    * \class FKSet
    * \brief Class for holding sets of FK tables, alongside information on compound observable computation
    */
    class FKSet
    {
      public:
        FKSet(SigmaOp, std::vector<FKTable*> const&); //!< Constructor
        virtual ~FKSet();                             //!< Destructor

        FKSet(FKSet const&); //!< Copy-constructor
        FKSet(FKSet const&, std::vector<int> const&); //!< Masked copy-constructor

        SigmaOp GetOperator()        const { return *fOperator;};

        int     const&  GetNSigma()  const { return fNSigma;   };
        int     const&  GetNDataFK() const { return fNDataFK; };
        bool    const&  IsHadronic() const { return fHadronic; };

        std::string const& GetDataName() const {return fDataName;};
        const FKTable* GetFK(size_t const& i) const {return fFK[i];}; //!< Fetch FK tables

        // parse operator string
        static SigmaOp parseOperator(std::string const& op);

      private:
        FKSet();                          //disable default constructor
        FKSet& operator=(const FKSet&);   //disable copy-assignment

        const SigmaOp fOperator;  //!< Operator acting on combination observables
        const int     fNSigma;    //!< Number of FastKernel Grids for this dataset
        const int     fNDataFK;   //!< Number of datapoints
        const bool    fHadronic;  //!< Hadronic Observables

        const std::string fDataName;    //!< Name of dataset
        FKTable **fFK;            //!< FastKernel tables
    };

}
