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
#include <memory>
#include "NNPDF/common.h"
#include "NNPDF/fastkernel.h"

namespace NNPDF
{
    /**
    * \class FKGenerator
    * \brief Class for filling FastKernel tables
    */
    class FKGenerator : public FKTable
    {
      private:
        FKGenerator();                          //!< Disable default constructor
        FKGenerator( const FKGenerator& );      //!< Disable default copy-constructor
        FKGenerator& operator=(const FKGenerator&); //!< Disable copy-assignment

		std::unique_ptr<double[]> fAccumulator;

      public:
        FKGenerator( std::istream& ); // Constructor
        virtual ~FKGenerator() = default;

        void Finalise();    //!< Copy Accumulator to fSigma

         // Fill the FKTable
        void Fill(  int const& d,     // Datapoint index
                    int const& ix1,   // First x-index
                    int const& ix2,   // Second x-index
                    int const& ifl1,  // First flavour index
                    int const& ifl2,  // Second flavour index
                    real const& isig  // FK Value
                  );

        // DIS version of Fill
        void Fill(  int const& d,     // Datapoint index
                    int const& ix,    // x-index
                    int const& ifl,   // flavour index
                    real const& isig  // FK Value
                  );
    };

}
