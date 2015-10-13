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
    /**
    * \class FKGenerator Class
    * \brief Class for filling FastKernel tables
    */
    class FKGenerator : public FKTable
    {
      private:
        FKGenerator();                          //!< Disable default constructor
        FKGenerator( const FKGenerator& );      //!< Disable default copy-constructor
        FKGenerator& operator=(const FKGenerator&); //!< Disable copy-assignment

      public:
        FKGenerator( std::istream& ); // Constructor
        virtual ~FKGenerator(); //!< Destructor

         // Fill the FKTable
        void Fill(  size_t const& d,     // Datapoint index
                    size_t const& ix1,   // First x-index
                    size_t const& ix2,   // Second x-index
                    size_t const& ifl1,  // First flavour index
                    size_t const& ifl2,  // Second flavour index
                    real const& isig  // FK Value
                  );

        // DIS version of Fill
        void Fill(  size_t const& d,     // Datapoint index
                    size_t const& ix,    // x-index
                    size_t const& ifl,   // flavour index
                    real const& isig  // FK Value
                  );
    };

}