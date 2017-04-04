// $Id$
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#include <stdio.h>
#include <iostream>
#include <cstdlib>

#include "NNPDF/fkgenerator.h"
#include "NNPDF/utils.h"
#include "NNPDF/exceptions.h"

namespace NNPDF
{

  FKGenerator::FKGenerator( std::istream& is ):
  FKTable(is),
  fAccumulator(new double[fNData*fDSz]())
  {
  }

  // Fill the FKTable - beware there is no error checking performed!
  void FKGenerator::Fill(   int const& d,     // Datapoint index
                            int const& ix1,   // First x-index
                            int const& ix2,   // Second x-index
                            int const& ifl1,  // First flavour index
                            int const& ifl2,  // Second flavour index
                            real const& fk  // FK Value
                          )
  {

    if (fk==0)
      return;

    if (d >= fNData)
      throw RangeError("FKGenerator::Fill","datapoint " + std::to_string(d) + "out of bounds.");

    if (ix1 >= fNx)
      throw RangeError("FKGenerator::Fill","xpoint " + std::to_string(ix1) + " out of bounds.");

    if (ix2 >= fNx)
      throw RangeError("FKGenerator::Fill","xpoint " + std::to_string(ix2) + " out of bounds.");

    if (ifl1 >= 14)
      throw RangeError("FKGenerator::Fill","flavour " + std::to_string(ifl1) + " out of bounds.");

    if (ifl2 >= 14)
      throw RangeError("FKGenerator::Fill","flavour " + std::to_string(ifl2) + " out of bounds.");

    // pointer to FKTable segment
    const size_t jFL = 14*ifl1 + ifl2;
    const size_t iSig = d*fDSz+jFL*fTx+ix1*fNx+ix2 ;

    // Assign FK Table
    fAccumulator[iSig] += fk;

    return;
  };

  // DIS version of Fill
  void FKGenerator::Fill( int const& d,     // Datapoint index
                          int const& ix,    // x-index
                          int const& ifl,   // flavour index
                          real const& fk  // FK Value
                        )
  {

    if (fk==0)
      return;

    if (d >= fNData)
      throw RangeError("FKGenerator::Fill","datapoint " + std::to_string(d) + " out of bounds.");

    if (ix >= fNx)
      throw RangeError("FKGenerator::Fill","xpoint " + std::to_string(ix) + " out of bounds.");

    if (ifl >= 14)
      throw RangeError("FKGenerator::Fill","flavour " + std::to_string(ifl) + " out of bounds.");

    // pointer to FKTable segment
    const size_t iSig = d*fDSz+ifl*fNx+ix;

    // Assign FK Table
    fAccumulator[iSig] += fk;

    return;
  };

  void FKGenerator::Finalise()
  {
    for (int i=0; i<fNData*fDSz; i++)
      fSigma[i] = fAccumulator[i];
  }

}
