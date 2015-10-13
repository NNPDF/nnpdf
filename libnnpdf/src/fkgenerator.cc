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

namespace NNPDF
{

  FKGenerator::FKGenerator( std::istream& is ):
  FKTable(is) 
  {

  };

  FKGenerator::~FKGenerator()
  {

  }

  // Fill the FKTable - beware there is no error checking performed!
  void FKGenerator::Fill(   size_t const& d,     // Datapoint index
                            size_t const& ix1,   // First x-index
                            size_t const& ix2,   // Second x-index
                            size_t const& ifl1,  // First flavour index
                            size_t const& ifl2,  // Second flavour index
                            real const& fk  // FK Value
                          )
  {

    if (fk==0)
      return;

    if (d >= fNData)
    {
      std::cerr << "FKGenerator::Fill Error: datapoint "<<d <<" out of bounds."<<std::endl;
      exit(-1);
    }

    if (ix1 >= fNx)
    {
      std::cerr << "FKGenerator::Fill Error: xpoint "<<ix1 <<" out of bounds."<<std::endl;
      exit(-1);
    }

    if (ix2 >= fNx)
    {
      std::cerr << "FKGenerator::Fill Error: xpoint "<<ix2 <<" out of bounds."<<std::endl;
      exit(-1);
    }

    if (ifl1 >= 14)
    {
      std::cerr << "FKGenerator::Fill Error: flavour "<<ifl1 <<" out of bounds."<<std::endl;
      exit(-1);
    }

    if (ifl2 >= 14)
    {
      std::cerr << "FKGenerator::Fill Error: flavour "<<ifl2 <<" out of bounds."<<std::endl;
      exit(-1);
    }

    // pointer to FKTable segment
    const size_t jFL = 14*ifl1 + ifl2;
    const size_t iSig = d*fDSz+jFL*fTx+ix1*fNx+ix2 ;

    // Assign FK Table
    fSigma[iSig] += fk;

    return;
  };

  // DIS version of Fill
  void FKGenerator::Fill( size_t const& d,     // Datapoint index
                          size_t const& ix,    // x-index
                          size_t const& ifl,   // flavour index
                          real const& fk  // FK Value
                        )
  {

    if (fk==0)
      return;

    if (d >= fNData)
    {
      std::cerr << "FKGenerator::Fill Error: datapoint "<<d <<" out of bounds."<<std::endl;
      exit(-1);
    }

    if (ix >= fNx)
    {
      std::cerr << "FKGenerator::Fill Error: xpoint "<<ix <<" out of bounds."<<std::endl;
      exit(-1);
    }

    if (ifl >= 14)
    {
      std::cerr << "FKGenerator::Fill Error: flavour "<<ifl <<" out of bounds."<<std::endl;
      exit(-1);
    }

    // pointer to FKTable segment
    const size_t iSig = d*fDSz+ifl*fNx+ix;

    // Assign FK Table
    fSigma[iSig] += fk;

    return;
  };

}
