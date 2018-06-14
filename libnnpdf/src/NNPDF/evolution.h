// $Id$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
#pragma once

#include <NNPDF/pdfset.h>
#include <vector>

using std::array;
using std::vector;

namespace NNPDF {

  /**
   *  \class EvolutionSubGrid
   *  \brief Performs PDF evolution from an initial scale to a grid of kinematics.
   *
   *  This class stores EvolutionOperator objects used to perform PDF evolution.
   *  It also provides convenience methods for performing said evolution.
   *
   *  Each subgrid spans scales between HQ thresholds. Therefore each subgrid
   *  outputs a fixed number of output flavours.
   *
   *  LHAPDF permits subgrids to have different x-grids. While we don't use that here,
   *  it means we don't have to worry too much about keeping x-grids synchronised between
   *  subgrids.
   */

    class EvolutionSubGrid
    {
        public:
            EvolutionSubGrid(); //!< Grid constructor

            // Evolve an input PDF Set to a point in the output (x,Q) grid
            vector<NNPDF::real> EvolPDF(const PDFSet& ipdf,
                                        const size_t iMember,
                                        const size_t ix_out,
                                        const size_t iQ_out) const;

            // Note PIDs must return 21 for gluon, not 0. And 22 for photon if present.
            vector<int>    GetPIDs()          const {return fPIDs;}; // Get a list of returned PIDs
            vector<double> GetEvolvedQ2grid() const {return fXout;}; // Get a list of evolved Q2 points
            vector<double> GetEvolvedXgrid()  const {return fXout;}; // Get a list of evolved x points

        private:
            const vector<int>    fPIDs;  //!< Evolved particle IDs
            const vector<double> fQ2;    //!< Evolved scales
            const vector<double> fXout;  //!< Evolved scale x-grid

            //const NNPDF::FKTable fEvolutionFactors;
    };
}
