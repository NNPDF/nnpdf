// $Id: positivity.h 3187 2015-08-23 11:08:42Z stefano.carrazza@mi.infn.it $
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk
#pragma once

#include <string>
#include <vector>

#include "common.h"
#include "fastkernel.h"
#include "commondata.h"

namespace NNPDF
{

  class PDFSet;

  /**
   * \class PositivitySet Constraint Class
   * \brief Positivity observables class manager
   */

  class PositivitySet : public CommonData, public FKTable
  {
  public:
    PositivitySet(CommonData const&, FKTable const&, real const& lambda); //!< Constructor
    PositivitySet(const PositivitySet&);  //!< Copy constructor
    virtual ~PositivitySet(); //!< Destructor.    

    void ComputeErf   (const PDFSet*, real*) const; //!< Compute error function for supplied PDF set
    
    void ComputeNViolated( const PDFSet*, int* res) const; //!< Number of violated points
    void ComputeNUnacceptable(const PDFSet*, int* res) const; //!< Number of unacceptable points (less than bounds)

    void SetBounds(const PDFSet*); //!< Set bounds
     
  private:

    PositivitySet(); // disable default constructor
    PositivitySet& operator=(const PositivitySet&); // disable copy-assignment

    void ComputePoints(const PDFSet*, real*) const; //!< Compute the positivity points

    const real        fLambda; //!< Lagrange Multiplier
    
    std::vector<real> fBounds;
  };

}
