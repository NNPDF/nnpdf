// $Id: dataset.h 2990 2015-06-08 19:06:21Z s1044006 $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "common.h"
#include <string>
#include <fstream>
#include <algorithm>

#include "fkset.h"
#include "commondata.h"

namespace NNPDF
{

  class ThPredictions;


  /**
   *  \class DataSet
   *  \brief Class for datasets loading and computing of observables
   */
  class DataSet : public CommonData, public FKSet
  {
   private:    
    bool fIsArtificial; //!< Flag to determine if data is artifical
    bool fIsT0;         //!< Flag to determine if covmat is T0

    // Data information
    double *fT0Pred;      //!< The t0 predictions - defaults to data in case of no t0-std::vector
    double **fCovMat;     //!< The covariance matrix
    double **fInvCovMat;  //!< The inverse cov matrix
    
    // private methods for constructor
    void GenCovMat();                       //!< Generate covariance matrix

    DataSet();                          //disable default constructor
    DataSet& operator=(const DataSet&); //disable copy-assignment

   public:
    DataSet(CommonData const&, FKSet const&); //!< Constructor
    virtual ~DataSet();                       //!< The destructor.    

    DataSet(const DataSet&);     //!< Copy constructor
    DataSet(const DataSet&, std::vector<int> const& );     //!< Masked Copy constructor

    // ****************   DataSet T0 Methods  **********************

    void   SetT0( ThPredictions  const&);               //!< Set T0 predictions
    bool   const& IsT0 ()  const { return fIsT0; }      //!< Return t0 covmat flag
    
    // ************************ Data Get Methods ******************************
    
    double const&  GetT0Pred(int i)    const { return fT0Pred[i];} //!< Return t0 prediction
    
    double** GetCovMat()         const { return fCovMat;   } //!< Return fCovMat
    double** GetInvCovMat()      const { return fInvCovMat;} //!< Return the inverse of CovMat
    double const& GetInvCovMat( int const& i, int const& j) const {return fInvCovMat[i][j];}; //!< Returns a iCov element

    bool const& IsArtificial()         const { return fIsArtificial;} //!< Returns the artificial flag
    
    // ****************   Update data values  ********************************

    void   MakeArtificial(); //!< Make an artificial data replica
    
    void   UpdateData(double* newdat);                //!< Update data
    void   UpdateData(double* newdat, double* norm);  //!< Update with a rescaling - also rescales additive uncertainties
    void   UpdateData(double* newdat, sysType* type); //!< Update data and systypes

    void SetArtificial( bool const& artificial) {fIsArtificial = artificial;};
    
    // ****************   Rescale Uncertainties   ****************************
    
    void RescaleErrors(const double mult);
    
    // **************************************************************************
  
  };
  
};
