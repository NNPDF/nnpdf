// $Id: experiments.h 2990 2015-06-08 19:06:21Z s1044006 $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include <vector>
#include <string>

#include "commondata.h"
#include "dataset.h"
#include "pdfset.h"

namespace NNPDF
{
  /**
    * \class Experiment
    * \brief Class for handling experiments
    */
  class Experiment
  {
  public:
    Experiment(std::vector<DataSet> const& sets, std::string const& expname); //!< Experiment constructor
    Experiment(Experiment const&, std::vector<DataSet> const& sets);   //!< Copy constructor with dataset replace
    Experiment(Experiment const&);   //!< Copy constructor
    ~Experiment();                  //!< Experiment destructor
    
    void MakeReplica();                                       //!< Shift the exp data and produces art. data
    void MakeClosure(const std::vector<ThPredictions>& predictions, bool const& noise);
    void MakeClosure(PDFSet* pdf, bool const& noise);         //!< Make fake exp data using theory predictions

    int GetNSet() const {return (int) fSets.size();}    //!< Return the number of sets in the experiment
    DataSet const&  GetSet(int i)   const {return fSets[i];}        //!< Return the ith DataSet  
    std::string const& GetExpName() const {return fExpName;}        //!< Return the experiment name
    std::string const& GetSetName(int i) const;
    
    int  GetNData() const { return fNData; }                  //!< Return the number of points in the experiment
    const double* GetData() const { return fData; }
    
    bool IsArtificial() const { return fIsArtificial; }       //!< Return the artificial flag
    bool IsClosure() const { return fIsClosure; }             //!< Return the artificial flag
    bool IsT0() const { return fIsT0; }                       //!< Return t0 covmat flag

    double** GetCovMat()         const { return fCovMat;    } //!< Return fCovMat
    double** GetSqrtCov()        const { return fSqrtCov;   } //!< Return the Cholesky decomposition of the covariance matrix
    
    void ExportCovMat(std::string);         //!< Export covariance matrix
    void ExportSqrtCov(std::string);        //!< Export Cholesky decomposition

    void SetT0(const PDFSet&); //!<Set T0 Predictions for each dataset in place and update internal structures

  private:

    Experiment();                               //disable default constructor
    Experiment& operator=(const Experiment&) ;  //disable copy-assignment

    std::string fExpName;                                    //!< Stores the name of the experiment

    std::vector<DataSet> fSets;                             //!< Vector of contained datasets

    int fNData;        //!< Number of data points
    int fNSys;      //!< Number of additive systematics correlations
    
    double *fData;       //!< The experimental data
    double *fT0Pred;     //!< The t0 predictions

    double **fCovMat;    //!< The covariance matrix
    double **fSqrtCov;   //!< The Cholesky decomposition of the covariance matrix
    
    double *fStat;       //!< The statistical errors
    sysError **fSys;    //!< The syscor
    
    int **fSetSysMap;    //!< Map for ordering of systematics in datasets
    
    bool fIsArtificial; //!< The artificial flag
    bool fIsClosure;    //!< If data is closure data 
    bool fIsT0;         //!< Flag to determine if covmat is T0

        
    void ClearLocalData();            //!< Clear data pulled from datasets
    void ClearLocalCovMat();          //!< Clear covariance matrices

    void PullData();                  //!< Pull experimental data from datasets
    void GenCovMat();                 //!< Generate covmat and inverse

  };  

}
