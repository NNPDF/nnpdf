// $Id: parametrisation.h 3196 2015-08-27 14:16:27Z stefano.carrazza@mi.infn.it $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include <cstdio>
#include <cstdlib>
#include <vector>

#include "common.h"

namespace NNPDF
{

  /**
   *  \class Parametrisation
   *  \brief Virtual parametrisaton base class
   */
  class Parametrisation
  {
  public:
    Parametrisation(std::string name);
    Parametrisation(Parametrisation const& );
    virtual ~Parametrisation();
    
    virtual void Compute(real*,real*) const = 0;

    void CopyPars(Parametrisation*);
    virtual Parametrisation* Duplicate() = 0;

    void ReadPars(std::ifstream&);
    void WritePars(std::ofstream&);
    
    real* GetParameters() {return fParameters;}
    int const& GetNParameters() const {return fNParameters;}
    
    std::string const& GetParamName() const {return fParamName;};
    
    virtual void InitParameters() = 0;
    virtual int GetNumNodeParams(int const& layer) const = 0; // sc - to improve
    void SetOutputNorm( real const& onorm ) {fOutputNorm = onorm;}; 
    
  protected:
    int   fNParameters;   //!< Total number of parameters
    real* fParameters;    //!< Parameters

    real fOutputNorm; //!< Output normalisation
    
  private:
    const std::string fParamName;
  };

  /**
   *  \class MultiLayerPerceptron
   *  \brief General MLP class
   */

  class MultiLayerPerceptron : public Parametrisation
  {
  public:
    MultiLayerPerceptron(std::vector<int> const& arch);         //!< Network constructor
    MultiLayerPerceptron(MultiLayerPerceptron const&);  //!< Network copy constructor
    ~MultiLayerPerceptron();                            //!< Network destructor
    
    void Compute(real*,real*) const;  //!< Returns a fArch[fNLayers-1] long array of output for a given input array
    Parametrisation* Duplicate();     //!< Returns a parametrisation based on an MLP
    
    const int*  GetArch() const {return fArch;};
    int   GetNumNodeParams(int const& layer) const;                   //!< Returns the number of parameters that belong to a specific node (including biases).
    real*       GetNodeParams   (int const& layer, int const& node);  //!< Returns a pointer to the fParameters coordinate representing the parameters for a specific node
    
    void InitParameters();   //!< Initialize (or reinitialize) parameters
    
  private:
    const int fNLayers;   //!< Number of layers
    int* fArch;           //!< Network architecture
    
    real** fWeightMatrix;         //!< Weights/Thresholds
    mutable real** fOutputMatrix; //!< Neuron Activation
  };

  /**
   *  \class ChebyshevPolynomial
   *  \brief Chebyshev Polynomial Parametrisation
   */

  class ChebyshevPolynomial : public Parametrisation
  {
  public:
    ChebyshevPolynomial(std::vector<int> const& order); //!< Parametrisation constructor - improve way to pass orders
    ChebyshevPolynomial(ChebyshevPolynomial const&);
    ~ChebyshevPolynomial();                    //!< Parametrisation destructor
    
    void Compute(real* ,real* ) const;
    
    Parametrisation* Duplicate();

    int GetNumNodeParams(int const& layer) const { return 0; }
    void InitParameters();   //!< Initialize (or reinitialize) parameters

  private:
    const int fNOrder;    //!< Highest order polynomial
    
    mutable real* fPolynomials;   //!< Used in the Compute method
  };

  /**
   *  \class MultiLayerPerceptron
   *  \brief General MLP class
   */

  class QuadMultiLayerPerceptron : public Parametrisation
  {
  public:
    QuadMultiLayerPerceptron(std::vector<int> const& arch);         //!< Network constructor
    QuadMultiLayerPerceptron(QuadMultiLayerPerceptron const&);  //!< Network copy constructor
    ~QuadMultiLayerPerceptron();                            //!< Network destructor

    void Compute(real*,real*) const;  //!< Returns a fArch[fNLayers-1] long array of output for a given input array
    Parametrisation* Duplicate();     //!< Returns a parametrisation based on an MLP

    const int*  GetArch() const {return fArch;}
    int   GetNumNodeParams(int const& layer) const;                   //!< Returns the number of parameters that belong to a specific node (including biases).

    void InitParameters();   //!< Initialize (or reinitialize) parameters

  private:
    const int fNLayers;   //!< Number of layers
    int* fArch;           //!< Network architecture

    real** fWeightMatrix;         //!< Weights/Thresholds
    mutable real** fOutputMatrix; //!< Neuron Activation
  };

}
