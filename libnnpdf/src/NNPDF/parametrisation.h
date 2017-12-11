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
#include <functional>
using actfunction = std::function<double(double const&)>;

#include "common.h"

namespace NNPDF
{
  /**
   *  \class Parametrisation
   *  \brief Abstract parametrisation base class
   */
  class Parametrisation
  {
  public:
    Parametrisation(std::string name, int nParameters);
    Parametrisation(Parametrisation const& );
    virtual ~Parametrisation();
    
    virtual void Compute(real*,real*) const = 0;
    virtual Parametrisation* Duplicate() = 0;
    void CopyPars(Parametrisation const* );

    void ReadPars(std::ifstream&);
    void WritePars(std::ofstream&);
    void SetPars(std::vector<real> const& param);
    
    real* GetParameters() {return fParameters;}
    int const& GetNParameters() const {return fNParameters;}
    
    std::string const& GetParamName() const {return fParamName;};
    
    virtual void InitParameters() = 0;
    void SetOutputNorm( real const& onorm ) {fOutputNorm = onorm;}; 
    
  protected:
    const int   fNParameters;   //!< Total number of parameters
    real* const fParameters;    //!< Parameters

    real fOutputNorm;         //!< Output normalisation
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
    MultiLayerPerceptron(std::vector<int> const& arch); //!< Network constructor
    MultiLayerPerceptron(MultiLayerPerceptron const&);  //!< Network copy constructor
    ~MultiLayerPerceptron();                            //!< Network destructor
    
    void InitParameters();            //!< Initialize (or reinitialize) parameters
    void Compute(real*,real*) const;  //!< Returns a fArch[fNLayers-1] long array of output for a given input array
    Parametrisation* Duplicate();     //!< Returns a parametrisation based on an MLP

    void SetActivationFunction(actfunction const& f) { fActFunction = f; } //!< set the activation function
   
    // Nodal GA
    const int*  GetArch() const {return fArch;}
    //!< Returns the number of parameters that belong to a specific node (including biases).
    int   GetNumNodeParams(int const& layer) const;
    //!< Returns a pointer to the fParameters coordinate representing the parameters for a specific node
    real*       GetNodeParams   (int const& layer, int const& node);  
    
  protected:
    const int fNLayers;   //!< Number of layers
    int* fArch;           //!< Network architecture
    
    real** fWeightMatrix;         //!< Weights/Thresholds
    mutable real** fOutputMatrix; //!< Neuron Activation

    actfunction fActFunction; //!< Activation function
  };

  /*!
   * \class SingleLayerPerceptron
   * \brief An implementation of a perceptron with a single hidden layer
   * This class can take an arbitrary number of extra parameters, for use
   * in derived classes which are, for example, fitting preprocessing exponents.
   */
  class SingleLayerPerceptron : public Parametrisation
  {
  public:
    SingleLayerPerceptron(std::vector<int> const& arch, unsigned int extra_pars = 0);
    virtual Parametrisation *Duplicate();
    virtual void Compute(real*,real*) const;  
    virtual void InitParameters();  
  protected:
    // Number of hidden nodes
    const int fNHidden; 
  };

  /*!
   * \class SingleLayerPerceptronPreproc
   * \brief A Single-hidden-layer perceptron with preprocessing 
   * Preprocessing is performed to the output of the network as
   *        output*x^{|alpha|}(1-x)^{|beta|}
   * where alpha and beta are the last and second-to-last parameters
   * respectively.
   */
  class SingleLayerPerceptronPreproc : public SingleLayerPerceptron
  {
  public:
    SingleLayerPerceptronPreproc(std::vector<int> const& arch):
    SingleLayerPerceptron(arch, 2){};    
    virtual Parametrisation *Duplicate();
    virtual void InitParameters();  
    void Compute(real*,real*) const;  
   protected:
    // This attribute scales the parameters corresponding to preprocessing
    // exponents in order to balance the training between them and NN params.
    const real fScaleFac = 0.2;
  };

}
