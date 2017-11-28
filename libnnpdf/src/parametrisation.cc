// $Id: parametrisation.cc 3195 2015-08-27 10:29:43Z stefano.carrazza@mi.infn.it $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#include <cstdlib>
#include <cmath>
#include <iostream>

#include <fstream>

#include "NNPDF/parametrisation.h"
#include "NNPDF/randomgenerator.h"
#include "NNPDF/exceptions.h"

using namespace NNPDF;

// ******************** Base ******************************
/**
 * @brief Parametrisation base class constructor
 * This is an abstract class handling parameter arrays and I/O
 * @param name a string specifying the type of parametrisation 
 * @param nParameters the number of parameters requested by the derived class 
 */
Parametrisation::Parametrisation( std::string name, int nParameters ):
fNParameters(nParameters),
fParameters(new real[fNParameters]()),
fOutputNorm(1),
fParamName(name)
{
}

/**
 * @brief Parametrisation base class copy constructor
 * @param o the parametrisation to be copied 
 *
 * This constructor copies the model parameters from another parametrisation.
 */
Parametrisation::Parametrisation( Parametrisation const& o):
fNParameters(o.fNParameters),
fParameters(new real[fNParameters]()),
fOutputNorm(o.fOutputNorm),
fParamName(o.fParamName)
{
  CopyPars(&o);
  return;
}

/**
 * @brief Parametrisation base class destructor
 *
 * Frees parametrisation array 
 */
Parametrisation::~Parametrisation()
{
    delete[] fParameters;
}

/**
 * @brief Parametrisation::CopyPars
 * @param t
 */
void Parametrisation::CopyPars(Parametrisation const* t)
{
  if (fNParameters != t->GetNParameters())
    throw EvaluationError("Parametrisation::CopyPars", "number of parameters does not match: " + std::to_string(fNParameters) + " vs " + std::to_string(t->GetNParameters()));
  
  for (int i=0; i<fNParameters; i++)
    fParameters[i] = t->fParameters[i];

  fOutputNorm = t->fOutputNorm;
  
  return;
}

/**
 * @brief Parametrisation::ReadPars
 * @param instream - input for parameters
 */
void Parametrisation::ReadPars(std::ifstream& instream)
{
  instream >> fOutputNorm;
  for (int i=0; i<fNParameters; i++)
    instream >> fParameters[i];
}

/**
 * @brief Parametrisation::WritePars
 * @param outstream - output for parameters
 */
void Parametrisation::WritePars(std::ofstream& outstream)
{
  outstream << " "<<fOutputNorm;
  for (int i=0; i<fNParameters; i++)
    outstream <<" "<< fParameters[i];
}

/**
 * @brief Parametrisation::SetPars
 * @param param - input for parameters
 */
void Parametrisation::SetPars(std::vector<real> const& param)
{
  for (int i=0; i<fNParameters; i++)
    fParameters[i] = param[i];
}

// ******************** MLP *********************************

// Helper function to count the number of neural network parameters for an architecture
int perceptron_count_parameters(std::vector<int> const& arch)
{
    int nparam = 0;
    for (size_t i=1; i<arch.size(); i++)
        nparam += arch[i]*(arch[i-1] + 1);
    return nparam;
}

// The sigmoid actiation function
double sigmoid(double const& h)
{
  return 1.0f/(1.0f+exp(h));
}

/**
 * @brief MultiLayerPerceptron::MultiLayerPerceptron
 * @param settings
 * @param rg
 */
MultiLayerPerceptron::MultiLayerPerceptron(std::vector<int> const& arch):
Parametrisation(std::string("MultiLayerPerceptron"), perceptron_count_parameters(arch)),
fNLayers(arch.size()),
fArch(0),
fWeightMatrix(0),
fOutputMatrix(0),
fActFunction(sigmoid)
{
  int  pSize = 0; // size of parameter array
  int* pMap = new int[fNLayers-1]; // map for weight matrix
  
  // Architecture information
  fArch = new int[fNLayers];
  for (int i=0; i<fNLayers; i++)
  {
    fArch[i] = arch[i];
    
    if (i > 0)
    {
      pMap[i-1] = pSize; // startpoint of this layer
      pSize+=fArch[i]*(1+fArch[i-1]); // size of this layer
    }
  }

  // Alloc parameter array
  if (fNParameters != pSize)
    throw EvaluationError("MultiLayerPerceptron::Constructor", "Error in the computation of nParameters");
  
  // Alloc WeightMatrix (map for parameters)
  fWeightMatrix = new real*[fNLayers - 1];

  // Alloc activation function (output) matrix
  fOutputMatrix = new real*[fNLayers];
  fOutputMatrix[0] = new real[fArch[0]+1];

  for (int i=1; i<fNLayers; i++)
  {
    // point to correct part of parameter array
    fWeightMatrix[i-1] = &fParameters[pMap[i-1]];
    // Alloc each activation function in the layer
    fOutputMatrix[i]   = new real[fArch[i]+1];
  }

  //Init
  fOutputMatrix[0][fArch[0]] =-1.0f;
  for (int i=1; i<fNLayers; i++)
  {
    for (int j=0; j<fArch[i]; j++)
      fOutputMatrix[i][j] = 1;
    
    // Threshold term
    fOutputMatrix[i][fArch[i]] = -1.0f;
  }

  InitParameters();
  
  delete[] pMap;
}

/**
 * @brief MultiLayerPerceptron::MultiLayerPerceptron
 * @param o
 */
MultiLayerPerceptron::MultiLayerPerceptron(MultiLayerPerceptron const& o):
Parametrisation(o),
fNLayers(o.fNLayers),
fArch(0),
fWeightMatrix(0),
fOutputMatrix(0),
fActFunction(o.fActFunction)
{
  int  pSize = 0; // size of parameter array
  int* pMap = new int[fNLayers-1]; // map for weight matrix
  
  // Architecture information
  fArch = new int[fNLayers];
  for (int i=0; i<fNLayers; i++)
  {
    fArch[i] = o.fArch[i];
    
    if (i > 0)
    {
      pMap[i-1] = pSize; // startpoint of this layer
      pSize+=fArch[i]*(1+fArch[i-1]); // size of this layer
    }
  }
  // Alloc WeightMatrix (map for parameters)
  fWeightMatrix = new real*[fNLayers - 1];
  
  // Alloc activation function (output) matrix
  fOutputMatrix = new real*[fNLayers];
  fOutputMatrix[0] = new real[fArch[0]+1];
  
  for (int i=1; i<fNLayers; i++)
  {
    // point to correct part of parameter array
    fWeightMatrix[i-1] = &fParameters[pMap[i-1]];
    // Alloc each activation function in the layer
    fOutputMatrix[i]   = new real[fArch[i]+1];
  }
  
  for (int i=0; i<fNLayers; i++) // Threshold term init
    fOutputMatrix[i][fArch[i]] = -1.0f;
  
  delete[] pMap;
}

/**
 * @brief MultiLayerPerceptron::~MultiLayerPerceptron
 */
MultiLayerPerceptron::~MultiLayerPerceptron()
{
  delete[] fOutputMatrix[0];
  for (int i=1; i<fNLayers; i++)
    delete[] fOutputMatrix[i];
  
  delete[] fOutputMatrix;
  
  delete[] fArch;
  delete[] fWeightMatrix;
}

/**
 * @brief Returns a Parametrisation with identical architecture and type
 */
Parametrisation* MultiLayerPerceptron::Duplicate()
{
  return new MultiLayerPerceptron(*this);
}

/**
 * @brief Initialisation of neural network parameters according to a normal distribution
 */
void MultiLayerPerceptron::InitParameters()
{
  for(int i=0; i<fNParameters; i++)
    fParameters[i] = RandomGenerator::GetRNG()->GetRandomGausDev(1.0);
}

/**
 * @brief MultiLayerPerceptron::Compute
 * @param in
 * @param out
 */
void MultiLayerPerceptron::Compute(real* in,real* out) const
{
  // setup input
  for (int i=0; i<fArch[0]; i++)
    fOutputMatrix[0][i] = in[i];
  
  for (int i=1; i<(fNLayers -1); i++)
    for (int j=0; j<fArch[i]; j++)
    {
      real h=0.0f;
      
      real *p = &fWeightMatrix[i-1][j*(1+fArch[i-1])]; // seems to help the compiler out
      for (int k=0; k<=fArch[i-1]; k++) // <= due to threshold term
        h-= (*(p+k))*fOutputMatrix[i-1][k];
      
      fOutputMatrix[i][j] = fActFunction(h);
    }
  
  // Linear in final layer - get 10% speedup by unrolling this from previous loop
  for (int j=0; j<fArch[fNLayers-1]; j++)
  {
    fOutputMatrix[fNLayers-1][j] = 0.0f;
    
    real *p = &fWeightMatrix[fNLayers-2][j*(1+fArch[fNLayers-2])];
    for (int k=0; k<=fArch[fNLayers-2]; k++)
      fOutputMatrix[fNLayers-1][j] += (*(p+k))*fOutputMatrix[fNLayers-2][k];
  }
  
  for (int i=0; i<fArch[fNLayers-1]; i++)
    out[i] = fOutputNorm*fOutputMatrix[fNLayers-1][i];
}

int MultiLayerPerceptron::GetNumNodeParams(int const& layer) const
{
  if (layer <=0 || layer >= fNLayers)
    throw RangeError("MultiLayerPerceptron::GetNumNodeParams", "layer requested (" + std::to_string(layer) + ") is out of bounds!");

  return fArch[layer-1] + 1;
}

real* MultiLayerPerceptron::GetNodeParams(int const& layer, int const& node)
{
  if (layer <=0 || layer >= fNLayers)
    throw RangeError("MultiLayerPerceptron::GetNodeParams","layer requested (" + std::to_string(layer) + ") is out of bounds!");
  
  if (node < 0 || node >= fArch[layer] )
    throw RangeError("MultiLayerPerceptron::GetNodeParams","node requested (" + std::to_string(node) + ") is out of bounds!");
  
  return &fWeightMatrix[layer-1][node];
}

// ******************** Single (hidden) layer network **************************

/**
 * @brief SingleLayerPerceptron::SingleLayerPerceptron
 * @details The constructor for a single (hidden) layer perceptron class
 * @param arch a vector specifying the NN architecture (must be of the form 2-N-1)
 * @param extra_pars number of additional parameters required (for derived classes)
 */
SingleLayerPerceptron::SingleLayerPerceptron(std::vector<int> const& arch, unsigned int extra_pars):
Parametrisation(std::string("SingleLayerPerceptron"), perceptron_count_parameters(arch) + extra_pars),
fNHidden(arch[1])
{
  if (arch.size() != 3)
    throw EvaluationError("SingleLayerPerceptron::Constructor", "Requested architecture invalid, must be of the form (2,N,1)" );
  if (arch[0] != 2)
    throw EvaluationError("SingleLayerPerceptron::Constructor", "Requested architecture invalid, must be of the form (2,N,1)" );
  if (arch[2] != 1)
    throw EvaluationError("SingleLayerPerceptron::Constructor", "Requested architecture invalid, must be of the form (2,N,1)" );
     
  InitParameters();
}

/**
 * @brief SingleLayerPerceptron::InitParameters
 * @details Initialises the parameters of the neural network, normally distributed.
 */
void SingleLayerPerceptron::InitParameters()
{
  for(int i=0; i<fNParameters; i++)
    fParameters[i] = RandomGenerator::GetRNG()->GetRandomGausDev(1.0);
}

/**
 * @copydoc MultiLayerPerceptron::Duplicate()
 */
Parametrisation* SingleLayerPerceptron::Duplicate()
{
  const std::vector<int> arch = {2, fNHidden, 1};
  return new SingleLayerPerceptron(arch);
}

/**
 * @brief Single layer neural network evaluation 
 * @details Computes the output of the neural network for a given input 
 * @param in, a 2-element array containing {x,log(x)}
 * @param out, a pointer to the location where the output should be written
 */
void SingleLayerPerceptron::Compute(real* in,real* out) const
{
    out[0] = 0;
    // Note no explicit use of fNParameters here, as it includes 'extra_pars'
    const real* final_node = fParameters + 3*fNHidden;
    for (int i=0; i<fNHidden; i++)
        {
                // Use SSE
                const real* node = fParameters + 3*i;
                const real h = in[0]*node[0] 
                             + in[1]*node[1]
                             +       node[2];
                const real sig = 1.0 / (1+exp(h));
                out[0] += final_node[i]*sig;
        }
    out[0] += final_node[fNHidden];
}

// ******************** Preprocessed SLP **************************

/**
 * @copydoc MultiLayerPerceptron::Duplicate()
 */
Parametrisation* SingleLayerPerceptronPreproc::Duplicate()
{
  const std::vector<int> arch = {2, fNHidden, 1};
  return new SingleLayerPerceptronPreproc(arch);
}

/**
 * @brief Initialise the parameters of the preprocessed single layer perceptron 
 * The final two parameters are used for preprocessing exponents
 * small-x: fNParameters - 1
 * large-x: fNParameters - 2
 */
void SingleLayerPerceptronPreproc::InitParameters()
{
    SingleLayerPerceptron::InitParameters();
    fParameters[fNParameters - 2] = RandomGenerator::GetRNG()->GetRandomUniform(0, 2.0)/fScaleFac;
    fParameters[fNParameters - 1] = RandomGenerator::GetRNG()->GetRandomUniform(0, 5.0)/fScaleFac;
}

/**
 * @brief Initialise the parameters of the preprocessed single layer perceptron 
 * Preprocessing is performed on the output of SingleLayerPerceptron::Compute
 */
void SingleLayerPerceptronPreproc::Compute(real* in,real* out) const
{
    // Implemented here as alpha, beta > 0. To set lower limits use normal preproc ranges.
    // Use last two parameters as preprocessing exponents
    SingleLayerPerceptron::Compute(in,out);
    out[0] *= pow(in[0],     fScaleFac*fabs(fParameters[fNParameters -2]));
    out[0] *= pow(1.0-in[0], fScaleFac*fabs(fParameters[fNParameters -1]));
}
