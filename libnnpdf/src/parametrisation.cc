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
 * @brief Parametrisation constructor
 */
Parametrisation::Parametrisation( std::string name ):
fNParameters(0),
fParameters(0),
fOutputNorm(1),
fParamName(name)
{
}

/**
 * @brief Parametrisation constructor
 * @param o the copy contructor
 */
Parametrisation::Parametrisation( Parametrisation const& o):
fNParameters(o.fNParameters),
fParameters(new real[fNParameters]),
fOutputNorm(o.fOutputNorm),
fParamName(o.fParamName)
{
  for (int i=0; i<fNParameters; i++)
    fParameters[i] = o.fParameters[i];
  
  return;
}

/**
 * @brief Parametrisation destructor
 */
Parametrisation::~Parametrisation()
{
}

/**
 * @brief Parametrisation::CopyPars
 * @param t
 */
void Parametrisation::CopyPars(Parametrisation* t)
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
/**
 * @brief MultiLayerPerceptron::MultiLayerPerceptron
 * @param settings
 * @param rg
 */
MultiLayerPerceptron::MultiLayerPerceptron(std::vector<int> const& arch):
Parametrisation(std::string("MultiLayerPerceptron")),
fNLayers(arch.size()),
fArch(0),
fWeightMatrix(0),
fOutputMatrix(0)
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
  fParameters = new real[pSize];
  fNParameters = pSize;
  
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
fOutputMatrix(0)
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
  delete[] fParameters;
}

/**
 * @brief MultiLayerPerceptron::Duplicate
 * @return
 */
Parametrisation* MultiLayerPerceptron::Duplicate()
{
  return new MultiLayerPerceptron(*this);
}

/**
 * @brief MultiLayerPerceptron::InitParameters
 */
void MultiLayerPerceptron::InitParameters()
{
  //Init
  for (int i=1; i<fNLayers; i++)
  {
    const int wsz = (fArch[i-1]+1)*fArch[i];    
    for (int j=0; j<wsz; j++)
      fWeightMatrix[i-1][j] = RandomGenerator::GetRNG()->GetRandomGausDev(1.0);
  }
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
      
      fOutputMatrix[i][j]=1.0f/(1.0f+exp(h));
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


// ******************** MLP with PREPROC **************************
/**
 * @brief MultiLayerPerceptronPreproc::MultiLayerPerceptronPreproc
 * @param settings
 * @param rg
 */
MultiLayerPerceptronPreproc::MultiLayerPerceptronPreproc(std::vector<int> const& arch):
MultiLayerPerceptron(arch)
{  
  // clenup constructor
  delete[] fOutputMatrix[0];
  for (int i=1; i<fNLayers; i++)
    delete[] fOutputMatrix[i];

  delete[] fOutputMatrix;

  delete[] fArch;
  delete[] fWeightMatrix;
  delete[] fParameters;

  // reallocated objects
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
  fParameters = new real[pSize];
  fNParameters = pSize;

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
}

/**
 * @brief MultiLayerPerceptronPreproc::MultiLayerPerceptronPreproc
 * @param o
 */
MultiLayerPerceptronPreproc::MultiLayerPerceptronPreproc(MultiLayerPerceptronPreproc const& o):
MultiLayerPerceptron(o)
{
}

/**
 * @brief MultiLayerPerceptronPreproc::InitParameters
 */
void MultiLayerPerceptronPreproc::InitParameters()
{
  MultiLayerPerceptron::InitParameters();
  fParameters[fNParameters-2] = RandomGenerator::GetRNG()->GetRandomUniform(-10,2);
  fParameters[fNParameters-1] = RandomGenerator::GetRNG()->GetRandomUniform(-10,10);
}

/**
 * @brief MultiLayerPerceptronPreproc::Compute
 * @param in
 * @param out
 */
void MultiLayerPerceptronPreproc::Compute(real* in,real* out) const
{
  MultiLayerPerceptron::Compute(in,out);

  // apply preprocessing alpha = fNParameters-2 and beta = fNParameters-1.
  for (int i=0; i< fArch[fNLayers-1]; i++)
    out[i] *= pow(in[i],-fParameters[fNParameters-2])*pow(1-in[i],fParameters[fNParameters-1]);
}

// ******************** Chebyshev *********************************
/**
 * @brief ChebyshevPolynomial::ChebyshevPolynomial
 * @param settings
 * @param rg
 */
ChebyshevPolynomial::ChebyshevPolynomial(std::vector<int> const& order):
Parametrisation(std::string("ChebyshevPolynomial")),
fNOrder(10),
fPolynomials(new real[fNOrder])
{
  fNParameters = fNOrder;
  fParameters = new real[fNOrder];
  
  if (fNOrder<2)
    throw InitError("ChebyshevPolynomial()", "Must have at least 2 orders");
  
  // init
  fPolynomials[0] = 1;
  for (int i=1; i<fNOrder; i++)
    fPolynomials[i] = 0.0f;

  // Init parameters
  InitParameters();
}

/**
 * @brief ChebyshevPolynomial::ChebyshevPolynomial
 * @param o
 */
ChebyshevPolynomial::ChebyshevPolynomial(ChebyshevPolynomial const& o):
Parametrisation(o),
fNOrder(o.fNOrder),
fPolynomials(new real[fNOrder])
{
  if (fNOrder<2)
    throw InitError("ChebyshevPolynomial()", "Must have at least 2 orders");
  
  // init
  fPolynomials[0] = 1;
  for (int i=1; i<fNOrder; i++)
    fPolynomials[i] = 0.0f;
}

/**
 * @brief ChebyshevPolynomial::~ChebyshevPolynomial
 */
ChebyshevPolynomial::~ChebyshevPolynomial()
{
  delete[] fParameters;
  delete[] fPolynomials;
}

/**
 * @brief ChebyshevPolynomial::Duplicate
 * @return
 */
Parametrisation* ChebyshevPolynomial::Duplicate()
{
  return new ChebyshevPolynomial(*this);
}

/**
 * @brief MultiLayerPerceptron::InitParameters
 */
void ChebyshevPolynomial::InitParameters()
{
  // init
  for (int i=0; i<fNOrder; i++)
    fParameters[i] = RandomGenerator::GetRNG()->GetRandomGausDev(1.0);

}  
/**
 * @brief ChebyshevPolynomial::Compute
 * @param input
 * @param output
 */
void ChebyshevPolynomial::Compute(real* input,real* output) const // this isnt particuarly fast actually - thanks to arbitary x values
{
  // use transform y = 1-2*sqrt(x) to move to -1 < y < 1
  const real y = 1- 2*sqrt(input[0]);

  fPolynomials[1] = y;
  for (int i=2; i<fNOrder; i++)
    fPolynomials[i] = 2*y*fPolynomials[i-1] - fPolynomials[i-2];
  
  // result
  output[0] = 0.0f;
  for (int i=0; i<fNOrder; i++)
    output[0] += fOutputNorm*fParameters[i]*fPolynomials[i];
}

// ******************** MLP *********************************
/**
 * @brief QuadMultiLayerPerceptron::QuadMultiLayerPerceptron
 * @param settings
 * @param rg
 */
QuadMultiLayerPerceptron::QuadMultiLayerPerceptron(std::vector<int> const& arch):
Parametrisation(std::string("QuadMultiLayerPerceptron")),
fNLayers(arch.size()),
fArch(0),
fWeightMatrix(0),
fOutputMatrix(0)
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
      pSize+=fArch[i]*(1+fArch[i-1]+fArch[i-1]*fArch[i-1]); // size of this layer
    }
  }

  // Alloc parameter array
  fParameters = new real[pSize];
  fNParameters = pSize;

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
QuadMultiLayerPerceptron::QuadMultiLayerPerceptron(QuadMultiLayerPerceptron const& o):
Parametrisation(o),
fNLayers(o.fNLayers),
fArch(0),
fWeightMatrix(0),
fOutputMatrix(0)
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
      pSize+=fArch[i]*(1+fArch[i-1]+fArch[i-1]*fArch[i-1]); // size of this layer
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
QuadMultiLayerPerceptron::~QuadMultiLayerPerceptron()
{
  delete[] fOutputMatrix[0];
  for (int i=1; i<fNLayers; i++)
    delete[] fOutputMatrix[i];

  delete[] fOutputMatrix;

  delete[] fArch;
  delete[] fWeightMatrix;
  delete[] fParameters;
}

/**
 * @brief MultiLayerPerceptron::Duplicate
 * @return
 */
Parametrisation* QuadMultiLayerPerceptron::Duplicate()
{
  return new QuadMultiLayerPerceptron(*this);
}

/**
 * @brief MultiLayerPerceptron::InitParameters
 */
void QuadMultiLayerPerceptron::InitParameters()
{
  //Init
  for (int i=1; i<fNLayers; i++)
  {
    const int wsz = (fArch[i-1]*fArch[i-1]+fArch[i-1]+1)*fArch[i];
    for (int j=0; j<wsz; j++)
      fWeightMatrix[i-1][j] = RandomGenerator::GetRNG()->GetRandomGausDev(1.0);
  }
}

/**
 * @brief MultiLayerPerceptron::Compute
 * @param in
 * @param out
 */
void QuadMultiLayerPerceptron::Compute(real* in,real* out) const
{
  // setup input
  for (int i=0; i<fArch[0]; i++)
    fOutputMatrix[0][i] = in[i];

  for (int i=1; i<(fNLayers -1); i++)
    for (int j=0; j<fArch[i]; j++)
    {
      real h=0.0f;

      // linear term + threshold
      real *p = &fWeightMatrix[i-1][j*(1+fArch[i-1]+fArch[i-1]*fArch[i-1])]; // seems to help the compiler out
      for (int k=0; k<=fArch[i-1]; k++) // <= due to threshold term
        h-= (*(p+k))*fOutputMatrix[i-1][k];

      // quadratic term
      for (int k = 0; k < fArch[i-1]; k++)
        for (int l = 0; l < fArch[i-1]; l++)
          h -= (*(p+fArch[i-1]+1+k*fArch[i-1]+l))*fOutputMatrix[i-1][k]*fOutputMatrix[i-1][l];

      fOutputMatrix[i][j]=1.0f/(1.0f+exp(h));
    }

  // Linear in final layer - get 10% speedup by unrolling this from previous loop
  for (int j=0; j<fArch[fNLayers-1]; j++)
  {
    fOutputMatrix[fNLayers-1][j] = 0.0f;

    // linear term
    real *p = &fWeightMatrix[fNLayers-2][j*(1+fArch[fNLayers-2]+fArch[fNLayers-2]*fArch[fNLayers-2])];
    for (int k=0; k<=fArch[fNLayers-2]; k++)
      fOutputMatrix[fNLayers-1][j] += (*(p+k))*fOutputMatrix[fNLayers-2][k];

    // quadratic term
    for (int k = 0; k < fArch[fNLayers-2]; k++)
      for (int l = 0; l < fArch[fNLayers-2]; l++)
        fOutputMatrix[fNLayers-1][j] += (*(p+fArch[fNLayers-2]+1+k*fArch[fNLayers-2]+l))*fOutputMatrix[fNLayers-2][k]*fOutputMatrix[fNLayers-2][l];
  }

  for (int i=0; i<fArch[fNLayers-1]; i++)
    out[i] = fOutputNorm*fOutputMatrix[fNLayers-1][i];
}

int QuadMultiLayerPerceptron::GetNumNodeParams(int const& layer) const
{
  if (layer <=0 || layer >= fNLayers)
    throw RangeError("QuadMultiLayerPerceptron::GetNumNodeParams","layer requested (" + std::to_string(layer) + ") is out of bounds!");

  return fArch[layer-1]*fArch[layer-1]+fArch[layer-1] + 1;
}
