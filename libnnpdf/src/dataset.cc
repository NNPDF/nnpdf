// $Id: dataset.cc 2809 2015-04-24 19:05:41Z stefano.carrazza@mi.infn.it $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <sstream>
#include <utility>
#include <cmath>
#include <sys/stat.h>
#include <iomanip>
#include <cstdlib>

#include "NNPDF/dataset.h"
#include "NNPDF/utils.h"
#include "NNPDF/fastkernel.h"
#include "NNPDF/thpredictions.h"
#include "NNPDF/randomgenerator.h"
#include "NNPDF/exceptions.h"

using namespace NNPDF;

/**
 * Initialization of variables.
 * \param setName std::string with the experimental set name
 * \param ptord perturbative order
 */
DataSet::DataSet(CommonData const& data, FKSet const& set):
  CommonData(data),
  FKSet(set),
  fIsArtificial(false),
  fIsT0(false),
  fT0Pred(new double[fNData]),
  fCovMat(new double*[fNData]),
  fInvCovMat(new double*[fNData])
{
  // Init covariance matrix and T0 Vector
  for (int i = 0; i < fNData; i++)
  {
    fT0Pred[i] = fData[i]; // Default T0 to data

    fCovMat[i] = new double[fNData];
    fInvCovMat[i] = new double[fNData];
    for (int j = 0; j < fNData; j++)
      fCovMat[i][j] = 0.0;
  }
  
  // Generate Covariance Matrix
  GenCovMat();

}

DataSet::DataSet(const DataSet& set):
  CommonData(set),
  FKSet(set),
  fIsArtificial(set.fIsArtificial),
  fIsT0(set.fIsT0),
  fT0Pred(new double[fNData]),
  fCovMat(new double*[fNData]),
  fInvCovMat(new double*[fNData])
{
  // Building Obs array
  for (int i = 0; i < fNData; i++)
    {
      fT0Pred[i] = set.fT0Pred[i];
      
      fCovMat[i] = new double[fNData];
      fInvCovMat[i] = new double[fNData];
      
      for (int j=0; j<fNData; j++)
      {
        fCovMat[i][j] = 0.0;
        fInvCovMat[i][j] = 0.0;
      }
    }

  // Generate covariance matrix
  GenCovMat();
}

DataSet::DataSet(const DataSet& set, std::vector<int> const& mask):
  CommonData(set, mask),
  FKSet(set, mask),
  fIsArtificial(set.fIsArtificial),
  fIsT0(set.fIsT0),
  fT0Pred(new double[fNData]),
  fCovMat(new double*[fNData]),
  fInvCovMat(new double*[fNData])
{
  // Building Obs array
  for (int i = 0; i < fNData; i++)
    {
      fT0Pred[i] = set.fT0Pred[mask[i]];
      
      fCovMat[i] = new double[fNData];
      fInvCovMat[i] = new double[fNData];
      
      for (int j=0; j<fNData; j++)
      {
        fCovMat[i][j] = 0.0;
        fInvCovMat[i][j] = 0.0;
      }
    }

  // Generate covariance matrix
  GenCovMat();
}

/**
 * Destroys variables
 */
DataSet::~DataSet()
{
  // Delete T0
  delete[] fT0Pred;

  if (fCovMat) {
    for (int i = 0; i < fNData; i++) {
      if (fCovMat[i])
        delete[] fCovMat[i];
    }
    delete[] fCovMat;
  }

  if (fInvCovMat) {
    for (int i = 0; i < fNData; i++) {
      if (fInvCovMat[i])
        delete[] fInvCovMat[i];
    }
    delete[] fInvCovMat;
  }
  
}

/**
 * Generate covariance matrix and inverse
 */
void DataSet::GenCovMat()
{  
  for (int i = 0; i < fNData; i++)
    for (int j = 0; j < fNData; j++)
    {
      double sig    = 0.0;
      double signor = 0.0;
      
      if (i == j)
        sig += fStat[i]*fStat[i]; // stat error
      
      for (int l = 0; l < fNSys; l++)
        if (i == j || fSys[i][l].name.compare("UNCORR")!=0)
          switch (fSys[i][l].type)
          {
            case ADD: sig += fSys[i][l].add*fSys[j][l].add; break; // additive systematics
            case MULT: signor += fSys[i][l].mult*fSys[j][l].mult; break; // multiplicative systematics
          }
              
      fCovMat[i][j] = sig + signor*fT0Pred[i]*fT0Pred[j]*1e-4;
    }
  
  InvertLU(fNData, fCovMat, fInvCovMat);
}

void DataSet::RescaleErrors(const double mult)
{
  for (int i = 0; i < fNData; i++)
  {
    fStat[i] *= mult;         // Rescale stat and sys to new data value
    
    for (int l = 0; l < fNSys; l++)
    {  
      fSys[i][l].add *=  mult;
      fSys[i][l].mult *= mult;
    }
  }
}

void DataSet::SetT0(ThPredictions const& t0pred)
{
  // switch flag
  fIsT0 = true;

  if (fNData != t0pred.GetNData())
    throw RangeError("DataSet::SetT0","number of datapoints in set and predictions do not match!");

  for (int i=0; i<fNData; i++)
    fT0Pred[i] = t0pred.GetObsCV(i);

  // Regenerate covariance matrix
  GenCovMat();
}


/**
 * @brief Dataset::MakeArtificial
 * Produces the shifts by using a MC sampling
 */
void DataSet::MakeArtificial()
{
  if (DataSet::Verbose)
    std::cout << "-- Generating replica data for " << fSetName << std::endl;

  double *rand  = new double[fNSys];  
  double *xnor  = new double[fNData];
  double *artdata = new double[fNData];
  sysType rST[2] = {ADD,MULT};
  
  NNPDF::RandomGenerator* rng = NNPDF::RandomGenerator::GetRNG();

  // Generate random systematics
  for (int l = 0; l < fNSys; l++)
  {
    rand[l] = rng->GetRandomGausDev(1.0);
    if (fSys[0][l].isRAND)
      fSys[0][l].type = rST[rng->GetRandomUniform(2)];
    for (int i = 1; i < fNData; i++)
      fSys[i][l].type = fSys[0][l].type;
  } 

  int genTries = 0; // Number of attempts at generation
  bool genRep = false;
  while ( genTries < 1E6 ) 
  {
    
    // Update the data
    for (int i = 0; i < fNData; i++)
    {
      const double xstat = rng->GetRandomGausDev(1.0)*fStat[i];   //Noise from statistical uncertainty
      
      double xadd = 0;
      xnor[i] = 1.0;
      
      for (int l = 0; l < fNSys; l++)
      {
        if (fSys[i][l].name.compare("UNCORR")==0)         //Noise from uncorrelated systematics
        {
          switch (fSys[i][l].type)
          {
            case ADD: xadd += rng->GetRandomGausDev(1.0)*fSys[i][l].add; break;
            case MULT: xnor[i] *= (1.0 + rng->GetRandomGausDev(1.0)*fSys[i][l].mult*1e-2); break;
          }
        }
        else                                              //Noise from correlated systematics
        {
          switch (fSys[i][l].type)
          {
            case ADD: xadd += rand[l]*fSys[i][l].add; break;
            case MULT: xnor[i] *= (1.0 + rand[l]*fSys[i][l].mult*1e-2); break;
          }
        }
      }

      // Generate the artificial data
      artdata[i] = xnor[i] * ( fData[i] + xadd + xstat );
      
      // Only generates positive artificial data (except for closure tests)
      if ( artdata[i] < 0 )
        break;

    }
    
    // Stop when all the datapoints are positive
    genRep = true;
    break;
  }

  if (!genRep)
    throw EvaluationError("DataSet::MakeArtificial", "Cannot generate positive datasets!");

  // Update data in set
  UpdateData(artdata);
  fIsArtificial = true;
  
  delete[] rand;
  delete[] xnor;
  delete[] artdata;

  // Note DO NOT Regenerate covariance matrices  
}



/**
 * Update data values - used by MakeArtificial and MakeReplica
 */
void DataSet::UpdateData(double *newdata)      // Update data only
{
  for (int i = 0; i < fNData; i++)
    fData[i] = newdata[i];                   // Set new data values
}

void DataSet::UpdateData(double *newdata, double* uncnorm)      // Update data only and rescale systematics by normalisation
{
  UpdateData(newdata);
  
  for (int i = 0; i < fNData; i++)
  {
    fStat[i] *= uncnorm[i];         // Rescale stat and sys to new data value
    
    for (int l = 0; l < fNSys; l++)
      fSys[i][l].add *=  uncnorm[i];
  }
}

void DataSet::UpdateData(double *newdata, sysType *type)       // Update data and change systypes
{
  UpdateData(newdata);
  
  for (int l = 0; l < fNSys; l++)
    for (int i = 0; i < fNData; i++)
      fSys[i][l].type = type[l]; 
}
