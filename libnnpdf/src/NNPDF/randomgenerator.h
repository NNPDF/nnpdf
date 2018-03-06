// $Id: randomgenerator.h 2120 2014-12-02 17:14:52Z s1044006 $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "common.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

#include <vector>

namespace NNPDF
{
  /**
   * @brief The RandomGenerator class
   * This class manages the random number generators
   */
  class RandomGenerator
  {
  private:
    RandomGenerator(int,unsigned long int);
    ~RandomGenerator();

    gsl_rng *fR;

  public:

    // Initialise and fetch RNG
    static void InitRNG(int method,unsigned long int seed);
    static RandomGenerator* GetRNG();

    // Verbosity control
    static bool Verbose;

    unsigned long int GetRandomInt();
    unsigned long int GetRandomUniform(unsigned long int);
    double GetRandomUniform();
    double GetRandomUniform(double, double);
    double GetRandomUniformPos();
    double GetRandomGausDev(const double);

    template<typename T>
    void ShuffleVector( std::vector<T>& vec )
    {return gsl_ran_shuffle (fR, &vec[0], vec.size(), sizeof(T));}

    // Set Methods
    void SetSeed(unsigned long int);
  };
}
