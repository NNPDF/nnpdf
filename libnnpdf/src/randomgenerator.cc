// $Id: randomgenerator.cc 2120 2014-12-02 17:14:52Z s1044006 $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#include <iostream>
#include "NNPDF/randomgenerator.h"
#include "NNPDF/exceptions.h"

namespace NNPDF
{
  // Singleton
  static RandomGenerator* rngInstance;

  // See GSL documentation
  const gsl_rng_type* rngs[] = { gsl_rng_ranlux,
                                 gsl_rng_cmrg,
                                 gsl_rng_mrg,
                                 gsl_rng_mt19937,
                                 gsl_rng_gfsr4,
                                 gsl_rng_ran0,
                                 gsl_rng_ran1,
                                 gsl_rng_ran2,
                                 gsl_rng_ran3,
                                 gsl_rng_rand,
                                 gsl_rng_ranlxd1,
                                 gsl_rng_ranlxd2,
                                 gsl_rng_ranlxs0,
                                 gsl_rng_ranlxs1,
                                 gsl_rng_ranlxs2,
                                 gsl_rng_taus,
                                 gsl_rng_taus2
                               };

  // Default Verbosity
  bool RandomGenerator::Verbose = true;


   // Management methods
  void RandomGenerator::InitRNG(int method,unsigned long int seed)
  {
    if (rngInstance != 0)
    delete rngInstance;

    rngInstance = new RandomGenerator(method, seed);
  }

  RandomGenerator* RandomGenerator::GetRNG()
  {
    if (rngInstance == 0)
      throw InitError("GetRNG()","RNG has not been initialised!");

    return rngInstance;
  }

  /**
   * @brief RandomGenerator::RandomGenerator
   */
  RandomGenerator::RandomGenerator(int method, unsigned long int seed)
  {
    // just export GSL_RNG_TYPE=mrg to change the generator
    //gsl_rng_env_setup();
    const gsl_rng_type *T = rngs[method];

    fR = gsl_rng_alloc(T);
    if (Verbose)
      std::cout << "- Random Generator allocated: " << gsl_rng_name(fR) << std::endl;
    gsl_rng_set(fR, seed);
  }

  /**
   * @brief RandomGenerator::~RandomGenerator
   */
  RandomGenerator::~RandomGenerator()
  {
    gsl_rng_free(fR);
  }

  /**
   * @brief RandomGenerator::SetSeed
   * @param s
   */
  void RandomGenerator::SetSeed(unsigned long int s)
  {
    gsl_rng_set(fR, s);
  }

  /**
   * @brief RandomGenerator::GetRandomInt
   * @return a random integer from the generator.
   * Min max are equally likely and depends on the
   * algorithm used.
   */
  unsigned long int RandomGenerator::GetRandomInt()
  {
    return gsl_rng_get(fR);
  }

  /**
   * @brief RandomGenerator::GetRandomUniform
   * @param n
   * @return an integer from 0 to n-1
   */
  unsigned long int RandomGenerator::GetRandomUniform(unsigned long int n)
  {
    return gsl_rng_uniform_int(fR, n);
  }

  /**
   * @brief RandomGenerator::GetRandomUniform
   * @return a double from [0,1)
   */
  double RandomGenerator::GetRandomUniform()
  {
    return gsl_rng_uniform(fR);
  }

  /**
   * @brief RandomGenerator::GetRandomUniform
   * @return a double from [a,b)
   */
  double RandomGenerator::GetRandomUniform(double a, double b)
  {
    return (b-a)*gsl_rng_uniform(fR) + a;
  }

  /**
   * @brief RandomGenerator::GetRandomUniformPos
   * @return a double from (0,1)
   */
  double RandomGenerator::GetRandomUniformPos()
  {
    return gsl_rng_uniform_pos(fR);
  }

  /**
   * @brief Gaussian Random number generator
   * using Box-Muller algorithm
   * @param x
   * @return
   */
  double RandomGenerator::GetRandomGausDev(const double sigma)
  {
    return gsl_ran_gaussian(fR, sigma);
  }
}
