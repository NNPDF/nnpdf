// $Id: utils.h 2825 2015-05-03 09:14:46Z stefano.carrazza@mi.infn.it $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "common.h"
#include <vector>
#include <string>
#include <fstream>

/** @defgroup utils Utils
 * \brief libnnpdf utility functions
 */

/*! \addtogroup utils
 *  @{
 */

namespace NNPDF
{
  /**
   * @brief untargz Decompress tar.gz files in istream.
   * @param filename the input filename
   * @return the istream object
   */
  std::vector<char> untargz(std::string const& filename);

  /**
   * @brief targz Store to disk the compressed data
   * @param filename the output file
   * @param data the stream buffer
   */
  void targz(std::string const& filename, std::stringstream const& data);

  // *******************  SWIG helpers *****************************

  /*!
   *  \class istream_proxy
   *  \brief Wrapper converting a filename to an std::ifstream for use with the SWIG
   */
  class istream_proxy
  {
  private:
    std::ifstream fStream; //!< Internal ifstream object
  public:
    istream_proxy(std::string const& filename): //!< Constructor takes a filename and initialises the internal ifstream
      fStream(filename.c_str()) {}

    std::istream& stream() {return fStream;} //!< Returns the internal ifstream
  };

  /*!
   *  \class ostream_proxy
   *  \brief Wrapper converting a filename to an std::ofstream for use with the SWIG
   */
  class ostream_proxy
  {
  private:
    std::ofstream fStream; //!< Internal ofstream object
  public:
    ostream_proxy(std::string const& filename): //!< Constructor takes a filename and initialises the internal ofstream
      fStream(filename.c_str()) {}

    std::ostream& stream() {return fStream;} //!< Returns the internal stream
  };

  // *******************  Numerical *****************************
  double integrate(double data[], size_t npoints, double h); //!< Basic simpson rule integrator

  // ******************* std::string Tokenisers **********************
  std::vector<std::string> split(std::string const& input);   //!< Split std::strings, return std::string std::vector.
  void split(std::vector<std::string>& results, std::string const& input);

  std::vector<real> rsplit(std::string const& input);  //!< Split std::strings, return real std::vector.
  void rsplit(std::vector<real>& results, std::string const& input);

  std::vector<int> isplit(std::string const& input);  //!< Split std::strings, return int std::vector.
  void isplit(std::vector<int>& results, std::string const& input);

  // ******************* Basic Stat Functions *********************
  real ComputeAVG(int const& n, const real *x);             //!< Compute average from x points
  real ComputeAVG(std::vector<real> const& x);              //!< Compute average from std::vector<double> x
  real ComputeStdDev(int const& n, const real *x);          //!< Compute the std deviation
  real ComputeStdDev(std::vector<real> const& x);           //!< Compute the std deviation
  real ComputeEigErr(int const& p, const real *x);          //!< Compute error in the Hessian method for PDFs
  real ComputeSymEigErr(int const& p, const real *x);       //!< Compute error in the symmetric Hessian method for PDFs
  real ComputeMom(int const& n, const real *x, int const& m);//!< Compute mth moment of distribution
  void Compute68cl(std::vector<real> const& x, real &up, real &dn);//!< Compute the 68% c.l.
  void Compute95cl(std::vector<real> const& x, real &up, real &dn);//!< Compute the 95% c.l.

  void CholeskyDecomposition(int const& n, double** const inmatrix, double** sqrtmat);

}

/*! @} */
