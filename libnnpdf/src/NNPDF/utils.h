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
#include <array>

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
  std::string untargz(std::string const& filename);

  /**
   * @brief targz Store to disk the compressed data
   * @param filename the output file
   * @param data the stream buffer
   */
  void targz(std::string const& filename, std::string const& data);

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

  /**
   * @brief The matrix class for the covmat related objects
   */
  template<class T>
  class matrix
  {
  public:
    //!< matrix constructor
    matrix(size_t row = 0, size_t col = 0): _size{{row,col}}
    {
      if (row*col != 0)
        _data.resize(row*col);
    }

    //!< resize matrix and fill with v
    void resize(size_t row, size_t col, T v)
    {
      _size = {{row,col}};
      _data.resize(row*col, v);
    }

    //!< clear matrix size and content
    void clear()
    {
      _size = {0,0};
      _data.clear();
    }

    //operators
    size_t const& size(size_t dim) const { return _size[dim]; } //!< Returns the (row,col) size pair.
    T&       operator()(size_t i, size_t j)       { return _data[i+_size[0]*j]; }
    T const& operator()(size_t i, size_t j) const { return _data[i+_size[0]*j]; }
    //TODO: Does this have to be const? In any case there
    //should be a const version.
    T const * data () const {return _data.data();} //!< Return the underlying buffer.

  private:
    std::array<size_t, 2> _size; //!< the dimension pair
    std::vector<T>   _data; //!< the data array
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

  void CholeskyDecomposition(int const& n, matrix<double> const& inmatrix, matrix<double> & sqrtmat);

}

/*! @} */
