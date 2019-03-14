// $Id: utils.h 2825 2015-05-03 09:14:46Z stefano.carrazza@mi.infn.it $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#pragma once

#include "common.h"
#include "NNPDF/exceptions.h"
#include <vector>
#include <string>
#include <fstream>
#include <array>
#include <iterator>
#include <initializer_list>

/** @defgroup utils Utils
 * \brief libnnpdf utility functions
 */

/*! \addtogroup utils
 *  @{
 */

namespace NNPDF
{

  // This is all so we can call both joinpath({"a", "b", "c"}) which
  // requires a special treatment and auto v = vector<string>{"a", "b", "c"}; joinpath(v). See
  // https://stackoverflow.com/questions/4757614/why-doesnt-my-template-accept-an-initializer-list#4763493
  namespace
  {
    template <typename T> inline std::string joinpath_inner(const T &list)
    {
      // This is implemented following
      // https://github.com/python/cpython/blob/05d68a8bd84cb141be9f9335f5b3540f15a989c4/Lib/posixpath.py#L75
      const auto sep = '/';
      auto path = std::string{""};

      for (const auto &it : list) {
        //TODO: Use this in C++17
        //if (std::empty(it)) {
        if (it.empty()) {
          continue;
        }
        if (*std::cbegin(it) == sep) {
          path = it;
        } else if (path.empty() || *std::crbegin(path) == sep) {
          path += it;
        } else {
          path += sep + it;
        }
      }
      return path;
    }
  }

  /**
   * @brief Produce a UNIX path from a container of string-like objects.
   * @param A iterable of objects that can be joined with strings and iterated over.
   * @return A string representing a joined path.
   *
   * The resuts should be the same as os.path.join in Python (on POSIX), except
   * that an empty string is returned from an empty iterable.
   */
  template <typename T> std::string joinpath(const T &list)
  {
    return joinpath_inner(list);
  }
  std::string joinpath(const std::initializer_list<std::string> &list);

  /**
   * @brief write_to_file Writes string to file handling failures.
   *
   * This function creates and writes data using a POSIX file
   * descriptor. The data buffer is written to the filesystem through
   * write and fsync methods. The implementation also throws
   * exceptions at multiple layers (opening, writing, flushing and
   * closing the file). This is quite useful on clusters like
   * LXPLUS/EOS where IO issues are non negligible. The downside of
   * this method is the requirement of loading in memory the entire
   * buffer as a string (and obviously the portability to non POSIX
   * systems).
   *
   * @param data the data string
   * @param filename the file name
   */
  void write_to_file(std::string const& filename, const std::string &data);

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
    T&       operator()(size_t i, size_t j)       { return _data[i*_size[1]+j]; }
    T const& operator()(size_t i, size_t j) const { return _data[i*_size[1]+j]; }

    // There is no doubt a better way to do this
    std::vector<T> operator*(std::vector<T> in) const {
        if (_size[0] != in.size())
            throw RangeError("matrix-vector product", "Mismatch of matrix and input vector dimension");
        std::vector<T> out(_size[0],0);
        for (size_t i=0; i<_size[0]; i++)
            for (size_t j=0; j<_size[1]; j++)
                out[i] += _data[i*_size[1]+j]*in[j];
        return out;
    }

    // Data access
    T *       data ()       {return _data.data();}  //!< Return the underlying buffer.
    T const * data () const {return _data.data();}  //!< Return the underlying buffer (const version).

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
}

/*! @} */
