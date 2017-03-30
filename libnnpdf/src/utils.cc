// $Id: utils.cc 2825 2015-05-03 09:14:46Z stefano.carrazza@mi.infn.it $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#include <cmath>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include <iterator>
#include <iostream> 
#include <limits>
#include <algorithm>

#include "NNPDF/utils.h"
#include "NNPDF/exceptions.h"

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"

#include <archive.h>
#include <archive_entry.h>

namespace NNPDF
{
  //__________________________________________________________________
  class archive_wrapper
  {
  public:
    enum mode {read, write};

    archive_wrapper(mode const& m):
      fmode(m),
      a(fmode == read ? archive_read_new() : archive_write_new()) {}

    operator struct archive*() {return a;}
    ~archive_wrapper()
    {
      if (fmode == read)
        archive_read_free(a);
      else
        {
          archive_write_finish_entry(a);
          archive_write_free(a);
        }
    }

  private:
    mode fmode;
    struct archive * a;
  };

  //__________________________________________________________________
  class archive_entry_wrapper
  {
  public:
    archive_entry_wrapper(): entry(archive_entry_new()) {}
    operator struct archive_entry*() {return entry;}
    ~archive_entry_wrapper()
    {
      archive_entry_free(entry);
    }

  private:
    struct archive_entry * entry;
  };

  //__________________________________________________________________
  std::vector<char> untargz(std::string const& filename)
  {
    auto a = archive_wrapper{archive_wrapper::read};
    struct archive_entry *entry;

    // allocate tar.gz decompressor
    if ((archive_read_support_filter_all(a) != ARCHIVE_OK) ||
        (archive_read_support_format_all(a) != ARCHIVE_OK) ||
        (archive_read_open_filename(a, filename.c_str(), 10240)))
        throw RuntimeException("untargz", "Cannot open file");

    std::vector<char> buf{};

    // if this operation fails, most likely you have a plain txt file.
    if (archive_read_next_header(a, &entry) != ARCHIVE_OK)
      {
        std::ifstream is(filename.c_str());
        if (is.fail())
          throw RuntimeException("untargz", "File not found " + filename);

        is.seekg(0, std::ios_base::end);
        const auto fileSize = size_t{is.tellg()};
        buf.resize(fileSize);

        is.seekg(0, std::ios_base::beg);
        is.read(&buf[0], fileSize);
      }
    else
      {
        // get the entry size
        auto entry_size = archive_entry_size(entry);
        if (entry_size == 0)
          throw RuntimeException("untargz", "Compression algorithm not enabled.");

        // read buffer
	char *throwaway = new char[entry_size];
	if ((archive_read_data(a, throwaway, entry_size) != entry_size) ||
	    (archive_read_data(a, throwaway, 1) != 0))
	  throw RuntimeException("untargz", archive_error_string(a));

	buf.assign(throwaway, throwaway + entry_size);
	delete[] throwaway;
      }

    return buf;
  }

  //____________________________________________________________________
  void targz(std::string const& filename, std::stringstream const& data)
  {    
    const auto strdata = data.str();
    auto a = archive_wrapper{archive_wrapper::write};
    if (a == NULL)
      throw RuntimeException("targz", "Empty archive write.");

    // Allocate the compressors and file
    if ((archive_write_add_filter_gzip(a)  != ARCHIVE_OK) ||
        (archive_write_set_format_ustar(a) != ARCHIVE_OK) ||
        (archive_write_open_filename(a, filename.c_str()) != ARCHIVE_OK))
        throw RuntimeException("targz", "Cannot allocate compressed file.");

    // Prepare entry
    auto entry = archive_entry_wrapper{};
    if (entry == NULL)
      throw RuntimeException("targz", "Empty archive entry");

    // Setup write file options
    struct timespec ts;
    if (clock_gettime(CLOCK_REALTIME, &ts) != 0)
      throw RuntimeException("targz", "Error retrieving clock time.");

    // Some options
    archive_entry_set_pathname(entry, filename.c_str());
    archive_entry_set_perm(entry, 0644);
    archive_entry_set_filetype(entry, AE_IFREG);
    archive_entry_set_size(entry, strdata.size());
    archive_entry_set_atime(entry, ts.tv_sec, ts.tv_nsec);
    archive_entry_set_birthtime(entry, ts.tv_sec, ts.tv_nsec);
    archive_entry_set_ctime(entry, ts.tv_sec, ts.tv_nsec);
    archive_entry_set_mtime(entry, ts.tv_sec, ts.tv_nsec);

    // Write header and data
    if (archive_write_header(a, entry) != ARCHIVE_OK)
      throw RuntimeException("targz", "Cannot write header.");

    if (archive_write_data(a, strdata.c_str(), strdata.size()) != (int) strdata.size())
      throw RuntimeException("targz", "Written length does not match data length");
  }

  // /very/ basic integrator
  double integrate(double data[], size_t npoints, double h)
  {
    double integral=0;
    
    integral+=data[0]+data[npoints-1];
    
    for ( size_t j=1; j<(npoints)/2 ; j++ )
      integral+=2*data[2*j -1];
    
    for (size_t j=1; j<(npoints)/2 + 1; j++)
      integral+=4*data[2*j - 2];
    
    return integral*h/3.0;
  }

  /**
   * Split string into string vector
   */
  std::vector<std::string> split(std::string const &input)
  {
  	std::stringstream strstr(input);
  	std::istream_iterator<std::string> it(strstr);
  	std::istream_iterator<std::string> end;
  	std::vector<std::string> results(it, end);
  	return results;
  }

  void split(std::vector<std::string>& results, std::string const& input)
  {
  	std::stringstream strstr(input);
  	std::istream_iterator<std::string> it(strstr);
  	std::istream_iterator<std::string> end;
    
  	results.assign(it, end);
  	return;
  }

  /**
   * Split std::string into real std::vector
   */
  std::vector<real> rsplit(std::string const& input)
  {
  	std::vector<real> results;
    char *buffer = new char[input.size() + 1];
    sprintf(buffer, input.c_str());
    char *token = strtok(buffer, " \t");
    while (token)
      {
        results.push_back(atof(token));
        token = strtok(NULL, " \t");
      }  
    delete[] buffer;
  	return results;
  }

  void rsplit(std::vector<real>& results, std::string const& input)
  {
    results.clear();
    char *buffer = new char[input.size() + 1];
    sprintf(buffer, input.c_str());
    char *token = strtok(buffer, " \t");
    while (token)
      {
        results.push_back(atof(token));
        token = strtok(NULL, " \t");
      }  
    delete[] buffer;
  	return;
  }

  /**
   * Split std::string into integer std::vector
   */

  std::vector<int> isplit(std::string const& input)
  {
  	std::stringstream strstr(input);
  	std::istream_iterator<int> it(strstr);
  	std::istream_iterator<int> end;
  	std::vector<int> results(it, end);
  	return results;
  }

  void isplit(std::vector<int>& results, std::string const& input)
  {
  	std::stringstream strstr(input);
  	std::istream_iterator<int> it(strstr);
  	std::istream_iterator<int> end;
    
  	results.assign(it, end);
  	return;
  }

  /**
    * Compute average
    * \param n number of points
    * \param x array with values
    * \return the average as real
    */
  real ComputeAVG(int const& n, const real *x)
  {
    real sum = 0.0;
    for (int i = 0; i < n; i++) 
      {
        sum += x[i];
      }

    return sum / n;
  }

  /**
    * Compute average
    * \param x std::vector<real> with values
    * \return the average as real
    */
  real ComputeAVG(std::vector<real> const& x)
  {
    if (x.size() != 0)
      {
        int n = (int) x.size();
        real sum = 0.0;
        for (int i = 0; i < n; i++)
  	sum += x[i];

        return sum / n;
      }
    
    return 0;
  }

  /**
    * Compute the standard deviation
    * \param n number of points
    * \param x array with values
    * \return the std dev as real
    */
  real ComputeStdDev(int const& n, const real *x)
  {
    real sum = 0.0;
    real avg = ComputeAVG(n, x);
    for (int i = 0; i < n; i++)
      sum += (x[i]-avg)*(x[i]-avg);

    sum /= n-1;

    return sqrt(sum);
  }

  /**
    * Compute the standard deviation
    * \param x std::vector<real> with values
    * \return the std dev as real
    */
  real ComputeStdDev(std::vector<real> const& x)
  {
    if (x.size() != 0)
      {
        real sum = 0.0;
        int n = (int) x.size();
        real avg = ComputeAVG(x);
        for (int i = 0; i < n; i++)
  	sum += (x[i]-avg)*(x[i]-avg);

        sum /= n-1;

        return sqrt(sum);
      }

    return 0;
  }

  /**
   * Compute the 68% c.l.
   */
  void Compute68cl(std::vector<real> const& x, real &up, real &dn)
  {
    up = 0;
    dn = 0;
    if (x.size() > 0)
      {
        int esc = (int) (x.size()*(1-0.68)/2);
        std::vector<real> xval(x);
        std::sort(xval.begin(),xval.end());
        up = xval[xval.size()-1-esc];
        dn = xval[esc];	
      }
  }

  /**
   * Compute the 95% c.l.
   */
  void Compute95cl(std::vector<real> const& x, real &up, real &dn)
  {
    up = 0;
    dn = 0;
    if (x.size() > 0)
      {
        int esc = (int) (x.size()*(1-0.95)/2);
        std::vector<real> xval(x);
        std::sort(xval.begin(),xval.end());
        up = xval[xval.size()-1-esc];
        dn = xval[esc];	
      }
  }

  /**
    * Compute the errors for the Hessian method
    * \param p number of pdfs
    * \param x array with values
    * \return the error for the Hessian method as real
    */
  real ComputeEigErr(int const& p, const real *x)
  {
    real err = 0;
    const int nvec = (p-1)/2.0;

    for (int i = 0; i < nvec; i++)
      err += pow(x[2*i+1]-x[2*i+2], 2); // Eigenstd::vector

    return sqrt(err)/2.0;
  }

  /**
    * Compute the errors for the symmetric Hessian method
    * \param p number of pdfs
    * \param x array with values
    * \return the error for the Hessian method as real
    */
  real ComputeSymEigErr(int const& p, const real *x)
  {
    real err = 0;
    for (int i = 1; i < p; i++)
      err += pow(x[i]-x[0], 2); // Eigenstd::vector

    return sqrt(err);
  }


  /**
   * Compute the mth moment of the distribution
   * \param n number of points
   * \param x array with values
   * \return the std dev as real
   */
  real ComputeMom(int const& n, const real *x, int const& m)
  {
    real sum = 0.0;
    real avg = ComputeAVG(n, x);
    for (int i = 0; i < n; i++)
      sum += pow(x[i]-avg,m);
    
    sum /= n;
    
    return sum;
  }

  // *********** Cholesky decomposition of a matrix ***************

  void CholeskyDecomposition(int const& n, double** const inmatrix, double** sqrtmat)
  {
    if (n <= 0)
      throw LengthError("CholeskyDecomposition","attempting a decomposition of an empty matrix!");
    gsl_matrix* mat = gsl_matrix_calloc(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        gsl_matrix_set(mat, i, j, inmatrix[i][j]);

    const int decomp = gsl_linalg_cholesky_decomp(mat);
    for (int i = 0; i < n; i++)
      for (int j = 0; j <= i; j++)
        sqrtmat[i][j] =  gsl_matrix_get(mat, i, j);
    gsl_matrix_free (mat);
  }


}
