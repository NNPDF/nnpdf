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
#include <fcntl.h>

#include "NNPDF/utils.h"
#include "NNPDF/exceptions.h"

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"

#include <archive.h>
#include <archive_entry.h>

namespace NNPDF
{
std::string joinpath(const std::initializer_list<std::string> &list)
{
  return joinpath_inner(list);
}
  //__________________________________________________________________
  class archive_wrapper_read
  {
  public:
    archive_wrapper_read(): a(archive_read_new()) {}
    operator struct archive*() {return a;}
    ~archive_wrapper_read()
    {
      archive_read_free(a);
    }
  private:
    struct archive * a;
  };

  //__________________________________________________________________
  class archive_wrapper_write
  {
  public:
    archive_wrapper_write(): a(archive_write_new()) {}
    operator struct archive*() {return a;}
    ~archive_wrapper_write()
    {
      archive_write_finish_entry(a);
      archive_write_free(a);
    }
  private:
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
  void write_to_file(std::string const& filename, std::string const& data)
  {
    int file = creat(filename.c_str(), S_IRUSR | S_IWUSR);
    if (file == -1)
      throw FileError("write_to_file::creat", "Error creating file " +
                      filename + ". Errno = " + std::string(strerror(errno)));

    int wr_err = write(file, data.c_str(), data.size());
    if (wr_err == -1)
      throw FileError("write_to_file::write", "Error writing to file " +
                      filename + ". Errno = " + std::string(strerror(errno)));

    int fs_err = fsync(file);
    if (fs_err == -1)
      throw FileError("write_to_file::write", "Error flushing file " +
                      filename + ". Errno = " + std::string(strerror(errno)));

    int close_err = close(file);
    if (close_err == -1)
      throw FileError("write_to_file::write", "Error closing file " +
                      filename + ". Errno = " + std::string(strerror(errno)));
  }

  //__________________________________________________________________
  std::string buf_from_file(std::string const & filename){
     std::string buf{};
     std::ifstream is(filename.c_str());
     if (is.fail()){
       throw RuntimeException("utils", "Could not open " + filename);
     }

     is.seekg(0, std::ios_base::end);
     const auto fileSize = static_cast<size_t>(is.tellg());
     buf.resize(fileSize);

     is.seekg(0, std::ios_base::beg);
     is.read(&buf[0], fileSize);

     is.close();
     return buf;
  }

  std::string untargz(std::string const& filename)
  {
    auto a = archive_wrapper_read{};

    //The lifetime of this is managed the archive struct
    struct archive_entry *entry;
    std::string buf{};

    auto filterr = archive_read_support_filter_gzip(a);

    if (!(filterr == ARCHIVE_OK || filterr == ARCHIVE_WARN))
    {
        throw RuntimeException("untargz", "Failed to setup support for gzip format");
    }

    auto readerr = archive_read_support_format_tar(a);
    if (readerr != ARCHIVE_OK){
        throw RuntimeException("untargz", "Failed to setup support for tar format");
    }

    // if these operations fail, most likely you have a plain txt file.
    // Apparently libarchive, when supporting all the formats can, at random, either fail to
    // recognize the format or fail to open the header. In fact, this is different for different
    // FKTables.
    // Could not understand the format
    if (archive_read_open_filename(a, filename.c_str(), 10240) != ARCHIVE_OK){
        buf = buf_from_file(filename);
    // Could not get the header
    } else if (archive_read_next_header(a, &entry) != ARCHIVE_OK)
      {
        buf = buf_from_file(filename);
      }
    else
      {
        // get the entry size
        auto entry_size = archive_entry_size(entry);
        if (entry_size == 0)
          throw RuntimeException("untargz", "Could not obtain the uncompressed size from the header in " + filename);

        // read buffer
        buf.resize(entry_size);
        char throwaway;
        if ((archive_read_data(a, &buf[0], entry_size) != entry_size) ||
            (archive_read_data(a, &throwaway, 1) != 0))
          throw RuntimeException("untargz", archive_error_string(a));
      }

    return buf;
  }

  //____________________________________________________________________
  void targz(std::string const& filename, std::string const& data)
  {
    auto a = archive_wrapper_write{};
    if (a == NULL)
      throw RuntimeException("targz", "Empty archive write.");

    // Allocate the compressors and file
    if ((archive_write_add_filter_gzip(a)  != ARCHIVE_OK) ||
        (archive_write_set_format_ustar(a) != ARCHIVE_OK) ||
        (archive_write_open_filename(a, filename.c_str()) != ARCHIVE_OK))
        throw RuntimeException("targz", "Cannot allocate compressed file " + filename);

    // Prepare entry
    auto entry = archive_entry_wrapper{};
    if (entry == NULL)
      throw RuntimeException("targz", "Empty archive entry");

    // Some options
    archive_entry_set_pathname(entry, filename.c_str());
    archive_entry_set_perm(entry, 0644);
    archive_entry_set_filetype(entry, AE_IFREG);
    archive_entry_set_size(entry, data.size());

    // Write header and data
    if (archive_write_header(a, entry) != ARCHIVE_OK)
      throw RuntimeException("targz", "Cannot write header.");

    if (archive_write_data(a, data.c_str(), data.size()) != (int) data.size())
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
    sprintf(buffer, "%s", input.c_str());
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
    sprintf(buffer, "%s", input.c_str());
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

}
