// $Id: utils.cc 2825 2015-05-03 09:14:46Z stefano.carrazza@mi.infn.it $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#include <cmath>
#include <cstdlib>
#include <sstream>
#include <iterator>
#include <iostream> 
#include <limits>
#include <algorithm>

#include "NNPDF/utils.h"
#include "NNPDF/exceptions.h"

namespace NNPDF
{

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
  	std::stringstream strstr(input);
  	std::istream_iterator<real> it(strstr);
  	std::istream_iterator<real> end;
  	std::vector<real> results(it, end);
  	return results;
  }

  void rsplit(std::vector<real>& results, std::string const& input)
  {
  	std::stringstream strstr(input);
  	std::istream_iterator<real> it(strstr);
  	std::istream_iterator<real> end;
    
  	results.assign(it, end);
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


  // *************************** INVERSION METHODS *******************************

  /**
   * @brief GetMatrixArray
   * @param n
   * @param covmat
   * @return
   */
  static double *GetMatrixArray(int const& n, double **covmat)
  {
    double *array = new double[n*n];

    int index = 0;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        {
          array[index] = covmat[i][j];
          index++;
        }

    return array;
  }

  /**
   * @brief SetMatrixArray
   * @param n
   * @param array
   * @param invcovmat
   */
  static void SetMatrixArray(int n, double *array, double **invcovmat)
  {
    int index = 0;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        {
          invcovmat[i][j] = array[index];
          index++;
        }
  }

  static bool DecomposeLUCrout(int n, double **lu, double **out, int *index, double &sign, double tol, int &nrZeros)
  {
    // Crout/Doolittle algorithm of LU decomposing a square matrix, with implicit partial
    // pivoting.  The decomposition is stored in fLU: U is explicit in the upper triag
    // and L is in multiplier form in the subdiagionals .
    // Row permutations are mapped out in fIndex. fSign, used for calculating the
    // determinant, is +/- 1 for even/odd row permutations. .

    double *pLU = GetMatrixArray(n, lu);
    double *scale = new double[n];

    sign = 1.0;
    nrZeros = 0;

    // Find implicit scaling factors for each row
    for (int i = 0; i < n; i++)
      {
        const int off_i = i*n;
        double max = 0.0;
        for (int j = 0; j < n; j++)
          {
            const double tmp = fabs(pLU[off_i + j]);
            if (tmp > max) max = tmp;
          }
        scale[i] = (max == 0.0 ? 0.0 : 1.0/max);
      }

    for (int j = 0; j < n; j++)
      {
        const int off_j = j*n;
        for (int i = 0; i < j; i++)
          {
            const int off_i = i*n;
            double r = pLU[off_i + j];
            for (int k = 0; k < i; k++)
              {
                const int off_k = k*n;
                r -= pLU[off_i + k]*pLU[off_k + j];
              }
            pLU[off_i + j] = r;
          }
        // Run down jth subdiag to form the residuals after the elimination of
        // the first j-1 subdiags.  These residuals divided by the appropriate
        // diagonal term will become the multipliers in the elimination of the jth.
        // subdiag. Find fIndex of largest scaled term in imax.

        double max = 0.0;
        int imax = 0;
        for (int i = j; i < n; i++)
          {
            const int off_i = i*n;
            double r = pLU[off_i + j];
            for (int k = 0; k < j; k++)
              {
                const int off_k = k*n;
                r -= pLU[off_i + k]*pLU[off_k + j];
              }
            pLU[off_i + j] = r;
            const double tmp = scale[i]*fabs(r);
            if (tmp >= max)
              {
                max = tmp;
                imax = i;
              }
          }

        // Permute current row with imax
        if (j != imax)
          {
            const int off_imax = imax*n;
            for (int k = 0; k < n; k++)
              {
                const double tmp = pLU[off_imax + k];
                pLU[off_imax + k] = pLU[off_j + k];
                pLU[off_j + k] = tmp;
              }
            sign = -sign;
            scale[imax] = scale[j];
          }
        index[j] = imax;

        // If diag term is not zero divide subdiag to form multipliers.
        if (pLU[off_j + j] != 0.0)
          {
            if (fabs(pLU[off_j + j]) < tol)
              nrZeros++;
            if (j != n-1)
              {
                const double tmp = 1.0/pLU[off_j + j];
                for (int i = j+1; i < n; i++)
                  {
                    const int off_i = i*n;
                    pLU[off_i+j] *= tmp;
                  }
              }
          }
        else
          {
            delete[] scale;
            throw EvaluationError("DecomposeLUCrout","matrix is singular");
          }
      }

    delete[] scale;

    SetMatrixArray(n, pLU, out);

    delete[] pLU;

    return true;
  }

  /**
   * @brief Invert
   * @param n
   * @param covmat
   * @param invcovmat
   */
  bool InvertLU(int n, double **covmat, double **invcovmat)
  {
    // Check for existence
    if (n == 0)
      {
        std::cerr << "Error inversion: Matrix dim should be > 0" << std::endl;
        return false;
      }

    int *index = new int[n];
    
    // Check for singular matrix
    double sign = 1.0;
    int nrZeros = 0;
    double tol = std::numeric_limits<double>::epsilon();
    if (!DecomposeLUCrout(n, covmat, invcovmat, index, sign, tol, nrZeros) || nrZeros > 0)
      throw EvaluationError("InvertLU","Error inversion: Matrix is singular");

    double *pLU = GetMatrixArray(n,invcovmat);

    // Form inv(U)
    for (int j = 0; j < n; j++)
      {
        const int off_j = j*n;

        pLU[off_j + j] = 1./pLU[off_j + j];
        const double mLU_jj = -pLU[off_j + j];

        // Compute elements 0:j-1 of j-th column
        double *pX = pLU + j;

        for (int k = 0; k < j; k++)
          {
            const int off_k = k*n;
            if ( pX[off_k] != 0.0)
              {
                const double tmp = pX[off_k];
                for (int i = 0; i < k; i++)
                  {
                    const int off_i = i*n;
                    pX[off_i] += tmp*pLU[off_i + k];
                  }
                pX[off_k] *= pLU[off_k + k];
              }
          }
        for (int k = 0; k < j; k++)
          {
            const int off_k = k*n;
            pX[off_k] *= mLU_jj;
          }
      }

    // Solve the equation inv(A)*L = inv(U) for inv(A)
    double *pWorkd = new double[n];
    for (int j = n-1; j >= 0; j--)
      {
        // copy current column j of L to work and replace with zeros
        for (int i = j+1; i < n; i++)
          {
            const int off_i = i*n;
            pWorkd[i] = pLU[off_i + j];
            pLU[off_i + j] = 0.0;
          }

        // Compute current column of inv(A)

        if (j  < n-1)
          {
            const double *mp = pLU + j + 1; // Matrix row ptr
            double *tp = pLU + j;           // Target std::vector ptr

            for (int irow = 0; irow < n; irow++)
              {
                double sum = 0.0;
                const double *sp = pWorkd + j + 1; // Source std::vector ptr
                for (int icol = 0; icol < n-1-j; icol++)
                  sum += *mp++ * *sp++;
                *tp = -sum + *tp;
                mp += j+1;
                tp += n;
              }
          }
      }

    delete[] pWorkd;

    // Apply column interchanges
    for (int j = n-1; j >= 0; j--)
      {
        const int jperm = index[j];
        if (jperm != j)
          {
            for (int i = 0; i < n; i++)
              {
                const int off_i = i*n;
                const double tmp = pLU[off_i + jperm];
                pLU[off_i + jperm] = pLU[off_i + j];
                pLU[off_i + j] = tmp;
              }
          }
      }

    delete[] index;

    // Recombine in invcovmat
    SetMatrixArray(n, pLU, invcovmat);

    delete[] pLU;

    return true;
  }

}
