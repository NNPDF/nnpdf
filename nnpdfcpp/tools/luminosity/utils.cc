/**
 * $Id$
 * Author: Stefano Carrazza, stefano.carrazza@mi.infn.it
 */

#include <cmath>
#include <cstdlib>
#include <sstream>
#include <iterator>

#include "utils.h"

/**
 * Split string int string vector
 */
vector<string> split(string &input)
{
	stringstream strstr(input);
	istream_iterator<string> it(strstr);
	istream_iterator<string> end;
	vector<string> results(it, end);
	return results;
}

/**
 * Split string int double vector
 */
vector<double> rsplit(string &input)
{
	stringstream strstr(input);
	istream_iterator<double> it(strstr);
	istream_iterator<double> end;
	vector<double> results(it, end);
	return results;
}

vector<int> isplit(string &input)
{
	stringstream strstr(input);
	istream_iterator<int> it(strstr);
	istream_iterator<int> end;
	vector<int> results(it, end);
	return results;
}

/**
  * Compute average
  * \param n number of points
  * \param x array with values
  * \return the average as double
  */
double ComputeAVG(int n, const double *x)
{
  double sum = 0.0;
  for (int i = 0; i < n; i++)
    sum += x[i];

  return sum / n;
}

/**
  * Compute average
  * \param n number of data points
  * \param pdf number of pdf
  * \param x array with values
  * \return the average as double
  */
double ComputeAVG(int n, int pdf, const double *x)
{
  double sum = 0.0;
  for (int i = 0; i < pdf; i++)
    sum += x[i*n];

  return sum / pdf;
}

/**
  * Compute average
  * \param x vector<double> with values
  * \return the average as double
  */
double ComputeAVG(vector<double> x)
{
  int n = (int) x.size();
  double sum = 0.0;
  for (int i = 0; i < n; i++)
    sum += x[i];

  return sum / n;
}

/**
  * Compute the standard deviation
  * \param n number of points
  * \param x array with values
  * \return the std dev as double
  */
double ComputeStdDev(int n, const double *x)
{
  double sum = 0.0;
  double avg = ComputeAVG(n, x);
  for (int i = 0; i < n; i++)
    sum += (x[i]-avg)*(x[i]-avg);

  sum /= n-1;

  return sqrt(sum);
}

/**
  * Compute the standard deviation
  * \param n number of points
  * \param pdf number of pdf
  * \param x array with values
  * \return the std dev as double
  */
double ComputeStdDev(int n, int pdf, const double *x)
{
  double sum = 0.0;
  double avg = ComputeAVG(n, pdf, x);
  for (int i = 0; i < pdf; i++)
    sum += (x[i*n]-avg)*(x[i*n]-avg);

  sum /= pdf-1;

  return sqrt(sum);
}

/**
  * Compute the standard deviation
  * \param x vector<double> with values
  * \return the std dev as double
  */
double ComputeStdDev(vector<double> x)
{
  double sum = 0.0;
  int n = (int) x.size();
  double avg = ComputeAVG(x);
  for (int i = 0; i < n; i++)
    sum += (x[i]-avg)*(x[i]-avg);

  sum /= n-1;

  return sqrt(sum);
}


/**
  * Compute the errors for the Hessian method
  * \param p number of pdfs
  * \param x array with values
  * \return the error for the Hessian method as double
  */
double ComputeEigErr(int p, const double *x)
{
  double err = 0;
  const int nvec = p/2;

  for (int i = 0; i < nvec; i++)
    err += pow(x[2*i]-x[2*i+1], 2); // Eigenvector

  return sqrt(err)/2.0;
}

/**
  * Compute the errors for the Hessian method
  * \param p number of pdfs
  * \param n number of points
  * \param x array with values
  * \return the error for the Hessian method as double
  */
double ComputeEigErr(int p, int n, const double *x)
{
  double err = 0;
  const int nvec = p/2;

  for (int i = 0; i < nvec; i++)
    err += pow(x[n*2*i]-x[n*(2*i+1)], 2.0); // Eigenvector

  return sqrt(err)/2.0;
}


/**
  * Compute the errors for the Hessian method
  * \param p number of pdfs
  * \param x array with values
  * \return the error for the Hessian method as double
  */
double ComputeSymEigErr(int p, double cv, const double *x)
{
  double err = 0;

  for (int i = 0; i < p; i++)
    err += pow(x[i]-cv, 2.0); // Eigenvector

  return sqrt(err);
}
