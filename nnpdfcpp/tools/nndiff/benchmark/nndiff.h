#include <NNPDF/common.h>
#include <vector>
using NNPDF::real;
using std::vector;

// n*x^(-a+1)*(1-x)^b*NN(x) 
real NNPDFval(real const& x_00, vector<real> const& params,
	      real const& a, real const& b, real const& n);

// d [n*x^(-a+1)*(1-x)^b*NN(x)] / d x
real NNPDFdev(real const& x_00, vector<real> const& params, 
	      real const& a, real const& b, real const& n);

// d2 [n*x^(-a+1)*(1-x)^b*NN(x)] / dx^2
real NNPDFdev2(real const& x_00, vector<real> const& params,
	       real const& a, real const& b, real const& n);
// NN(x)
real nnval(real const& x, vector<real> const& params);

// d NN(x) / dx
real nnder(real const& x, vector<real> const& params);

// alphaeff
real alphaeff(real const& x_00, vector<real> const& params, real const& a);

// betaeff
real betaeff(real const& x_00, vector<real> const& params, real const& b);
