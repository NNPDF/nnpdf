#include <NNPDF/common.h>
#include <vector>
using NNPDF::real;
using std::vector;

real NNPDFval(real const& x_00, vector<real> const& params,
	      real const& a, real const& b, real const& n);

real NNPDFdev(real const& x_00, vector<real> const& params, 
	      real const& a, real const& b, real const& n);

real NNPDFdev2(real const& x_00, vector<real> const& params,
	       real const& a, real const& b, real const& n);

real alphaeff(real const& x_00, vector<real> const& params,
	      real const& a, real const& b, real const& n);

real betaeff(real const& x_00, vector<real> const& params,
	      real const& a, real const& b, real const& n);
