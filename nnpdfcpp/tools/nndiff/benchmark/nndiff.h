#include <NNPDF/common.h>
#include <vector>
using NNPDF::real;
using std::vector;

real nnval(real x_00, vector<real> const& params);

real nndiff(real x_00, vector<real> const& params);

real NNPDFdev(real const& x, vector<real> const& params, 
	      real const& a, real const& b, real const& n);
