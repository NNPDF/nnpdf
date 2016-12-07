#ifdef SSE_CONV
#define REALDOUBLE 0
%numpy_typemaps(NNPDF::real, NPY_FLOAT , int)
#else
#define REALDOUBLE 1
%numpy_typemaps(NNPDF::real, NPY_DOUBLE , int)
#endif
