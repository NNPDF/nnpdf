#ifdef SSE_CONV
%numpy_typemaps(NNPDF::real, NPY_FLOAT , int)
#else
%numpy_typemaps(NNPDF::real, NPY_DOUBLE , int)
#endif
