%module(package="NNPDF") experiments
 %{
#define SWIG_FILE_WITH_INIT
#include "NNPDF/exceptions.h"
#include "NNPDF/experiments.h"
/*This needs to be here because we need to
 * be able to access the full definition of ThPredictions at compile
 * time. For some reason, an %import doesn't do the trick.
 */
#include "NNPDF/thpredictions.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 

%include "include/numpy.i"

%init %{
    import_array();
%}

%import "thpredictions.i"

%import "dataset.i"
%import "pdfset.i"
%import "utils.i"

/* Parse the header file to generate wrappers */

%ignore std::vector<NNPDF::DataSet>::vector(size_type);
%ignore std::vector<NNPDF::DataSet>::resize;
%template (vector_dataset) std::vector<NNPDF::DataSet>;

%ignore std::vector<NNPDF::ThPredictions>::vector(size_type);
%ignore std::vector<NNPDF::ThPredictions>::resize;
%template (vector_thpredictions) std::vector<NNPDF::ThPredictions>;


%template (vector_experiment_pointer) std::vector<NNPDF::Experiment*>;

%apply (double** ARGOUTVIEWM_ARRAY2, int* DIM1, int* DIM2) {(double** datamat, int* n, int* m)}
%apply (double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** cv, int* n)}

%feature("autodoc", "3");

%include "include/excepthandler.i"

%include "NNPDF/experiments.h"

%feature("docstring") NNPDF::Experiment::get_covmat
"Return a copy of the experiment covariance matrix."
%feature("docstring") NNPDF::Experiment::get_cv
"Return a copy of the central values for the experiment."

%extend NNPDF::Experiment{

void get_covmat(double ** datamat, int* n, int* m){
    int len = $self->GetNData();
    auto data = $self->GetCovMat().data();
    auto size = len*len;
    auto result = (double*) malloc(sizeof(double)*size);
    std::copy(data, data + size, result);
    *datamat = result;
    *m = *n = len;
}


void get_sqrtcovmat(double ** datamat, int* n, int* m){
    int len = $self->GetNData();
    auto data = $self->GetSqrtCov().data();
    auto size = len*len;
    auto result = (double*) malloc(sizeof(double)*size);
    std::copy(data, data + size, result);
    *datamat = result;
    *m = *n = len;
}

void get_cv (double **cv, int* n){
    int len = $self->GetNData();
    auto result = (double*) malloc(sizeof(double)*len);
    const double *data = self->GetData();
    for (int i = 0; i < len; i++){
        result[i] = data[i];
    }
    *cv = result;
    *n = len;
}

}
