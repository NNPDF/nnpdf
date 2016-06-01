%module(package="NNPDF") dataset
 %{
#define SWIG_FILE_WITH_INIT
#include "../../src/NNPDF/exceptions.h"
#include "../../src/NNPDF/dataset.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 

%include "include/numpy.i"

%init %{
    import_array();
%}

%import "commondata.i"
%import "fkset.i"
%import "pdfset.i"
/* Parse the header file to generate wrappers */

%feature("autodoc", "3");

%include "include/excepthandler.i"

%template(vector_int) std::vector<int>;

/* We copy the arrays for every reason. It's too dangerous to pass by
 * reference with too little benefict. */
%apply (double** ARGOUTVIEWM_ARRAY2, int* DIM1, int* DIM2) {(double** datamat, int* n, int* m)}
%apply (double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** cv, int* n)}

%ignore NNPDF::swap;
%ignore NNPDF::DataSet::operator=;
%ignore NNPDF::DataSet::DataSet(DataSet &&);
%include "../../src/NNPDF/dataset.h"

%feature("docstring") NNPDF::DataSet::get_covmat
"Return a copy of the experiment covariance matrix."
%feature("docstring") NNPDF::DataSet::get_cv
"Return a copy of the central values for the experiment."

%extend NNPDF::DataSet{

void get_covmat(double ** datamat, int* n, int* m){
    int len = $self->GetNData();
    double** cov = $self->GetCovMat();
    double * result = new double[len*len];
    for (int i = 0; i < len; i++){
        for (int j = 0; j < len; j++){
            result[len*i + j] = cov[i][j];
        }
    }
    *datamat = result;
    *m = *n = len;
}


void get_invcovmat(double ** datamat, int* n, int* m){
    int len = $self->GetNData();
    double** invcov = $self->GetInvCovMat();
    double * result = new double[len*len];
    for (int i = 0; i < len; i++){
        for (int j = 0; j < len; j++){
            result[len*i + j] = invcov[i][j];
        }
    }
    *datamat = result;
    *m = *n = len;
}

void get_cv (double **cv, int* n){
    int len = $self->GetNData();
    double * result = new double[len];
    double *data = self->GetData();
    for (int i = 0; i < len; i++){
        result[i] = data[i];
    }
    *cv = result;
    *n = len;
}

%pythoncode{

def __len__(self):
    return self.GetNData();

}

}
