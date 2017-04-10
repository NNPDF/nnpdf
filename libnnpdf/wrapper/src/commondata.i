%module(package="NNPDF") commondata
 %{
#define SWIG_FILE_WITH_INIT
#include "NNPDF/exceptions.h"
#include "NNPDF/commondata.h"
 %}

%include "std_string.i" 
%include "std_vector.i" 
%include "std_map.i" 

%include "include/numpy.i"

%init %{
    import_array();
%}
%include "common.i"
%template(map_str_vector_str) std::map<std::string, std::vector<std::string> >;
/* Parse the header file to generate wrappers */

%feature("autodoc", "3");

%include "include/excepthandler.i"

%apply (double** ARGOUTVIEWM_ARRAY2, int* DIM1, int* DIM2) {(double** datamat, int* n, int* m)}
%apply (double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** cv, int* n)}

%ignore NNPDF::swap;
%ignore NNPDF::CommonData::operator=;
%ignore NNPDF::CommonData::CommonData(CommonData &&);
%include "NNPDF/commondata.h"

%extend NNPDF::CommonData{

void get_kintable(double** datamat, int* n, int* m){
    int len = $self->GetNData();
    double* result = (decltype(result)) malloc(sizeof(*result)*len*3);
    for (int i = 0; i < len; i++){
        int pos = 3*i;
        result[pos] = $self->GetKinematics(i, 0);
        result[pos + 1] = $self->GetKinematics(i, 1);
        result[pos + 2] = $self->GetKinematics(i, 2);
    }
    *n = len;
    *m = 3;
    *datamat = result;
}

void get_cv (double **cv, int* n){
    auto len = $self->GetNData();
    auto * result = (double*) malloc(sizeof(double)*len);
    auto *data = self->GetData();
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

