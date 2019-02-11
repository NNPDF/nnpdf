%module(package="NNPDF") nnpdf
 %{
#define SWIG_FILE_WITH_INIT
#include "NNPDF/exceptions.h"
#include "NNPDF/common.h"
#include "NNPDF/utils.h"
#include "NNPDF/randomgenerator.h"
#include "NNPDF/pathlib.h"
#include "NNPDF/logger.h"
#include "NNPDF/lhapdfset.h"
#include "LHAPDF/LHAPDF.h"
#include "NNPDF/dataset.h"
#include "NNPDF/experiments.h"
#include "NNPDF/thpredictions.h"
#include "NNPDF/positivity.h"
#include "NNPDF/chisquared.h"
 %}


%include "std_string.i"
%include "std_vector.i"
%include "std_map.i"

%include "include/numpy.i"

%init %{
    import_array();
%}

/* Parse the header file to generate wrappers */

%feature("autodoc", "3");

%include "include/excepthandler.i"
%include "include/real_typemap.i"

%template(vector_str) std::vector<std::string>;
%template(vector_double) std::vector<double>;
%template(vector_int) std::vector<int>;

%typemap(out) NNPDF::matrix<double> {
    auto size = $1.size(0)*$1.size(1);
    /*Have to cast to long explicitly in clang*/
    npy_intp dims[2] = {(long) $1.size(0),  (long) $1.size(1)};
    auto data = (double*) malloc(sizeof(double)*size);
    std::copy($1.data(), $1.data()+size,data);
    $result = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, data);
}

/* Simpler headers without a lot of dependencies and headers */

%include "NNPDF/common.h"
%include "NNPDF/utils.h"

%template(matrix_double) NNPDF::matrix<double>;

%include "NNPDF/randomgenerator.h"
%include "NNPDF/pathlib.h"
%include "NNPDF/logger.h"


%template(map_str_vector_str) std::map<std::string, std::vector<std::string> >;
/* Parse the header file to generate wrappers */

%feature("autodoc", "3");

%include "include/excepthandler.i"

/* Commondata */

/* typemaps for numpy */

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

/* PDF stuff */

%include "NNPDF/pdfset.h"
%rename(single_replica) NNPDF::LHAPDFSet::LHAPDFSet(std::string const &,int const &);

%apply (NNPDF::real** ARGOUTVIEWM_ARRAY4, int* DIM1, int* DIM2, int* DIM3,\
        int* DIM4)\
{(NNPDF::real** datamat, int* nrep_out, int* nf_out, int* nx_out,\
int* nq_out)}

%apply (NNPDF::real* IN_ARRAY1, int DIM1) {(NNPDF::real* xmat,  int nx_in)}
%apply (int*         IN_ARRAY1, int DIM1) {(int*         flmat, int nf_in)}
%apply (NNPDF::real* IN_ARRAY1, int DIM1) {(NNPDF::real* qmat,  int nq_in)}
%include "NNPDF/lhapdfset.h"

%rename(SetLHAPDFVerbosity) LHAPDF::setVerbosity(int v);
namespace LHAPDF{
    int verbosity();
    void setVerbosity(int v);
}

%extend NNPDF::LHAPDFSet{

void grid_values(
                 int* flmat, int nf_in,
                 NNPDF::real* xmat, int nx_in,
                 NNPDF::real* qmat, int nq_in,
                 NNPDF::real** datamat,
                 int* nrep_out, int* nf_out,
                 int* nx_out, int* nq_out
                )
{
    int nrep = self->GetMembers();
    *nrep_out = nrep;
    *nx_out = nx_in;
    *nf_out = nf_in;
    *nq_out = nq_in;
    NNPDF::real* result = (decltype(result)) malloc(sizeof(*result)*nrep*nx_in*nf_in*nq_in);
    for (int irep=0; irep<nrep; irep++){
        for (int ifl=0; ifl<nf_in; ifl++){
            for (int ix=0; ix<nx_in; ix++){
                for(int iq=0; iq<nq_in; iq++){
                    int index = ((irep*nf_in + ifl)*nx_in + ix)*nq_in + iq;
                    result[index] = self->xfxQ(xmat[ix], qmat[iq], irep, flmat[ifl]);
                }
            }
        }
    }
    *datamat = result;
}
}

/* FKtable stuff */

%apply (double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** data, int* n)}
%apply (NNPDF::real** ARGOUTVIEWM_ARRAY1, int* DIM1) {(NNPDF::real** data, int* n)}
%apply (int** ARGOUTVIEWM_ARRAY1, int* DIM1) {(int** data, int* n)}

%include "NNPDF/fastkernel.h"
%template (vector_fktable_p) std::vector<NNPDF::FKTable*>;
%ignore NNPDF::swap;
%ignore NNPDF::FKSet::operator=;
%ignore NNPDF::FKSet::FKSet(FKSet &&);
%include "NNPDF/fkset.h"


%extend NNPDF::FKTable{

void get_xgrid(double **data, int* n){
  int len = $self->GetNx();
  double* result = (double*) malloc(sizeof(double)*len);
  double* xgrid = $self->GetXGrid();
  std::copy(xgrid, xgrid+len, result);
  *data = result;
  *n = len;
}

void get_sigma(NNPDF::real **data, int* n){
  int len = $self->GetDSz()*$self->GetNData();
  NNPDF::real* result = (decltype(result)) malloc(sizeof(*result)*len);
  NNPDF::real* sigma = $self->GetSigma();
  std::copy(sigma, sigma+len, result);
  *data = result;
  *n = len;
}

void get_flmap(int **data, int* n){
  int len = $self->IsHadronic() ? 2*$self->GetNonZero() : $self->GetNonZero();
  int* result = (int*) malloc(sizeof(int)*len);
  int* flmap = $self->GetFlmap();
  std::copy(flmap, flmap+len, result);
  *data = result;
  *n = len;
}

}

/* Dataset */

/* We copy the arrays for every reason. It's too dangerous to pass by
 * reference with too little benefict. */
%apply (double** ARGOUTVIEWM_ARRAY2, int* DIM1, int* DIM2) {(double** datamat, int* n, int* m)}

%ignore NNPDF::swap;
%ignore NNPDF::DataSet::operator=;
%ignore NNPDF::DataSet::DataSet(DataSet &&);
%include "NNPDF/dataset.h"

%feature("docstring") NNPDF::DataSet::get_covmat
"Return a copy of the experiment covariance matrix."
%feature("docstring") NNPDF::DataSet::get_cv
"Return a copy of the central values for the experiment."

%extend NNPDF::DataSet{

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


%pythoncode{

def __len__(self):
    return self.GetNData();

}

}


/*Experiments*/

%ignore std::vector<NNPDF::DataSet>::vector(size_type);
%ignore std::vector<NNPDF::DataSet>::resize;
%template (vector_dataset) std::vector<NNPDF::DataSet>;

%ignore std::vector<NNPDF::ThPredictions>::vector(size_type);
%ignore std::vector<NNPDF::ThPredictions>::resize;
%template (vector_thpredictions) std::vector<NNPDF::ThPredictions>;


%template (vector_experiment_pointer) std::vector<NNPDF::Experiment*>;

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
/*Theory predictions */
%apply (NNPDF::real** ARGOUTVIEWM_ARRAY1, int* DIM1) {(NNPDF::real** data, int* n)}
%apply (NNPDF::real** ARGOUTVIEWM_ARRAY2, int* DIM1, int* DIM2){(NNPDF::real** data, int*m, int* n)}
%ignore NNPDF::ThPredictions::operator=;
%ignore NNPDF::ThPredictions::ThPredictions(ThPredictions &&);
%include "NNPDF/thpredictions.h"
%extend NNPDF::ThPredictions{

void get_data(NNPDF::real **data, int* m, int* n){
    *m = $self->GetNData();
    *n = $self->GetNPdf();
    int len = (*m) * (*n);
    NNPDF::real * result = (decltype(result)) malloc(sizeof(*result)*len);
    NNPDF::real* obs = $self->GetObs();
    std::copy(obs, obs+len ,result);
    *data = result;
}

void get_cv (NNPDF::real  **data, int* n){
    int len = $self->GetNData();
    NNPDF::real * result = (decltype(result)) malloc(sizeof(*result)*len);
    for (int i = 0; i < len; i++){
        result[i] = self->GetObsCV(i);
    }
    *data = result;
    *n = len;
}

void get_error (NNPDF::real  **data, int* n){
    int len = $self->GetNData();
    NNPDF::real * result = (decltype(result)) malloc(sizeof(*result)*len);
    for (int i = 0; i < len; i++){
        result[i] = self->GetObsError(i);
    }
    *data = result;
    *n = len;
}

%pythoncode{

def __len__(self):
    return self.GetNData()


def __iter__(self):
    for i in range(len(self)):
        yield self.GetObsCV(i), self.GetObsError(i)

}
}

/*Positivity*/
%apply int *OUTPUT {int *res};
%apply (NNPDF::real** ARGOUTVIEWM_ARRAY2, int* DIM1, int* DIM2)\
{(NNPDF::real** result, int* ndata, int* npdf)}

%include "NNPDF/positivity.h"



%include "NNPDF/chisquared.h"
