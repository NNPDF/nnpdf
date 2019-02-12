// $Id: thpredictions.cc 3184 2015-08-21 19:54:00Z stefano.carrazza@mi.infn.it $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cstring>

#include "NNPDF/common.h"
#include "NNPDF/experiments.h"
#include "NNPDF/thpredictions.h"
#include "NNPDF/dataset.h"
#include "NNPDF/fkset.h"
#include "NNPDF/fastkernel.h"
#include "NNPDF/pdfset.h"
#include "NNPDF/utils.h"
#include "NNPDF/timer.h"
#include "NNPDF/nnmpi.h"
#include "NNPDF/exceptions.h"
using namespace NNPDF;

// Default verbosity level
bool ThPredictions::Verbose = true;

#ifdef SSE_CONV
  #include <pmmintrin.h>

  static inline void convolute(const real* __restrict__ x, const real* __restrict__ y, real& retval, int const& n)
  {
    __m128 acc = _mm_setzero_ps();

    const __m128* a = (const __m128*) x;
    const __m128* b = (const __m128*) y;

    for (int i=0; i<n/4; i++)
      acc = _mm_add_ps(acc, _mm_mul_ps(*(a+i),*(b+i)));

    // Turns out, hadd is a waste of time, 2*shuff +2*add faster on ECDF
    const __m128 shuffle1 = _mm_shuffle_ps(acc, acc, _MM_SHUFFLE(1,0,3,2));
    const __m128 tmp1     = _mm_add_ps(acc, shuffle1);
    const __m128 shuffle2 = _mm_shuffle_ps(tmp1, tmp1, _MM_SHUFFLE(2,3,0,1));
    const __m128 tmp2     = _mm_add_ps(tmp1, shuffle2);

    _mm_store_ss(&retval,tmp2);
    _mm_empty();

    return;

    /*
    __m256 acc = _mm256_setzero_ps();

    const __m256* a = (const __m256*) x;
    const __m256* b = (const __m256*) y;

    for (int i=0; i<n/8; i++)
      acc = _mm256_add_ps(acc, _mm256_mul_ps(*(a+i),*(b+i)));

    __m256 t1 = _mm256_hadd_ps(acc,acc);
    __m256 t2 = _mm256_hadd_ps(t1,t1);
    __m128 t3 = _mm256_extractf128_ps(t2,1);
    __m128 t4 = _mm_add_ss(_mm256_castps256_ps128(t2),t3);
    retval = _mm_cvtss_f32(t4);

    _mm_empty();
    */


  }

#else

  // Basic convolution
  static inline void convolute(const real* __restrict__ pdf,const real* __restrict__ sig,real& retval,int const& n)
  {
    for (int i = 0; i < n; i++)
      retval += pdf[i]*sig[i];
  }

#endif

/**
 * Constructor for observables only
 */
ThPredictions::ThPredictions(const PDFSet *pdfset, const Experiment *exp):
fObs(NULL),
fNpdf(pdfset->GetMembers()),
fNData(exp->GetNData()),
fPDFName(pdfset->GetSetName()),
fSetName(exp->GetExpName()),
fEtype(pdfset->GetEtype())
{
  // New theory array
  fObs = new real[fNData*fNpdf];

  Timer timer=Timer();
  timer.start();

  int index = 0;
  for (int i=0; i<exp->GetNSet(); i++)
  {
    ThPredictions::Convolute(pdfset,&(exp->GetSet(i)),fObs+index);
    index+=pdfset->GetMembers()*exp->GetSet(i).GetNData();
  }

  fTconv=timer.stop();
}

/**
 * Constructor for observables only
 */
ThPredictions::ThPredictions(const PDFSet *pdfset, const FKTable *fktab):
fObs(NULL),
fNpdf(pdfset->GetMembers()),
fNData(fktab->GetNData()),
fPDFName(pdfset->GetSetName()),
fSetName(fktab->GetDataName()),
fEtype(pdfset->GetEtype())
{
  // New theory array
  fObs = new real[fNData*fNpdf];

  Timer timer=Timer();
  timer.start();

  // Calculate observables
  ThPredictions::Convolute(pdfset,fktab,fObs);

  fTconv=timer.stop();
}

/**
 * Constructor for observables only
 */
ThPredictions::ThPredictions(const PDFSet *pdfset, const FKSet *fkset):
fObs(NULL),
fNpdf(pdfset->GetMembers()),
fNData(fkset->GetNDataFK()),
fPDFName(pdfset->GetSetName()),
fSetName(fkset->GetDataName()),
fEtype(pdfset->GetEtype())
{
  // New theory array
  fObs = new real[fNData*fNpdf];

  Timer timer=Timer();
  timer.start();

  // Calculate observables
  ThPredictions::Convolute(pdfset,fkset,fObs);

  fTconv=timer.stop();
}

/**
 * Constructor for predictions with different beam PDFs
 */
ThPredictions::ThPredictions(const PDFSet *pdf1, const PDFSet *pdf2, const FKTable* fktab):
fObs(NULL),
fNpdf(pdf1->GetMembers()),
fNData(fktab->GetNData()),
fPDFName(pdf1->GetSetName()+"_"+pdf2->GetSetName()),
fSetName(fktab->GetDataName()),
fEtype(pdf1->GetEtype())
{

  if (pdf1->GetMembers() != pdf2->GetMembers())
    throw RangeError("ThPredictions::ThPredictions","PDFs for different beams must have the same number of members!");

  if (pdf1->GetEtype() != pdf2->GetEtype())
    throw RangeError("ThPredictions::ThPredictions","PDFs for different beams must have the same error type!");

  // New theory array
  fObs = new real[fNData*fNpdf];

  Timer timer=Timer();
  timer.start();

  // Calculate observables
  ThPredictions::Convolute(pdf1,pdf2,fktab,fObs);

  fTconv=timer.stop();
}


/**
 * Copy Constructor
 */
ThPredictions::ThPredictions(const ThPredictions& o):
fObs(new real[o.fNData*o.fNpdf]),
fTconv(o.fTconv),
fNpdf(o.fNpdf),
fNData(o.fNData),
fPDFName(o.fPDFName),
fSetName(o.fSetName),
fEtype(o.fEtype)
{
// Copy predictions
 for (int i=0; i<fNpdf*fNData; i++)
  fObs[i] = o.fObs[i];
}

/**
 * Empty Constructor
 */
ThPredictions::ThPredictions(std::string pdfname, std::string setname, int nPDF, int nDat, PDFSet::erType erty):
fObs(new real[nPDF*nDat]()),
fTconv(0),
fNpdf(nPDF),
fNData(nDat),
fPDFName(pdfname),
fSetName(setname),
fEtype(erty)
{
}


void NNPDF::swap(ThPredictions& lhs, ThPredictions& rhs)
{
  using std::swap;
  swap(lhs.fObs, rhs.fObs);
  swap(lhs.fTconv, rhs.fTconv);
  swap(lhs.fNpdf, rhs.fNpdf);
  swap(lhs.fNData, rhs.fNData);
  swap(lhs.fPDFName, rhs.fPDFName);
  swap(lhs.fSetName, rhs.fSetName);
  swap(lhs.fEtype, rhs.fEtype);
}

ThPredictions& ThPredictions::operator=(ThPredictions other){
  using std::swap;
  swap(*this, other);
  return *this;
}

ThPredictions::ThPredictions(ThPredictions && other):
fObs(nullptr),
fTconv(0),
fNpdf(0),
fNData(0),
fPDFName(""),
fSetName(""),
fEtype(PDFSet::erType::ER_NONE)
{
  swap(*this, other);
}

/**
 * The destructor
 */
ThPredictions::~ThPredictions()
{
  if (fObs) delete[] fObs;
}


ThPredictions& ThPredictions::operator+=(const ThPredictions& o)
{
  if (o.fEtype != fEtype)
    throw EvaluationError("ThPredictions::operator+=","Cannot sum predictions with different error values");

  if (o.fNpdf != fNpdf || o.fNData != fNData)
    throw EvaluationError("ThPredictions::operator+=","Cannot sum predictions with different numbers of datapoints/PDF members");

  if (fEtype != PDFSet::erType::ER_MC && fEtype != PDFSet::erType::ER_MC68 && fEtype != PDFSet::erType::ER_MCT0 && ThPredictions::Verbose == true)
    std::cerr << "ThPredictions::operator+= Warning: Observable summation undefined for ErrorTypes other than Monte-Carlo" <<std::endl;

  // Increment predictions
  for (int i=0; i<fNpdf*fNData; i++)
    fObs[i] += o.fObs[i];
  
  fTconv += o.fTconv;
  if (fPDFName.compare(o.fPDFName) != 0)
    {std::stringstream pdf; pdf << "(" << fPDFName <<"+"<<o.fPDFName << ")"; fPDFName = pdf.str();}
  if (fSetName.compare(o.fSetName) != 0)
    {std::stringstream set; set << "(" << fSetName <<"+"<<o.fSetName << ")"; fSetName = set.str();}
  
  return *this;
}

ThPredictions& ThPredictions::operator-=(const ThPredictions& o)
{
  if (o.fEtype != fEtype)
    throw EvaluationError("ThPredictions::operator-=","Cannot subtract predictions with different error values");

  if (o.fNpdf != fNpdf || o.fNData != fNData)
    throw EvaluationError("ThPredictions::operator-=","Cannot subtract predictions with different numbers of datapoints/PDF members");

  if ( ( fEtype != PDFSet::erType::ER_MC && fEtype != PDFSet::erType::ER_MC68 && fEtype != PDFSet::erType::ER_MCT0 ) && ThPredictions::Verbose == true)
    std::cerr << "ThPredictions::operator-= Warning: Observable subtraction undefined for ErrorTypes other than Monte-Carlo" <<std::endl;

  // Increment predictions
  for (int i=0; i<fNpdf*fNData; i++)
    fObs[i] -= o.fObs[i];

  fTconv += o.fTconv;
  if (fPDFName.compare(o.fPDFName) != 0)
    {std::stringstream pdf; pdf << "(" << fPDFName <<"-"<<o.fPDFName << ")"; fPDFName = pdf.str();}
  if (fSetName.compare(o.fSetName) != 0)
    {std::stringstream set; set << "(" << fSetName <<"-"<<o.fSetName << ")"; fSetName = set.str();}

  return *this;
}

ThPredictions& ThPredictions::operator/=(const ThPredictions& o)
{
  if (o.fEtype != fEtype)
    throw EvaluationError("ThPredictions::operator/=","Cannot divide predictions with different error values");

  if (o.fNpdf != fNpdf || o.fNData != fNData)
    throw EvaluationError("ThPredictions::operator/=","Cannot divide predictions with different numbers of datapoints/PDF members");

  if ( (fEtype != PDFSet::erType::ER_MC && fEtype != PDFSet::erType::ER_MC68 && fEtype != PDFSet::erType::ER_MCT0 ) && ThPredictions::Verbose == true)
    std::cerr << "ThPredictions::operator/= Warning: Observable division undefined for ErrorTypes other than Monte-Carlo" <<std::endl;

  // Increment predictions
  for (int i=0; i<fNpdf*fNData; i++)
    fObs[i] /= o.fObs[i];

  fTconv += o.fTconv;
  if (fPDFName.compare(o.fPDFName) != 0)
    {std::stringstream pdf; pdf << "(" << fPDFName <<"/"<<o.fPDFName << ")"; fPDFName = pdf.str();}
  if (fSetName.compare(o.fSetName) != 0)
    {std::stringstream set; set << "(" << fSetName <<"/"<<o.fSetName << ")"; fSetName = set.str();}

  return *this;
}

ThPredictions& ThPredictions::operator*=(const ThPredictions& o)
{
  if (o.fEtype != fEtype)
    throw EvaluationError("ThPredictions::operator/=","Cannot multiply predictions with different error values");

  if (o.fNpdf != fNpdf || o.fNData != fNData)
    throw EvaluationError("ThPredictions::operator/=","Cannot multiply predictions with different numbers of datapoints/PDF members");

  if ( (fEtype != PDFSet::erType::ER_MC && fEtype != PDFSet::erType::ER_MC68 && fEtype != PDFSet::erType::ER_MCT0 ) && ThPredictions::Verbose == true)
    std::cerr << "ThPredictions::operator/= Warning: Observable multiplication undefined for ErrorTypes other than Monte-Carlo" <<std::endl;

  // Increment predictions
  for (int i=0; i<fNpdf*fNData; i++)
    fObs[i] *= o.fObs[i];

  fTconv += o.fTconv;
  if (fPDFName.compare(o.fPDFName) != 0)
    {std::stringstream pdf; pdf << "(" << fPDFName <<"*"<<o.fPDFName << ")"; fPDFName = pdf.str();}
  if (fSetName.compare(o.fSetName) != 0)
    {std::stringstream set; set << "(" << fSetName <<"*"<<o.fSetName << ")"; fSetName = set.str();}

  return *this;
}


const ThPredictions ThPredictions::operator+(const ThPredictions &o) const
{  ThPredictions res = *this;  res += o; return res;   }

const ThPredictions ThPredictions::operator-(const ThPredictions &o) const
{  ThPredictions res = *this;  res -= o; return res;   }

const ThPredictions ThPredictions::operator/(const ThPredictions &o) const
{  ThPredictions res = *this;  res /= o; return res;   }

const ThPredictions ThPredictions::operator*(const ThPredictions &o) const
{  ThPredictions res = *this;  res *= o; return res;   }

  /**
 * Static FK level convolution function
 */
void ThPredictions::Convolute(const PDFSet *pdfset, const FKTable *fk, real* theory)
{
  const int Npdf   = pdfset->GetMembers();
  const int NData  = fk->GetNData();
  const int Dsz    = fk->GetDSz();
  const int Psz    = sizeof(real)*Dsz*Npdf;

  real* sigma = fk->GetSigma();
  real *pdf = 0;
  int err = posix_memalign(reinterpret_cast<void **>(&pdf), 16, Psz);
  std::memset(pdf,0,Psz);
  if (err != 0)
    throw RangeError("ThPredictions::Convolute","ThPredictions posix_memalign " + std::to_string(err));

  // Fetch PDF array
  GetNZPDF(pdfset, fk, pdf);

  // Switcher between standard sequencial (+optional openmp) convolution
  // or openmpi point-to-point communication implementation
#ifndef OPENMPI
  // Calculate observables
  #pragma omp parallel for
  for (int i = 0; i < NData; i++)
    for (int n = 0; n < Npdf; n++)
    {
      theory[i*Npdf + n] = 0;
      convolute(pdf+Dsz*n,sigma+Dsz*i,theory[i*Npdf + n],Dsz);
    }
#else
  int offset, chunksize;
  MPI::DecomposeDomain(NData, 0, &offset, &chunksize); // compule data for master
  MPI::SendJobs(NData, Npdf, Dsz, theory, pdf, sigma); // prepare slaves
  for (int i = 0; i < chunksize; i++)
    for (int n = 0; n < Npdf; n++)
    {
      theory[i*Npdf + n] = 0;
      convolute(pdf+Dsz*n,sigma+Dsz*i,theory[i*Npdf + n],Dsz);
    }
  MPI::RecvJobs(NData, Npdf, theory);
#endif

  // Delete pdfs
  free(reinterpret_cast<void *>(pdf));
  return;
}

  /**
 * Static FK level convolution function - differing beams - NO CHECK ON MATCHING HERE
 */
void ThPredictions::Convolute(const PDFSet *pdf1, const PDFSet *pdf2, const FKTable *fk, real* theory)
{
  const int Npdf   = pdf1->GetMembers();
  const int NData  = fk->GetNData();
  const int Dsz    = fk->GetDSz();
  const int Psz    = sizeof(real)*Dsz*Npdf;

  real* sigma = fk->GetSigma();
  real *pdf = 0;
  int err = posix_memalign(reinterpret_cast<void **>(&pdf), 16, Psz);
  std::memset(pdf,0,Psz);

  if (err != 0)
    throw RangeError("ThPredictions::Convolute","ThPredictions posix_memalign " + std::to_string(err));

  // Fetch PDF array
  GetNZPDF(pdf1, pdf2, fk, pdf);

  // Switcher between standard sequencial (+optional openmp) convolution
  // or openmpi point-to-point communication implementation
#ifndef OPENMPI
  // Calculate observables
  #pragma omp parallel for
  for (int i = 0; i < NData; i++)
    for (int n = 0; n < Npdf; n++)
    {
      theory[i*Npdf + n] = 0;
      convolute(pdf+Dsz*n,sigma+Dsz*i,theory[i*Npdf + n],Dsz);
    }
#else
  int offset, chunksize;
  MPI::ResetWorld();
  MPI::DecomposeDomain(NData, 0, &offset, &chunksize); // compule data for master
  MPI::SendJobs(NData, Npdf, Dsz, theory, pdf, sigma); // prepare slaves
  for (int i = 0; i < chunksize; i++)
    for (int n = 0; n < Npdf; n++)
    {
      theory[i*Npdf + n] = 0;
      convolute(pdf+Dsz*n,sigma+Dsz*i,theory[i*Npdf + n],Dsz);
    }
  MPI::RecvJobs(NData, Npdf, theory);
#endif

  // Delete pdfs
  free(reinterpret_cast<void *>(pdf));
  return;
}

/**
 * Static convolution function - required for a fast fit
 */
void ThPredictions::Convolute(const PDFSet *pdfset, const FKSet *fkset, real* theory)
{
  const int NSigma = fkset->GetNSigma();

  if (NSigma == 1)
    ThPredictions::Convolute(pdfset,fkset->GetFK(0), theory );
  else
  {
    SigmaOp   op = fkset->GetOperator();
    const int NData = fkset->GetNDataFK();
    const int Npdf  = pdfset->GetMembers();

    std::vector<real*> results;
    for (int i=0; i<NSigma; i++)
      results.push_back(new real[NData*Npdf]);

    for (int i=0; i<NSigma; i++)
      ThPredictions::Convolute(pdfset, fkset->GetFK(i), results[i] );

    op( NData*Npdf, results, theory );

    for (int i=0; i<NSigma; i++)
      delete[] results[i];

    results.clear();
  }

  return;
}

/**
 * Compute the central value of observables
 * for PDFs from MC or Hessian approach
 * \param idat the data point [0,N)
 * \return the mean value from all replicas
 */
real ThPredictions::GetObsCV(int const& idat) const
{
  real avg = 0;

  if (fEtype == PDFSet::erType::ER_MC || fEtype == PDFSet::erType::ER_MC68)
    {
      real *iObs = fObs + idat*fNpdf;
      avg = ComputeAVG(fNpdf, iObs);
    }
  else
    {
      avg = fObs[idat*fNpdf];
    }

  return avg;
}

/**
 * Compute the error of observables central values
 * for PDFs from MC or Hessian approach
 * \param idat the data point [0,N)
 * \return the error value from all replicas
 */
real ThPredictions::GetObsError(int const& idat) const
{
  real err = 0;
  if (fEtype == PDFSet::erType::ER_MC || fEtype == PDFSet::erType::ER_MC68)
    {
      real *iObs = fObs + idat*fNpdf;
      err = ComputeStdDev(fNpdf, iObs);
    }
  else
    {
      real *iObs = fObs + idat*fNpdf;
      err = (fEtype == PDFSet::erType::ER_SYMEIG) ? ComputeSymEigErr(fNpdf, iObs) : ComputeEigErr(fNpdf, iObs);
      if (fEtype == PDFSet::erType::ER_EIG90)
        err /= 1.64485;
    }

  return err;
}

/**
 * Export ThPredictions observables to ostream
 */
void ThPredictions::Print(std::ostream& out, bool latex) const
{
  if (latex)
  {
    // Remove underscores for Latex processing
    std::string pdfname = fPDFName;
    std::replace(pdfname.begin(), pdfname.end(), '_', ' ');

    out << "\\begin{table}[htdp]" << std::endl
        << "\\caption{default}" << std::endl
        << "\\begin{center}" << std::endl
        << "\\begin{tabular}{|c|c|}" << std::endl
        << "\\hline" << std::endl
        << "\\multicolumn{2}{|c|}{"<<fSetName<<"} \\\\" << std::endl
        << "\\hline" << std::endl
        << "Datapoint & "<<pdfname<<" \\\\" << std::endl
        << "\\hline" << std::endl;

    for (int i=0; i<fNData; i++)
      out <<i+1<<" & $"<< GetObsCV(i) <<" \\pm "<<GetObsError(i)<<"$ \\\\"<< std::endl;

    out << "\\hline" << std::endl
        << "\\end{tabular}" << std::endl
        << "\\end{center}" << std::endl
        << "\\label{default}" << std::endl
        << "\\end{table}" << std::endl;
  }
  else
  {
    out << "THPREDICTIONS "<<fSetName<<"  "<<fPDFName<<std::endl;
    for (int i=0; i<fNData; i++)
      out <<std::left<<std::setw(5)<<i<<std::setw(20)<< GetObsCV(i) <<std::setw(20)<<GetObsError(i)<<std::endl;
  }
}


/*
 * For the xgrid, nonzero flavours in data, construct the pdf array.
 * \param data the dataset pointer
 * \param pdf the pdf array
 * \param grid the index for the FK grid
 * \return pdf the array
 */

void ThPredictions::GetNZPDF(const PDFSet* pdfset, const FKTable *fk, real* pdf)
{
  if (!fk)
    throw RangeError("ThPredictions::GetNZPDF","FKTable null");

  if (!pdfset)
    throw RangeError("ThPredictions::GetNZPDF","PDFSet null");

  // Get PDF set data
  const int NPDF     = pdfset->GetMembers();

  const double* xgrid  = fk->GetXGrid();
  const int* NZFlmap = fk->GetFlmap();
  const double Q20     = fk->GetQ20();

  const int Nfl      = 14;
  const int NonZero  = fk->GetNonZero();
  const int Nx       = fk->GetNx();
  const int Tx       = fk->GetTx();
  const int DSz      = fk->GetDSz();

  if ( fk->IsHadronic() )
  {
    // Hadronic process
    real* EVLN = new real[Nx*Nfl];

    for (int n = 0; n < NPDF; n++)
    {
      for (int i = 0; i < Nx; i++)
        pdfset->GetPDF(xgrid[i],Q20,n,&EVLN[i*Nfl]);

      for (int fl=0; fl<NonZero; fl++)
      {
        const int fl1 = NZFlmap[2*fl];
        const int fl2 = NZFlmap[2*fl+1];
        const int idx = n*DSz + fl*Tx;
        for (int i = 0; i < Nx; i++)
          for (int j = 0; j < Nx; j++)
            pdf[ idx + i*Nx + j ] = EVLN[i*Nfl+fl1]*EVLN[j*Nfl+fl2];
      }
    }

    delete[] EVLN;

  } else {
    // DIS
    real* EVLN = new real[Nfl];
    for ( int n = 0; n < NPDF; n++)
      for ( int i = 0; i < Nx; i++)
      {
        pdfset->GetPDF(xgrid[i],Q20,n,EVLN);
        for (int fl=0; fl<NonZero; fl++)
          pdf[ DSz*n + Nx*fl + i ] = EVLN[NZFlmap[fl]];
      }

    delete[] EVLN;
  }

}

/*
 * Build the PDF array with two different beam PDFs
 */

void ThPredictions::GetNZPDF(const PDFSet* pdf1, const PDFSet* pdf2, const FKTable *fk, real* pdf)
{
  if (!fk)
    throw RangeError("ThPredictions::GetNZPDF","FKTable null");

  if (!pdf1 || !pdf2)
    throw RangeError("ThPredictions::GetNZPDF","PDFSet null");

  if (!fk->IsHadronic())
    throw EvaluationError("ThPredictions::GetNZPDF","Dual-beam convolution requires an hadronic FK table!");

  // Get PDF set data
  const int NPDF     = pdf1->GetMembers();

  const double* xgrid  = fk->GetXGrid();
  const int* NZFlmap = fk->GetFlmap();
  const double Q20     = fk->GetQ20();

  const int Nfl      = 14;
  const int NonZero  = fk->GetNonZero();
  const int Nx       = fk->GetNx();
  const int DSz      = fk->GetDSz();

  // Hadronic process
  real* EVLN1 = new real[Nx*Nfl];
  real* EVLN2 = new real[Nx*Nfl];

  for (int n = 0; n < NPDF; n++)
  {
    for (int i = 0; i < Nx; i++)
    {
      pdf1->GetPDF(xgrid[i],Q20,n,&EVLN1[i*Nfl]);
      pdf2->GetPDF(xgrid[i],Q20,n,&EVLN2[i*Nfl]);
    }

    for (int fl=0; fl<NonZero; fl++)
    {
      const int fl1 = NZFlmap[2*fl];
      const int fl2 = NZFlmap[2*fl+1];

      for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Nx; j++)
          pdf[n*DSz + fl*Nx*Nx + i*Nx + j] = EVLN1[i*Nfl+fl1]*EVLN2[j*Nfl+fl2];
    }

  }

  delete[] EVLN1;
  delete[] EVLN2;

}

