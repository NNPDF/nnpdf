// $Id$
//
// NNPDF++ 2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#include "NNPDF/common.h"
#include "NNPDF/nnmpi.h"
#include "NNPDF/thpredictions.h"
using namespace NNPDF;

#ifdef OPENMPI
#include <mpi/mpi.h>
#ifdef SSE_CONV
#define NNMPI_REAL MPI_FLOAT
#else
#define NNMPI_REAL MPI_DOUBLE
#endif
#endif

#define UNUSED(expr) (void)(expr)

namespace NNPDF {

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
  }
#else
  // Basic convolution
  static inline void convolute(const real* __restrict__ pdf,const real* __restrict__ sig,real& retval,int const& n)
  {
    for (int i = 0; i < n; i++)
      retval += pdf[i]*sig[i];
  }
#endif

  MPI *MPI::mpiInstance = NULL;

  MPI::MPI():
    fnumbertasks(1),
    ftaskid(0),
    fworld(1)
  {
  }

  MPI::~MPI()
  {
  }

  void MPI::Finalize()
  {
#ifdef OPENMPI
    if (getInstance()->TaskID() == 0)
      {
        MPI::ResetWorld();
        int v = -1;
        for (int d = 1; d < getInstance()->NumTasks(); d++)
            MPI_Send(&v, 1, MPI_INT, d, 0, MPI_COMM_WORLD);
      }
    MPI_Finalize();
#endif
  }

  void MPI::Init()
  {
#ifdef OPENMPI
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD,&getInstance()->fworld);
    MPI_Comm_rank(MPI_COMM_WORLD,&getInstance()->ftaskid);
    getInstance()->fnumbertasks = getInstance()->fworld;
    if (getInstance()->ftaskid > 0) getInstance()->SlaveConvoluteLoop();
#endif
  }

  void MPI::DecomposeDomain(int domain_size, int taskid, int *subdomain_start, int *subdomain_size)
  {
#ifdef OPENMPI
    if (taskid == 0)
      if (getInstance()->fnumbertasks > domain_size)
        getInstance()->fnumbertasks = domain_size;

    *subdomain_start = domain_size / getInstance()->fnumbertasks * taskid;
    *subdomain_size = domain_size / getInstance()->fnumbertasks;
    if (taskid == getInstance()->fnumbertasks - 1)
      *subdomain_size += domain_size % getInstance()->fnumbertasks;
#else
  UNUSED(domain_size);
  UNUSED(taskid);
  UNUSED(subdomain_start);
  UNUSED(subdomain_size);
#endif
  }

  void MPI::SlaveConvoluteLoop()
  {
#ifdef OPENMPI
    MPI_Status status;
    int offset, chunksize, Npdf, Dsz, NData;
    real *theory;
    real *pdf;
    real *sigma;

    while(true)
      {
        // Perform the receiver
        MPI_Recv(&offset, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        if (offset == -1) break;

        MPI_Recv(&chunksize, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&Npdf, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
        MPI_Recv(&Dsz, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, &status);
        MPI_Recv(&NData, 1, MPI_INT, 0, 4, MPI_COMM_WORLD, &status);

        theory = new real[NData*Npdf];
        MPI_Recv(&theory[offset*Npdf], chunksize*Npdf, NNMPI_REAL, 0, 5, MPI_COMM_WORLD, &status);

        pdf = new real[Dsz*Npdf];
        MPI_Recv(pdf, Dsz*Npdf, NNMPI_REAL, 0, 6, MPI_COMM_WORLD, &status);

        sigma = new real[Dsz*NData];
        MPI_Recv(sigma, Dsz*NData, NNMPI_REAL, 0, 7, MPI_COMM_WORLD, &status);

        for (int i = offset; i < (offset+chunksize); i++)
          for (int n = 0; n < Npdf; n++)
          {
            theory[i*Npdf + n] = 0;
            convolute(pdf+Dsz*n,sigma+Dsz*i,theory[i*Npdf + n],Dsz);
          }

        MPI_Send(&theory[offset*Npdf], chunksize*Npdf, NNMPI_REAL, 0, 0, MPI_COMM_WORLD);

        delete[] theory;
        delete[] pdf;
        delete[] sigma;
      }
#endif
  }

  void MPI::SendJobs(int NData, int Npdf, int Dsz, real *theory, real *pdf, real *sigma)
  {
#ifdef OPENMPI
    int offset, chunksize;
    for (int d = 1; d < getInstance()->NumTasks(); d++)
      {
        MPI::DecomposeDomain(NData, d, &offset, &chunksize);
        MPI_Send(&offset, 1, MPI_INT, d, 0, MPI_COMM_WORLD);
        MPI_Send(&chunksize, 1, MPI_INT, d, 1, MPI_COMM_WORLD);

        MPI_Send(&Npdf, 1, MPI_INT, d, 2, MPI_COMM_WORLD);
        MPI_Send(&Dsz, 1, MPI_INT, d, 3, MPI_COMM_WORLD);
        MPI_Send(&NData, 1, MPI_INT, d, 4, MPI_COMM_WORLD);

        MPI_Send(&theory[offset*Npdf],chunksize*Npdf,NNMPI_REAL, d, 5, MPI_COMM_WORLD);

        MPI_Send(pdf, Dsz*Npdf, NNMPI_REAL, d, 6, MPI_COMM_WORLD);
        MPI_Send(sigma, Dsz*NData, NNMPI_REAL, d, 7,MPI_COMM_WORLD);

        //printf("Sent %d elements to task %d offset= %d\n",chunksize,d,offset);
      }
#else
  UNUSED(NData);
  UNUSED(Npdf);
  UNUSED(Dsz);
  UNUSED(theory);
  UNUSED(pdf);
  UNUSED(sigma);
#endif
  }

  void MPI::RecvJobs(int NData, int Npdf, real *theory)
  {
#ifdef OPENMPI
    MPI_Status status;
    int offset, chunksize;
    for (int d = 1; d < getInstance()->NumTasks(); d++)
      {
        MPI::DecomposeDomain(NData, d, &offset, &chunksize);
        MPI_Recv(&theory[offset*Npdf], chunksize*Npdf, NNMPI_REAL, d, 0, MPI_COMM_WORLD, &status);
      }
#else
  UNUSED(NData);
  UNUSED(Npdf);
  UNUSED(theory);
#endif
  }

  void MPI::ResetWorld()
  {
#ifdef OPENMPI
    getInstance()->fnumbertasks = getInstance()->fworld;
#endif
  }

}
