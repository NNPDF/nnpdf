//
// The NNPDF++ Collaboration - 2012
// Author: Stefano Carrazza, stefano.carrazza@mi.infn.it
//

#include <iostream>
#include <stdlib.h>
#include <cmath>

#include "evolqed.h"
#include "gsl/gsl_sf_psi.h"

//! Quark charges (u,d,s,c,b,t)
const double eu =  2.0/3.0;
const double ed = -1.0/3.0;
const double e[] = {eu, ed, ed, eu, ed, eu};
const double enns[] = { eu, ed, ed, eu, eu, ed, ed, eu, ed, eu};

//! Constant for anomalous dimension calculation
const double emc = 0.5772156649;

//! LAPACK Fortran wrapper for zgeev
extern "C" { 
  void zgeev_(char *jobvl, char *jobvr, int *n, dcomplex *a,
	      int *lda, dcomplex *w, dcomplex *vl,
	      int *ldvl, dcomplex *vr, int *ldvr,
	      dcomplex *work, int *lwork, double *rwork,
	      int *info);

  void zgetrf_(int *M, int *N, dcomplex* A, 
	       int* lda, int*IPIV, int* INFO);

  void zgetri_(int *N, dcomplex *A, int* lda, int* IPIV,
	       dcomplex *WORK, int* lwork, int* INFO);
}

/**
  * The constructor
  */
EvolQED::EvolQED(int nf, bool running, int algorithm):
  fNF(nf),
  fAlgorithm(algorithm),
  fNINT(50),
  fNTRUNC(10),
  fAlphaRunning(running),
  fQ20Ref(1.777*1.777), // tau mass
  fAlphaQ20Ref(1/133.4) // at tau mass
{
  if (nf < 3 || nf > 6)
    {
      cerr << "EvolQED: Error Nf outside available range." << endl;
      exit(-1);
    }
}

/**
  * Cleanup variables and destroys the class
  */
EvolQED::~EvolQED()
{
}

/**
  * Computes the running of alpha QED at LO
  * (see equation 2.9 in the QED notes)
  * \param Q2 the current scale to compute
  * \return the running alpha QED at Q2
  */
double EvolQED::AlphaRunning(double Q2)
{
  return fAlphaQ20Ref/(1 - Beta0()*fAlphaQ20Ref*log(Q2/fQ20Ref));
}

/**
  * Computes the Beta0 coefficient in QED
  * (see equation 2.10 in the QED doc)
  */
double EvolQED::Beta0()
{
  double b0 = 0;

  b0 += 3; // e,mu,tau charges

  for (int nf = 0; nf < fNF; nf++)
    b0 += 3.0*pow(e[nf], 2.0); // 3 colors for each quark

  return 2.0/(3.0*M_PI)*b0;
}

/**
  * Compute the Anomalous dimension Gqq at LO
  * \param N the Mellin variable
  * (see equation 3.23 in the QED doc)
  */
dcomplex EvolQED::getADqq(dcomplex N)
{
  dcomplex sum(0,0);
  sum = dcomplex(emc,0.0) + Digamma(N + 1.0);

  return 3.0/2.0 - 2.0*sum + 1.0/(N*(N + 1.0));
}

/**
  * Compute the Anomalous dimension Gfq at LO
  * \param N the Mellin variable
  * (see equation 3.21 in the QED doc)
  */
dcomplex EvolQED::getADfq(dcomplex N)
{
  return (N*N + N + 2.0)/(N*(N*N - 1.0));
}

/**
  * Compute the Anomalous dimension Gff at LO
  * \param N the Mellin variable
  * (see equation 3.22 in the QED doc)
  */
dcomplex EvolQED::getADff()
{
  // N not used at LO
  return -2.0/3.0;
}

/**
  * Compute the Anomalous dimension Gqf at LO
  * \param N the Mellin variable
  * (see equation 3.20 in the QED doc)
  */
dcomplex EvolQED::getADqf(dcomplex N)
{
  return (N*N + N + 2.0)/(N*(N + 1.0)*(N + 2.0));
}

/**
  * Computes the Evolution Factors in N space
  * \param N the Mellin variable
  * \param Q2I the initial Q2
  * \param Q2F the final Q2
  * \return EFNNS the non-singlet solution, 1D array dim(10)
  * \return EFNSG the singlet solution, 2D array dim(3,3)
  */
void EvolQED::EvolFactNQED(dcomplex N, double Q2I, double Q2F,
                           dcomplex* EFNNS, dcomplex** EFNSG)
{
  // Initializing variables
  double alphaQ2I = AlphaRunning(Q2I);
  double alphaQ2F = AlphaRunning(Q2F); 

  // Anomalous dimensions
  dcomplex Gqq = getADqq(N);
  dcomplex Gqf = getADqf(N);
  dcomplex Gfq = getADfq(N);
  dcomplex Gff = getADff();

  ///////////////////////////////////////
  // Building the Non-singlet solution
  // equation 4.6 from QED notes, see special QED basis

  for (int n = 0; n < 10; n++)
    {
      EFNNS[n] = dcomplex(0.0, 0.0);

      // Filtering conditions
      if (fNF == 3 && (n == 0 || (n > 1 && n < 4) || n > 6)) continue;
      if (fNF == 4 && ( (n > 1 && n < 4) || n > 7) ) continue;
      if (fNF == 5 && ( n == 3 || n > 8) ) continue;

      const double e2 = enns[n]*enns[n];

      if (fAlphaRunning == true)
        EFNNS[n] = exp( e2*Gqq*log(alphaQ2F/alphaQ2I)/(2.0*M_PI*Beta0()) );
      else
        EFNNS[n] = exp( e2*Gqq*fAlphaQ20Ref*log(Q2F/fQ20Ref)/(2.0*M_PI) );
    }

  ////////////////////////////////
  // Building the Singlet solution
  // see equation 4.2

  int Nc = 3; // Quark color multiplicity
  double nfup = 3;
  double nfdn = 3;

  if (fNF == 3) { nfup = 1; nfdn = 2; } // u,d,s
  if (fNF == 4) { nfup = 2; nfdn = 2; } // u,d,s,c
  if (fNF == 5) { nfup = 2; nfdn = 3; } // u,d,s,c,b

  double delta  = (nfup - nfdn)/(double)fNF;
  double etap   = 0.5*(eu*eu + ed*ed);
  double etam   = 0.5*(eu*eu - ed*ed);
  double thetap = 2.0*Nc*fNF*(delta*etap + etam);
  double thetam = 2.0*Nc*fNF*(etap + delta*etam);

  // e_sigma
  double esgm2 = nfup*eu*eu + nfdn*ed*ed;

  // Initializing Matrices and Vectors for the diagonalization
  dcomplex **A    = new dcomplex*[3];

  for (int i = 0; i < 3; i++)
    A[i] = new dcomplex[3];

  // Building DGLAP equation 4.2
  A[0][0] = esgm2*Gff;  A[0][1] = etap*Gfq; A[0][2] = etam*Gfq;
  A[1][0] = thetam*Gqf; A[1][1] = etap*Gqq; A[1][2] = etam*Gqq;
  A[2][0] = thetap*Gqf; A[2][1] = etam*Gqq; A[2][2] = etap*Gqq;

  // Solving using numerical diagonalization
  if (fAlgorithm == 0)
    {
      cout << " - DGLAP SG solution by numerical diagonalization" << endl;

      dcomplex *D     = new dcomplex[3];
      dcomplex **P    = new dcomplex*[3];
      dcomplex **PInv = new dcomplex*[3];

      for (int i = 0; i < 3; i++)
        {
          D[i] = dcomplex(0.0, 0.0);
          P[i] = new dcomplex[3];
          PInv[i] = new dcomplex[3];
          for (int j = 0; j < 3; j++)
            P[i][j] = PInv[i][j] = dcomplex(0.0, 0.0);
        }

      // Getting eigenvalues and eigenvectors
      ComputeEigenValuesVectors(3, A, D, P, PInv);

      // Solution in the diagonal basis
      for (int n = 0; n < 3; n++)
        {
          if (fAlphaRunning == true)
            D[n] = exp( D[n]*log(alphaQ2F/alphaQ2I)/(2.0*M_PI*Beta0()) );
          else
            D[n] = exp( D[n]*fAlphaQ20Ref*log(Q2F/fQ20Ref)/(2.0*M_PI) );
        }

      // Final Rotation P*D*PInv
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          {
            dcomplex sum(0.0, 0.0);
            for (int k = 0; k < 3; k++)
              sum += P[i][k]*D[k]*PInv[k][j];

            EFNSG[i][j] = sum;
          }

      for (int i = 0; i < 3; i++)
        {
          if (P[i]) delete[] P[i];
          if (PInv[i]) delete[] PInv[i];
        }

      delete[] P;
      delete[] PInv;
      delete[] D;

    }
  else if (fAlgorithm == 1)
    {
      cout << " - DGLAP SG solution by path ordering" << endl;

      // Path-ordering solution
      const int nint = fNINT;
      const int ntrunc = fNTRUNC;
      double da = (alphaQ2F - alphaQ2I) / (double) nint;

      dcomplex **PSG = new dcomplex*[3];
      dcomplex **SPSG = new dcomplex*[3];
      dcomplex **ONE = new dcomplex*[3];
      dcomplex ***C  = new dcomplex**[ntrunc];

      for (int i = 0; i < 3; i++)
        {
          PSG[i]  = new dcomplex[3];
          SPSG[i] = new dcomplex[3];
          ONE[i]  = new dcomplex[3];
        }

      for (int i = 0; i < ntrunc; i++)
        {
          C[i] = new dcomplex*[3];
          for (int j = 0; j < 3; j++)
            C[i][j] = new dcomplex[3];
        }

      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          {
            PSG[i][j] = A[i][j] * (1.0/alphaQ2F + 1.0/alphaQ2I)/ (2.0*M_PI*Beta0());
            SPSG[i][j] = da * PSG[i][j] / 2.0;
          }

      double AK = alphaQ2I;

      for (int k = 1; k < nint; k++)
        {
          AK += da;
          for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
              {
                PSG[i][j] = A[i][j] / (2.0*M_PI*Beta0()) / AK;
                SPSG[i][j] = SPSG[i][j] + da * PSG[i][j];
              }
        }

      double as_check = fabs(AK+da-alphaQ2F);
      if (as_check > 1e-7)
        {
          cerr << "In src/modules/evol" << endl;
          exit(-1);
        }

      // Now we exponentiate the matrix 3x3 SPSG[i][j]
      // I finally went for the expansion:
      //
      // exp(A) = 1 + C + C^2/2 + C^3/6

      for (int i = 0; i < 3; i++)
        {
          ONE[i][i] = dcomplex(1.0, 0.0);
          for (int j = 0; j < 3; j++)
            {
              for (int l = 0; l < ntrunc; l++)
                C[l][i][j] = dcomplex(0.0,0.0);
            }
        }

      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          {
            C[0][i][j] = SPSG[i][j];
            EFNSG[i][j] = ONE[i][j] + C[0][i][j];
          }

      for (int l = 1; l < ntrunc; l++)
        for (int i = 0; i < 3; i++)
           for (int j = 0; j < 3; j++)
            {
              for (int k = 0; k < 3; k++)
                C[l][i][j] += C[l-1][i][k]*SPSG[k][j] / (double) (l+1);
              EFNSG[i][j] += C[l][i][j];
            }

      for (int i = 0; i < 3; i++)
        {
          if (PSG[i]) delete[] PSG[i];
          if (SPSG[i]) delete[] SPSG[i];
          if (ONE[i]) delete[] ONE[i];
        }

      for (int l = 0; l < ntrunc; l++)
        if (C[l]) delete[] C[l];

      delete[] C;
      delete[] PSG;
      delete[] SPSG;
      delete[] ONE;

    }
  else
    {
      cerr << "Error algorithm not valid." << endl;
      exit(-1);
    }

  for (int i = 0; i < 3; i++)
    if (A[i]) delete[] A[i];

  delete[] A;
}

/**
  * Computes the eigenvalues of matrix A
  * using the numerical method zgeev from Lapack
  */
void EvolQED::ComputeEigenValuesVectors(int n, dcomplex **A, dcomplex *D, dcomplex **V, dcomplex **VInv)
{
  // Diagonalization by LAPACK zgeev
  // zgeev, complex algorithm for diagonalization
  char jobvl, jobvr;
  int lda, ldvl, ldvr, lwork, info;
  double *rwork;
  dcomplex *a, *vl, *vr, *work;

  jobvl = 'N'; // Don't calculate the left eigenvectors of A
  jobvr = 'V'; // Calculate the right eigenvectors of A

  lda = n; // The leading dimension of the matrix A

  // zgeev_ctof convert the matrix A from double pointer C
  // to single pointer Fortran form.
  a = new dcomplex[n*lda];
  for (int i = 0; i < n; i++)
    for (int j = 0; j < lda; j++)
      a[i + j*lda] = A[i][j];

  // we need to define the matrices for the eigenvectors
  ldvl = n;
  vl = new dcomplex[ldvl*n];
  ldvr = n;
  vr = new dcomplex[ldvr*n];
  lwork = 4*n;
  work = new dcomplex[lwork];
  rwork = new double[2*n];

  zgeev_(&jobvl, &jobvr, &n, a, &lda, D, vl, &ldvl, vr,
         &ldvr, work, &lwork, rwork, &info);

  dcomplex *vinv = new dcomplex[n*n];
  // Now we get the right eigenvectors into the output array
  for (int i = 0; i < ldvr; i++)
    for (int j = 0; j < n; j++)
      V[i][j] = vinv[i + j*n] = vr[i + j*n];
  
  // Inverting the V matrix
  int *ipiv = new int[n+1];
  int lwork2 = n*n, info2;
  dcomplex *work2 = new dcomplex[lwork2];

  zgetrf_(&n, &n, vinv, &n, ipiv, &info2);
  zgetri_(&n, vinv, &n, ipiv, work2, &lwork2, &info2);
  
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      VInv[i][j] = vinv[i + j*n];
  
  // Clean up the memory
  delete[] a;
  delete[] vl;
  delete[] vr;
  delete[] work;
  delete[] rwork;

  delete[] ipiv;
  delete[] work2;
  delete[] vinv;
  
}

/**
  * Computes the Digamma = PSI function
  * Using the gsl implementation, see gsl_sf_complex_psi_e()
  */
dcomplex EvolQED::Digamma(dcomplex x)
{
  gsl_sf_result r, i;
  int status = gsl_sf_complex_psi_e(x.real(), x.imag(), &r, &i);

  if (status != 0)
    {
      cerr << " Error in gsl PSI function." << endl;
      exit(-1);
    }

  return dcomplex(r.val, i.val);
}

/* C wrapper interfaces to C++ routine */
EvolQED *EvolQED__new(int nf, bool run)
{
  return new EvolQED(nf, run);
}

void EvolQED__EvolFactNQED(EvolQED *This, dcomplex N, double Q2I, double Q2F,
                           dcomplex* EFNNS, dcomplex** EFNSG)
{
  This->EvolFactNQED(N, Q2I, Q2F, EFNNS, EFNSG);
}

void EvolQED__delete(EvolQED *This)
{
  delete This;
}

void EvolQED__SetNF(EvolQED *This, int nf)
{
  This->SetNF(nf);
}

void EvolQED__ActivateRunning(EvolQED *This, bool run)
{
  This->ActivateRunning(run);
}

void EvolQED__SetRefCoupling(EvolQED *This, double Q20, double AlphaQ20)
{
  This->SetRefCoupling(Q20, AlphaQ20);
}

int EvolQED__GetNF(EvolQED *This)
{
  return This->GetNF();
}

double EvolQED__GetQ20Ref(EvolQED *This)
{
  return This->GetQ20Ref();
}

double EvolQED__GetAlphaQ20(EvolQED *This)
{
  return This->GetAlphaQ20();
}


