//
// The NNPDF++ Collaboration - 2012
// Author: Stefano Carrazza, stefano.carrazza@mi.infn.it
//

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include "evolqed.h"

using namespace std;

/**
 * Testing the QED DGLAP solution
 */

int main()
{
  const int  nf = 3;
  const bool run = true;
  const double Q2I = 2.0;
  const double Q2F = 10000.0;

  cout << "Nf = " << nf << endl;
  cout << "Running alphaQED: " << run << endl;
  cout << endl;

  EvolQED *evol = new EvolQED(nf, run, 1);

  dcomplex *solns = new dcomplex[10];
  dcomplex **solsg = new dcomplex*[3];
  for (int i = 0; i < 3; i++)
    solsg[i] = new dcomplex[3];

  struct timespec start, finish;
  double elapsed;

  clock_gettime(CLOCK_MONOTONIC, &start);

  for (double N = 2.0; N < 3.0; N = N+1.0)
    {
      evol->EvolFactNQED(N, Q2I, Q2F, solns, solsg);

      cout << "Mellin variable N = " << N << endl;
      cout << "Q2I = " << Q2I << " GeV^2" << endl;
      cout << "Q2F = " << Q2F << " GeV^2" << endl;

      cout << "\nSinglet solution:" << endl;
      cout.precision(15);

      for (int i = 0; i < 3; i++) {
        cout << fixed << "[\t";
        for (int j = 0; j < 3; j++)
          cout << solsg[i][j].real() << "\t";
        cout << "]" << endl;
       }

      cout << endl;

      cout << "Non singlet solution:" << endl;
      for (int n = 0; n < 10; n++)
        cout << solns[n].real() << endl;

      cout << endl;
    }

  // Benchmarking computation time
  clock_gettime(CLOCK_MONOTONIC, &finish);
  elapsed = (finish.tv_sec - start.tv_sec);
  elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
  printf("Time elapsed: %f seconds.\n", elapsed);

  for (int i = 0; i < 3; i++)
    if (solsg[i]) delete[] solsg[i];
  delete[] solsg;

  delete[] solns;
  delete evol;

  return 0;
}
