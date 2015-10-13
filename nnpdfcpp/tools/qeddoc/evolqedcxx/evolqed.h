//
// The NNPDF++ Collaboration - 2012
// Author: Stefano Carrazza, stefano.carrazza@mi.infn.it
//

#ifndef EVOLQED_H
#define EVOLQED_H

/**
 *  \class EvolQED
 *  \brief Class for QED evolution.
 */

#include <complex>
using namespace std;

typedef complex<double> dcomplex;

class EvolQED {

private:
  int    fNF;             //!< The number of quark flavours, typically 3, 4, 5 or 6
  int    fAlgorithm;      //!< The algorithm (0 == numerical diag., 1 == path ordering)
  int    fNINT;           //!< Number of interactions for path ordering
  int    fNTRUNC;         //!< Path ordering exponential truncation
  bool   fAlphaRunning;   //!< Activate/Deactivate the qed alpha running
  double fQ20Ref;         //!< The reference Q20 (tau mass by default)
  double fAlphaQ20Ref;    //!< The reference Alpha_QED(Q20) by default Alpha_QED(tau mass)

  //! Compute the alpha running
  double AlphaRunning(double Q2);

  //! Compute the Beta function
  double Beta0();

  //! Compute the Anomalous dimension
  dcomplex getADqq(dcomplex N);
  dcomplex getADfq(dcomplex N);
  dcomplex getADff();
  dcomplex getADqf(dcomplex N);

  //! Digamma complex function for getADqq()
  dcomplex Digamma(dcomplex x);

  //! Analytical Diagonalization algorithm
  void ComputeEigenValuesVectors(int n, dcomplex **A,
                                 dcomplex *D, dcomplex **V, dcomplex **VInv);

public:
  EvolQED(int nf, bool running, int algorithm = 0); //!< Possibles nf = [3,4,5,6]
  ~EvolQED();                                       //!< The destructor

  //! Compute the QED solution factors for a fixed N
  void EvolFactNQED(dcomplex N, double Q2I, double Q2F,
                    dcomplex* EFNNS, dcomplex** EFNSG);

  // Set Methods
  void SetNF(int nf) { fNF = nf; } //!< Sets the number of flavours

  //! Activate and deactivate the running coupling in QED
  void ActivateRunning(bool run) { fAlphaRunning = run; }

  //! Set the reference scale and coupling constant
  void SetRefCoupling(double Q20, double AlphaQ20) { fQ20Ref = Q20; fAlphaQ20Ref = AlphaQ20; }

  // Get Methods
  int    GetNF()       const { return fNF;          } //!< Returns the current number of flavours
  double GetQ20Ref()   const { return fQ20Ref;      } //!< Returns the current reference scale
  double GetAlphaQ20() const { return fAlphaQ20Ref; } //!< Returns the current coupling constant  
};

/* C wrapper interfaces to C++ routines */
extern "C" {
  EvolQED *EvolQED__new(int nf, bool run);
  void EvolQED__delete(EvolQED *This);

  void EvolQED__EvolFactNQED(EvolQED *This, dcomplex N, double Q2I, double Q2F,
                              dcomplex* EFNNS, dcomplex** EFNSG);

  void EvolQED__SetNF(EvolQED *This, int nf);
  void EvolQED__ActivateRunning(EvolQED *This, bool run);
  void EvolQED__SetRefCoupling(EvolQED *This, double Q20, double AlphaQ20);

  int    EvolQED__GetNF(EvolQED *This);
  double EvolQED__GetQ20Ref(EvolQED *This);
  double EvolQED__GetAlphaQ20(EvolQED *This);
}


#endif // EVOLQED_H
