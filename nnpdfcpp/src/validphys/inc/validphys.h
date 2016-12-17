// $Id: validphys.h 1981 2014-07-30 12:49:26Z stefano.carrazza@mi.infn.it $
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#pragma once

#include "nnpdfsettings.h"
#include "datautils.h"
#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include <NNPDF/positivity.h>

using std::setw;

enum {CUR,REF,CTEQ,MSTW};

/**
  * Check for 4-sigma deviations from the mean
  */
void CheckForBadReplicas(PDFSet *pdf, vector<Experiment*> exps, vector<ExperimentResult*> res)
{
  // check for bad replicas
  Chi2Results global;
  global.fDOF = 0;
  global.fChi2Cent = 0;
  global.fChi2Avg = 0;

  global.fMembers = pdf->GetMembers();
  global.fChi2Mem = new real[global.fMembers];
  for (int n=0; n < global.fMembers; n++)
    global.fChi2Mem[n] = 0.0;

  for (int i=0; i < (int) res.size(); i++)
  {
    if (!exps[i]->GetNSet()) continue;
    global.fChi2Avg += res[i]->GetChi2Results().fChi2Avg;
    global.fChi2Cent+= res[i]->GetChi2Results().fChi2Cent;

    global.fDOF+= res[i]->GetChi2Results().fDOF;

    for (int n=0; n < global.fMembers; n++)
      global.fChi2Mem[n]+=res[i]->GetChi2Results().fChi2Mem[n];
  }

  real globalAVG = ComputeAVG(global.fMembers, global.fChi2Mem);
  real globalSTD = ComputeStdDev(global.fMembers, global.fChi2Mem);

  cout <<endl<< "Checking for 4-Sigma deviations from mean"<<endl;
  for (int i=0; i<global.fMembers; i++)
    if (global.fChi2Mem[i] > globalAVG + 4*globalSTD )
      cerr << Colour::FG_RED << "Replica " << i <<" chi2 is too large: "<<global.fChi2Mem[i]/(real)global.fDOF<< Colour::FG_DEFAULT << endl;

  cout << "All replicas tested and verified"<<endl;
  cout << "Global average: "<< globalAVG/(real)global.fDOF<<" STD: "<<globalSTD/(real)global.fDOF<<endl;

}

/**
  * Print the chi2 for experiments and datasets
  */
void printchi2(ExperimentResult *res, Experiment *exps)
{
  const float eDOF = exps->GetNData();
  cout << endl << Colour::FG_RED
       << "Experiment: " << Colour::FG_DEFAULT
       << setw(16) << exps->GetExpName()
       << "\t"
       << "Npts:    " << res->GetExperiment()->GetNData()
       << "\t"
       << "central chi2:    " << setw(8) << res->GetChi2Cent()/eDOF
       << "\t"
       << "average chi2:    " << setw(8) << res->GetChi2Avg()/eDOF
      << endl;

  for (int j = 0; j < exps->GetNSet(); j++)
  {
    DataSetResult *dr = res->GetSetResult(j);
    ThPredictions *tt = dr->GetTheory();
    const float dDOF = dr->GetChi2Results().fDOF;
    cout << Colour::FG_BLUE
         << "Dataset: " << Colour::FG_DEFAULT
         <<  setw(16) << tt->GetSetName()
         << "\t"
         << "Npts:    " << tt->GetNData()
         << "\t"
         << "central chi2:    " << setw(8) << dr->GetChi2Results().fChi2Cent/dDOF
         << "\t"
         << "average chi2:    " << setw(8) << dr->GetChi2Results().fChi2Avg/dDOF
         << "\t"
   << endl;
  }

}

/**
  * Print the total chi2 for experiments and datasets
  */
void printchi2tot(vector<ExperimentResult*> res)
{
  // compute the chi2 for exp and datasets
  real chi2exp = 0, chi2expavg = 0, chi2sets = 0, chi2setsavg = 0;
  int dofexp = 0, dofsets = 0;

  for (int i = 0; i < (int) res.size(); i++)
    {
      chi2exp += res[i]->GetChi2Cent();
      chi2expavg += res[i]->GetChi2Avg();
      dofexp += res[i]->GetDOF();

      for (int j = 0; j < res[i]->GetExperiment()->GetNSet(); j++)
        {
          chi2sets += res[i]->GetSetResult(j)->GetChi2Cent();
          chi2setsavg += res[i]->GetSetResult(j)->GetChi2Avg();
          dofsets += res[i]->GetSetResult(j)->GetDOF();
        }
    }

  cout << "=> TOTAL Chi2 for experiments:\tNpts:    "
       << dofexp
       << "\t"
       << "central chi2:    " << setw(8) << chi2exp / dofexp
       << "\t"
       << "average chi2:    " << setw(8) << chi2expavg / dofexp << endl;

  cout << "=> TOTAL Chi2 for datasets:\tNpts:    "
       << dofsets
       << "\t"
       << "central chi2:    " << setw(8) << chi2sets / dofsets
       << "\t"
       << "average chi2:    " << setw(8) << chi2setsavg / dofsets << endl;

}

/**
 * @brief CheckPositivityPoints
 * @param pos
 * @param pdf
 */
void CheckPositivityPoints(PositivitySet const& pos, const PDFSet* pdf)
{
  // loop over all positivity observables
  int *res = new int[pdf->GetMembers()];
  pos.ComputeNUnacceptable(pdf,res);
  for (int i = 0; i < pdf->GetMembers(); i++)
    {
      int violated = res[i];
      if (violated != 0)
        cerr << Colour::FG_RED << "Replica " << i+1 << " violates " << violated << " point(s)" << Colour::FG_DEFAULT << endl;
    }
  delete[] res;

  cout << " ----------------- Positivity Test End ----------------- " << endl;
  cout << endl;
}
