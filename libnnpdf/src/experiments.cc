// $Id: experiments.cc 3222 2015-09-09 16:52:06Z stefano.carrazza@mi.infn.it $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <cmath>

#include "NNPDF/experiments.h"
#include "NNPDF/pdfset.h"
#include "NNPDF/dataset.h"
#include "NNPDF/thpredictions.h"
#include "NNPDF/utils.h"
#include "NNPDF/randomgenerator.h"
#include "NNPDF/exceptions.h"

using namespace std;
using namespace NNPDF;

/**
  * Constructor
  */
Experiment::Experiment(std::vector<DataSet> const& sets, std::string const& expname):
fExpName(expname),
fNData(0),
fNSys(0),
fData(NULL),
fT0Pred(NULL),
fCovMat(NULL),
fInvCovMat(NULL),
fStat(NULL),
fSys(NULL),
fSetSysMap(NULL),
fIsArtificial(false),
fIsClosure(false),
fIsT0(false)
{
  // Check
  if (sets.size() == 0)
    throw LengthError("Experiment::Experiment","No DataSets in constructor!");

  fIsT0 = sets[0].IsT0();

  // Create datasets
  for (size_t i=0; i< sets.size(); i++)
  {
    fSets.push_back(sets[i]);
    if (sets[i].IsT0() != fIsT0)
      throw UserError("Experiment::Experiment", "Supplied datasets must all be either T0 or EXP");
  }
  
  // Pull data from datasets
  PullData();
  GenCovMat();
}

/**
  * Copy-Constructor
  */
Experiment::Experiment(Experiment const& exp):
fExpName(exp.fExpName),
fNData(exp.fNData),
fNSys(exp.fNSys),
fData(NULL),
fT0Pred(NULL),
fCovMat(NULL),
fInvCovMat(NULL),
fStat(NULL),
fSys(NULL),
fSetSysMap(NULL),
fIsArtificial(exp.fIsArtificial),
fIsClosure(exp.fIsClosure),
fIsT0(exp.fIsT0)
{
  // Copy datasets
  for (size_t i=0; i< exp.fSets.size(); i++)
    fSets.push_back(exp.fSets[i]);
  
  // Pull data from datasets
  PullData();
  GenCovMat();
}

/**
  * Copy-Constructor
  */
Experiment::Experiment(Experiment const& exp, std::vector<DataSet> const& sets):
fExpName(exp.fExpName),
fNData(exp.fNData),
fNSys(exp.fNSys),
fData(NULL),
fT0Pred(NULL),
fCovMat(NULL),
fInvCovMat(NULL),
fStat(NULL),
fSys(NULL),
fSetSysMap(NULL),
fIsArtificial(exp.fIsArtificial),
fIsClosure(exp.fIsClosure),
fIsT0(exp.fIsT0)
{
  // Copy datasets
  for (size_t i=0; i< sets.size(); i++)
    fSets.push_back(sets[i]);
  
  // Pull data from datasets
  PullData();
  GenCovMat();
}

/**
  * Destructor
  */
Experiment::~Experiment()
{
  ClearLocalData();  
  ClearLocalCovMat();
  fSets.clear();
}

std::string const& Experiment::GetSetName(int i) const { return fSets[i].GetSetName();};

/*
 * Clears data pulled from DataSets
 */
void Experiment::ClearLocalData()
{ 
  // Already clear
  if (!fData)
    return;

  delete[] fData;
  delete[] fStat;
  delete[] fT0Pred;
  
  for (int s = 0; s < GetNSet(); s++)
    delete[] fSetSysMap[s];
    
  delete[] fSetSysMap;
  
  for (int i = 0; i < fNData; i++)
    delete[] fSys[i];
  delete[] fSys;
  
  fNData = 0;
  fNSys = 0; 
}

/*
 * Clears covariance matrices
 * This method is intentionally separate from ClearLocalData
 */
void Experiment::ClearLocalCovMat()
{
  // Already clear
  if (!fCovMat)
    return;

  for (int i = 0; i < fNData; i++)
  {
    delete[] fCovMat[i];
    delete[] fInvCovMat[i];
  }

  delete[] fCovMat;
  delete[] fInvCovMat;
}

/**
 * @brief Experiment::MakeReplica
 * Produces the shifts by using a MC sampling
 * Modifies the experiment and correlated datasets
 */
void Experiment::MakeReplica()
{
  cout << "-- Generating replica data for " << fExpName << endl;

  if (fNData == 0)
    throw RangeError("Experiment::MakeReplica","you must run ReadData before making a replica.");

  RandomGenerator* rng = RandomGenerator::GetRNG();

  double *rand  = new double[fNSys];  
  double *xnor  = new double[fNData];
  double *artdata = new double[fNData];
  sysType rST[2] = {ADD,MULT};

  // generate procType array for ease of checking
  std::string *proctype = new std::string[fNData];
  int counter = 0;
  for (int s = 0; s < GetNSet(); s++)
    for (int i=0; i< GetSet(s).GetNData(); i++)
      proctype[counter++] = GetSet(s).GetProc(i);


  // Generate positive defined replicas
  bool isArtNegative = true;
  while (isArtNegative)
  {
    isArtNegative = false;
    for (int l = 0; l < fNSys; l++)
    {
      rand[l] = rng->GetRandomGausDev(1.0);
      if (fSys[0][l].isRAND)
        fSys[0][l].type = rST[rng->GetRandomUniform(2)];
      for (int i = 1; i < fNData; i++)
        fSys[i][l].type = fSys[0][l].type;
    } 
    
    for (int i = 0; i < fNData; i++) // should rearrange to update set-by-set -- nh
    {
      double xstat = rng->GetRandomGausDev(1.0)*fStat[i];   //Noise from statistical uncertainty
      
      double xadd = 0;
      xnor[i] = 1.0;
      
      for (int l = 0; l < fNSys; l++)
      {
        if (fSys[i][l].name.compare("UNCORR")==0)         //Noise from uncorrelated systematics
        {
          switch (fSys[i][l].type)
          {
            case ADD: xadd += rng->GetRandomGausDev(1.0)*fSys[i][l].add; break;
            case MULT: xnor[i] *= (1.0 + rng->GetRandomGausDev(1.0)*fSys[i][l].mult*1e-2); break;
          }
        }
        else                                              //Noise from correlated systematics
        {
          switch (fSys[i][l].type)
          {
            case ADD: xadd += rand[l]*fSys[i][l].add; break;
            case MULT: xnor[i] *= (1.0 + rand[l]*fSys[i][l].mult*1e-2); break;
          }
        }
      }
      artdata[i] = xnor[i] * ( fData[i] + xadd + xstat );
      
      // Only generates positive artificial data (except for closure tests and asymmetry data)
      if (artdata[i] < 0 && !fIsClosure && proctype[i].find("ASY") == std::string::npos )
      {
        isArtNegative = true;
        break;
      }
    }
  }

  // Update data in set
  int index = 0;
  for (int s = 0; s < GetNSet(); s++)
  {
    sysType *type = new sysType[fSets[s].GetNSys()];
    for (int l = 0; l < fSets[s].GetNSys(); l++)
      type[l] = fSys[0][fSetSysMap[s][l]].type; 
    fSets[s].UpdateData(artdata+index, type);
    fSets[s].SetArtificial(true);
    index+=fSets[s].GetNData();
    delete[] type;
  }
  
  // Update local data and recompute covariance matrix
  PullData();
  //GenCovMat();   //Use standard t0 covariance matrix for all replicas
  
  delete[] rand;
  delete[] xnor;
  delete[] artdata;
  delete[] proctype;
  
  // Now the fData is artificial
  fIsArtificial = true;
}


void Experiment::MakeClosure(const vector<ThPredictions>& predictions, bool const& noise)
{
  cout << "-- Generating closure data for " << fExpName << endl;

  // Set closure flag
  fIsClosure = true;

  for (size_t s = 0; s < predictions.size(); s++)
  {
    auto & theory = predictions[s];
    auto & set = fSets[s];
    auto newdata  =  vector<double>(set.GetNData());

    for (int i = 0; i < set.GetNData(); i++)
      newdata[i] = theory.GetObsCV(i);

    set.UpdateData(newdata.data()); // MakeClosure treated as shifts rather than normalisations

  }

  // Rebuild uncertainty breakdowns
  PullData();

  // Add fluctations in data according to uncertainties.
  if (noise)
    MakeReplica();

  // Rebuild covariance matrix
  GenCovMat();
}

void Experiment::MakeClosure(PDFSet* pdf, bool const& noise)
{
  vector<ThPredictions> predictions;
  for (auto& set: fSets)
  {
    predictions.emplace_back(pdf,&set);

  }
  MakeClosure(predictions, noise);
}

/**
 * Pulls data from held datasets and allocates covariance matrices
 */
void Experiment::PullData()
{
  // Clear local arrays
  ClearLocalData();
  
  fNData = 0;
  fNSys = 0;

  for (int s = 0; s < GetNSet(); s++)
  {
    fNData += fSets[s].GetNData();
    fNSys += fSets[s].GetNSys();
  }
  
  fData = new double[fNData];
  fStat = new double[fNData];
  fT0Pred = new double[fNData];
  
  int idat = 0;
  for (int s = 0; s < GetNSet(); s++)
    for (int i = 0; i < fSets[s].GetNData(); i++)
    {
      fData[idat]=fSets[s].GetData(i);
      fStat[idat]=fSets[s].GetStat(i);
      fT0Pred[idat]=fSets[s].GetT0Pred(i);
      idat++;
    }
  
  // Uncertainty breakdown
  sysError** TempSys = new sysError*[fNData];
  for (int i=0; i<fNData; i++)
    TempSys[i] = new sysError[fNSys];
    
  fSetSysMap = new int*[GetNSet()];
  for (int s = 0; s < GetNSet(); s++)
    fSetSysMap[s] = new int[fSets[s].GetNSys()];
      
  int DataIndex = 0;
  int SysIndex = 0;
  vector<sysError> sysdef;
  vector<int> sysnumber;

  for (int s = 0; s < GetNSet(); s++)
  {
    for (int l = 0; l < fSets[s].GetNSys(); l++)
      {
        int sysplace = l+SysIndex;
        sysError testsys = fSets[s].GetSys(0,l);
        if (testsys.name.compare("UNCORR")!=0 && testsys.name.compare("CORR")!=0)        // Check for specially named systematics
        {
          for (size_t j = 0; j < sysdef.size(); j++)
            if (sysdef[j].name.compare(testsys.name)==0)              // Check whether a systematic of that name has been seen yet
              {
                sysplace = sysnumber[j];
                if (testsys.type != sysdef[j].type || testsys.isRAND != sysdef[j].isRAND)
                  throw RangeError("Experiment::PullData","Systematic " + testsys.name + " definition not consistant between datasets");
              }
                
          if (sysplace == l+SysIndex)
          {
            sysdef.push_back(testsys);              // If new, add to list of systematics
            sysnumber.push_back(sysplace);
          }
          else 
          {
            SysIndex--;                             // If a repeat, need to decrement SysIndex and fNSys
            fNSys--;
          }
        }
        
        for (int i = 0; i < fSets[s].GetNData(); i++)
          TempSys[i+DataIndex][sysplace] = fSets[s].GetSys(i,l);         // Read sys
        fSetSysMap[s][l] = sysplace;                                      // Build map (needed for randomization)
        
        for (int i = 0; i < fNData; i++)             //Need to match type, random and name of zero systematics
        {
          TempSys[i][sysplace].type = testsys.type;
          TempSys[i][sysplace].isRAND = testsys.isRAND;
          TempSys[i][sysplace].name = testsys.name;
        }
      }     
    DataIndex+=fSets[s].GetNData();
    SysIndex+=fSets[s].GetNSys();
  }
  
  fSys = new sysError*[fNData];           // Build correct size systemtics array
  for (int i=0; i<fNData; i++)
  {
    fSys[i] = new sysError[fNSys];  
    for (int l = 0; l < fNSys; l++)   
      fSys[i][l] = TempSys[i][l];
    delete[] TempSys[i];
  }
  
  delete[] TempSys;
  
  return;
  
}

/**
 * Generate covariance matrix and inverse
 */
void Experiment::GenCovMat()
{
  ClearLocalCovMat();
  
  // Allocate arrays 
  fCovMat = new double*[fNData];
  fInvCovMat = new double*[fNData];
  
  for (int i=0; i<fNData; i++)
  {
    fCovMat[i] = new double[fNData];
    fInvCovMat[i] = new double[fNData];
    
    for (int j=0; j<fNData; j++)
    {
      fCovMat[i][j] = 0;
      fInvCovMat[i][j] = 0;
    }
  }


  for (int i = 0; i < fNData; i++)
    for (int j = 0; j < fNData; j++)
    {
      double sig    = 0.0;
      double signor = 0.0;
      
      if (i == j)
        sig += fStat[i]*fStat[i]; // stat error

      for (int l = 0; l < fNSys; l++)
        if (i == j || fSys[i][l].name.compare("UNCORR")!=0)
          switch (fSys[i][l].type)
          {
            case ADD: sig += fSys[i][l].add*fSys[j][l].add; break; // additive systematics
            case MULT: signor += fSys[i][l].mult*fSys[j][l].mult; break; // multiplicative systematics
          }
          
      fCovMat[i][j] = sig + signor*fT0Pred[i]*fT0Pred[j]*1e-4;
    }
    
  InvertLU(fNData, fCovMat, fInvCovMat);
}

void Experiment::ExportCovMat(string filename)
{
  ofstream outCovMat(filename.c_str());
  outCovMat << setprecision(5);
  outCovMat << scientific;
  for (int i=0; i<fNData; i++)
  {
    for (int j=0; j<fNData; j++)
      outCovMat << fCovMat[i][j] << "\t";
    outCovMat <<endl;
  }
  outCovMat.close();
  return;
}

void Experiment::ExportInvCovMat(string filename)
{
  ofstream outCovMat(filename.c_str());
  outCovMat << setprecision(5);
  outCovMat << scientific;
  for (int i=0; i<fNData; i++)
  {
    for (int j=0; j<fNData; j++)
      outCovMat << fInvCovMat[i][j] << "\t";
    outCovMat <<endl;
  }
  outCovMat.close();
  return;
}
