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
#include <numeric>

#include "NNPDF/experiments.h"
#include "NNPDF/chisquared.h"
#include "NNPDF/pdfset.h"
#include "NNPDF/dataset.h"
#include "NNPDF/thpredictions.h"
#include "NNPDF/utils.h"
#include "NNPDF/randomgenerator.h"
#include "NNPDF/exceptions.h"

using namespace std;
namespace NNPDF{

matrix<double> read_total_covmat(const std::string filename) {
  // Reads in covariance matrix for an experiment from pandas dataframe
  ifstream f1;
  f1.open(filename.c_str(), ios::in);

  if(!f1.good()) {
    throw FileError("experiments", "Cannot read covariance matrix file from " + filename);
  }

  // Skip first four lines (headers)
  for (int i = 0; i < 4; i++) {
    read.ignore(fNData, '\n')
  }

  double mat[fNData][fNData];
  string line, exp_name, dataset_name;
  double row_number;
  for (int i = 0; i < fNData; i++) {
    getline(f1, line);
    istringstream lstream(line);
    // Skip over first three elements of line
    lstream >> exp_name >> dataset_name >> row_number;
    for (int j = 0; j < fNData; j++) {
      lstream >> mat[i][j];
    }
  }

  f1.close();

  return mat;
}

/**
  * Constructor
  */
Experiment::Experiment(std::vector<DataSet> const& sets, std::string const& expname):
fExpName(expname),
fNData(0),
fNSys(0),
fData(NULL),
fT0Pred(NULL),
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
fCovMat(exp.fCovMat),
fSqrtCov(exp.fSqrtCov),
fRepCovMat(exp.fRepCovMat),
fSqrtRepCov(exp.fSqrtRepCov),
fFitCovMat(exp.fFitCovMat),
fSqrtFitCov(exp.fSqrtFitCov),
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
}

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

  fNSys = 0;
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
        if (fSys[i][l].name.compare("THEORYCORR")==0) continue;   // Skip theoretical uncertainties
        if (fSys[i][l].name.compare("THEORYUNCORR")==0) continue; // Skip theoretical uncertainties
        if (fSys[i][l].name.compare("SKIP")==0) continue;         // Skip uncertainties
        if (fSys[i][l].name.compare("UNCORR")==0)                 // Noise from uncorrelated systematics
        {
          switch (fSys[i][l].type)
          {
            case ADD: xadd += rng->GetRandomGausDev(1.0)*fSys[i][l].add; break;
            case MULT: xnor[i] *= (1.0 + rng->GetRandomGausDev(1.0)*fSys[i][l].mult*1e-2); break;
            case UNSET: throw RuntimeException("Experiment::MakeReplica", "UNSET systype encountered");
          }
        }
        else                                                      // Noise from correlated systematics
        {
          switch (fSys[i][l].type)
          {
            case ADD: xadd += rand[l]*fSys[i][l].add; break;
            case MULT: xnor[i] *= (1.0 + rand[l]*fSys[i][l].mult*1e-2); break;
            case UNSET: throw RuntimeException("Experiment::MakeReplica", "UNSET systype encountered");
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

void Experiment::SetT0(const PDFSet& pdf){
  for (auto & set : fSets){
    set.SetT0(pdf);
  }
  fIsT0 = true;
  PullData();
  GenCovMat();
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
  fNData = 0; fNSys = 0;

  // List of types that do not need to be correlated across the experiment
  const std::array<std::string,5> special_types = {
      "UNCORR",
      "CORR",
      "THEORYUNCORR",
      "THEORYCORR",
      "SKIP"
  };

  // We need to correlate all subset systematics that share a name.
  // Therefore we now loop over the subsets, and form a total systematics vector,
  // taking into account these correlations. Futhermore we build a map (fSetSysMap)
  // which maps each dataset systematic to a correponding systematic in the experiment
  // (i.e an index of total_systematics)
  vector<sysError> total_systematics;
  fSetSysMap = new int*[GetNSet()];
  for ( int i = 0; i< GetNSet(); i++ )
  {
    DataSet& subset = fSets[i];
    fNData += subset.GetNData();               // Count number of datapoints
    fSetSysMap[i] = new int[subset.GetNSys()]; // Initialise map
    for (int l = 0; l < subset.GetNSys(); l++)
    {
        const sysError& testsys = subset.GetSys(0,l);
        // Check for systematics that are uncorrelated across the experiment
        if ( find(special_types.begin(), special_types.end(), testsys.name) != special_types.end() )
        {
            // This systematic is an unnamed/special type, add it to the map and the total systematics vector
            fSetSysMap[i][l] = total_systematics.size();
            total_systematics.emplace_back(testsys);
        }
        else
        {
          // This is a named systematic, we need to cross-check against existing named systematics
          bool found_systematic = false;
          for ( size_t k=0; k < total_systematics.size(); k++ )
          {
            sysError& existing_systematic = total_systematics[k];
            if (testsys.name == existing_systematic.name )
              {
                // This already exists in total_systematics, it's not a new systematic
                found_systematic = true;
                fSetSysMap[i][l] = k;
                // Perform some consistency checks
                if (testsys.type != existing_systematic.type || testsys.isRAND != existing_systematic.isRAND)
                  throw RangeError("Experiment::PullData","Systematic " + testsys.name + " definition not consistant between datasets");
                break;
              }
          }
          // If the systematic doesn't already exist in the list, add it
          if ( found_systematic == false )
          {
            fSetSysMap[i][l] = total_systematics.size();
            total_systematics.emplace_back(testsys);
          }
        }
    }
  }

  // Initialise data arrays
  fData   = new double[fNData];
  fStat   = new double[fNData];
  fT0Pred = new double[fNData];
  fSqrtWeights = vector<double>();
  fSqrtWeights.reserve(fNData);


  fNSys   = total_systematics.size();
  fSys    = new sysError*[fNData];

  // Fill the experiment arrays with the values from its subsets
  int idat = 0;
  for (int s = 0; s < GetNSet(); s++)
    for (int i = 0; i < fSets[s].GetNData(); i++)
    {
      DataSet& subset = fSets[s];
      fData[idat]   = subset.GetData(i);
      fStat[idat]   = subset.GetStat(i);
      fT0Pred[idat] = subset.GetT0Pred(i);
      //The cast is so that it does not complain about an ambiguous resolution,
      //since apparently there is sqrt(float) and sqrt(long double) but not sqrt(double)
      fSqrtWeights.push_back(std::sqrt((long double)subset.GetWeight()));

      // Loop over experimental systematics, find if there is a corresponding dataset systematic
      // and set it accordingly.
      fSys[idat] = new sysError[fNSys];
      for (int l = 0; l < fNSys; l++)
      {
        fSys[idat][l].name   = total_systematics[l].name;
        fSys[idat][l].type   = total_systematics[l].type;
        fSys[idat][l].isRAND = total_systematics[l].isRAND;
        // Find the dataset systematic corresponding to experimental systematic 'l'
        // If it doesn't exist (dataset 'subset' doesn't have this systematic) this
        // will be subset.GetNSys()
        int* map = fSetSysMap[s];
        const int dataset_systematic = find(map, map+subset.GetNSys(), l) - map;
        if (dataset_systematic == subset.GetNSys())
        {
            // This subset doesn't have this systematic, set it to zero
            fSys[idat][l].mult = 0;
            fSys[idat][l].add  = 0;
        } else {
            // Copy the systematic
            fSys[idat][l].mult = subset.GetSys(i, dataset_systematic).mult;
            fSys[idat][l].add  = subset.GetSys(i, dataset_systematic).add;
        }
      }
      idat++;
    }

  return;
}

/**
 * Generate covariance matrix and inverse
 */
void Experiment::GenCovMat()
{
  fCovMat.clear();
  fCovMat.resize(fNData, fNData, 0);

  for (int i = 0; i < fNData; i++) {
    // Diagonal case

    fCovMat(i, i) = fStat[i] * fStat[i]; // stat error
  }
  for (int l = 0; l < fNSys; l++) {
    auto &compsys = fSys[0][l];
    if (compsys.name == "SKIP") {
      continue;
    }
    bool iscorrelated =
        (compsys.name != "UNCORR" && compsys.name != "THEORYUNCORR");
    for (int i = 0; i < fNData; i++) {
      double diagsig = 0.0;
      double diagsignor = 0.0;
      auto &sys = fSys[i][l];
      switch (compsys.type) {
      case UNSET: throw RuntimeException("Experiment::GenCovMat", "UNSET systype encountered");
      case ADD:
        diagsig += sys.add * sys.add;
        break; // additive systematics
      case MULT:
        diagsignor += sys.mult * sys.mult;
        break; // multiplicative systematics
      }
      fCovMat(i,i) += diagsig + diagsignor * fT0Pred[i] * fT0Pred[i] * 1e-4;

      // No need to loop over the nondiagonal parts
      if (!iscorrelated) {
        continue;
      }
      for (int j = 0; j < i; j++) {
        auto &othersys = fSys[j][l];
        // Hopefully easy enough for the compiler to fuse this up
        decltype(sys.add) res;
        switch (compsys.type) {
        case UNSET: throw RuntimeException("Experiment::GenCovMat", "UNSET systype encountered");
        case ADD:
          res = sys.add * othersys.add;
          break; // additive systematics
        case MULT:
          res = sys.mult * othersys.mult * fT0Pred[i] * fT0Pred[j] * 1e-4;
          break; // multiplicative systematics
	default:
	  throw NNPDF::RuntimeException("Experiment::GenCovMat", "sys type not recognized");
	  break;
        }
        fCovMat(i, j) += res;
        fCovMat(j, i) += res;
      }
    }
  }
  //Hopefully the compiler will merge the loops?
  for(int i=0; i<fNData; i++){
      for(int j=0; j<=i; j++){
          auto w = fSqrtWeights[i]*fSqrtWeights[j];
          fCovMat(i,j) /= w;
          if(i!=j){
              fCovMat(j,i) /= w;
          }
      }
  }

  fSqrtCov = ComputeSqrtMat(fCovMat);
}

/**
* Read in covariance matrix for replica generation from file, and generate covariance matrix and its square root
*/
void Experiment::LoadRepCovMat(string filename)
{
  fRepCovMat.clear();
  fRepCovMat.resize(fNData, fNData, 0);  

  fRepCovMat = read_total_covmat(filename);
  fSqrtRepCov = ComputeSqrtMat(fRepCovMat);
}

/**
* Read in covariance matrix to be used in fit from file, and generate covariance matrix and its square root
*/
void Experiment::LoadFitCovMat(string filename)
{
  fFitCovMat.clear();
  fFitCovMat.resize(fNData, fNData, 0);  

  fFitCovMat = read_total_covmat(filename);
  fSqrtFitCov = ComputeSqrtMat(fRepFitMat);
}

void Experiment::ExportCovMat(string filename)
{
  ofstream outCovMat(filename.c_str());
  outCovMat << setprecision(5);
  outCovMat << scientific;
  for (int i=0; i<fNData; i++)
  {
    for (int j=0; j<fNData; j++)
      outCovMat << fCovMat(i,j) << "\t";
    outCovMat <<endl;
  }
  outCovMat.close();
  return;
}

void Experiment::ExportSqrtCov(string filename)
{
  ofstream outCovMat(filename.c_str());
  outCovMat << setprecision(5);
  outCovMat << scientific;
  for (int i=0; i<fNData; i++)
  {
    for (int j=0; j<fNData; j++)
      outCovMat << fSqrtCov(i, j) << "\t";
    outCovMat <<endl;
  }
  outCovMat.close();
  return;
}

//___________________________________________________
vector<Experiment *> pseudodata(vector<Experiment *> const &exps,
                                unsigned long int dataseed, int replica)
{
  // make a copy of the experiments
  vector<Experiment *> output;
  output.reserve(exps.size());
  for (auto &e : exps) {
    output.emplace_back(new Experiment(*e));
  }

  // select the appropriate random seed, using dataseed
  // as initial condition and replica for the filtering selection
  unsigned long int seed = 0;
  RandomGenerator::GetRNG()->SetSeed(dataseed);
  for (int i = 0; i < replica; i++) {
    seed = RandomGenerator::GetRNG()->GetRandomInt();
  }
  RandomGenerator::GetRNG()->SetSeed(seed);
  for (auto e : output)
  {
    // take exps and MakeReplica
    e->MakeReplica();

    // keep rng flow in sync with nnfit by calling tr/val random seed shuffle
    for (auto &set : e->DataSets()) {
      // Creating Masks
      vector<int> mask(set.GetNData());
      std::iota(mask.begin(), mask.end(), 0);
      RandomGenerator::GetRNG()->ShuffleVector(mask);
    }
  }

  return output;
}

} // namespace: NNPDF
