// $Id: nnfit.cc 1760 2014-05-06 14:56:31Z s0673800 $
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

/**
 * nnfit - creates replica and apply fit
 */

#include "nnfit.h"

#include <signal.h>
#include <NNPDF/thpredictions.h>
#include <NNPDF/logger.h>
#include <NNPDF/lhapdfset.h>
#include <NNPDF/parametrisation.h>
#include <NNPDF/timer.h>
#include <NNPDF/chisquared.h>

#include "loadutils.h"
#include "datautils.h"
#include "fitbases.h"
#include "fitpdfset.h"
#include "minimizer.h"
#include "stopping.h"
#include "fastaddchi2.h"
#include "exportgrid.h"
#include "evolgrid.h"

// Signal catcher
static void catch_int(int) {
  if ( state == FIT_ITER )
  {
    cout << Colour::FG_RED << "\n\n----------------- SIGINT - Interrupting fit -----------------\n" << Colour::FG_DEFAULT << endl;
    state = FIT_ABRT;
  }
  else
    exit(-1);
}

// Set the RNG seed from replica id
void SetSeed(int const& replica)
{
  unsigned long int seed = 0;
  for (int i = 0; i < replica; i++)
    seed = RandomGenerator::GetRNG()->GetRandomInt();
  RandomGenerator::GetRNG()->SetSeed(seed);
}

/**
 * \param argv the filename containing the configuration
 */
int main(int argc, char **argv)
{
  MPI::Init();

  if (MPI::TaskID() == 0) // master slave
    {
      // Fitting status and signal handler
      state = FIT_INIT;
      signal(SIGINT, catch_int);

      // Read configuration filename from arguments
      string folder;
      int replica = 0;
      if (argc > 2)
        {
          replica = atoi(argv[1]);
          folder.assign(argv[2]);
        }
      else
        {
          cout << Colour::FG_RED << "\nusage: nnfit [replica number] [configuration folder]\n" << Colour::FG_DEFAULT << endl;
          exit(-1);
        }

      if (replica <= 0)
        {
          cerr << Colour::FG_RED << "\nError replica must be > 0" << Colour::FG_DEFAULT << endl;
          exit(-1);
        }

      // Creates the configuration class
      NNPDFSettings settings(folder);
      settings.VerifyConfiguration();

      // Creating output folder
      CreateResultsFolder(settings, replica);

      // Initialise log manager
      stringstream logPath;
      logPath << settings.GetResultsDirectory()<<"/nnfit/replica_"<<replica<<"";
      LogManager::InitPath(logPath.str());

      // Initialise log chi^2 breakdown
      LogManager::AddLogger("Chi2Breakdown","chi2exps.log");
      if (settings.Get("closuretest","printpdf4gen").as<bool>())
        {
          LogManager::AddLogger("PDFgenerations","PDFgenerations.grid");
          LogManager::AddLogger("Chi2RealData","chi2expsrealdata.grid");
        }

      // Use 'dataseed' if it exists in the yml file otherwise
      // use as the default 'seed' (set by the NNPDFSettings constructor)
      if (settings.Exists("fitting","dataseed"))
        RandomGenerator::GetRNG()->SetSeed(settings.Get("fitting", "dataseed").as<unsigned long int>());
      SetSeed(replica);

      // Creating folders
      cout << "\n- MC generation of replica " << replica << ":" << endl;

      // Read data and perform training-validation split
      cout << Colour::FG_YELLOW << "-----------------  Performing training - validation split ----------------- " << Colour::FG_DEFAULT << endl;

      vector<Experiment*> training;
      vector<Experiment*> validation;
      vector<PositivitySet> pos;
      LoadAllDataAndSplit(settings, training, validation, pos);

      // Fit Basis
      std::unique_ptr<FitBasis> fitbasis(getFitBasis(settings, NNPDFSettings::getFitBasisType(settings.Get("fitting","fitbasis").as<string>()), replica));

      // If 'dataseed' exists then reset the RNG to 'seed' for the GA
      if (settings.Exists("fitting","dataseed"))
        {
          RandomGenerator::GetRNG()->SetSeed(settings.Get("fitting", "seed").as<unsigned long int>());
          SetSeed(replica);
        }

      // Minimiser control
      std::unique_ptr<Minimizer> minim;
      switch (NNPDFSettings::getFitMethod(settings.Get("fitting","fitmethod").as<string>()))
      {
        case MIN_GA:
          minim = std::make_unique<GAMinimizer>(settings);
          cout  << Colour::FG_BLUE << "Minimiser: Genetic Algorithm" << Colour::FG_DEFAULT << endl;
          break;

        case MIN_NGA:
          minim = std::make_unique<NGAMinimizer>(settings);
          cout  << Colour::FG_BLUE << "Minimiser: Genetic Algorithm w/ nodal mutations" << Colour::FG_DEFAULT << endl;
          break;

        case MIN_NGAFT:
          minim = std::make_unique<NGAFTMinimizer>(settings);
          cout  << Colour::FG_BLUE << "Minimiser: Genetic Algorithm w/ fixed threshold term NN(x)-NN(1)" << Colour::FG_DEFAULT << endl;
          break;

        case MIN_CMAES:
          minim = std::make_unique<CMAESMinimizer>(settings);
          cout  << Colour::FG_BLUE << "Minimiser: CMA-ES" << Colour::FG_DEFAULT << endl;
          break;

        default:
          cout << Colour::FG_RED << "ERROR: Invalid Minimiser" << Colour::FG_DEFAULT <<endl;
          exit(-1);
          break;
      }

      // Parametrisation control
      std::unique_ptr<FitPDFSet> fitset;
      switch (NNPDFSettings::getParamType(settings.Get("fitting","paramtype").as<string>()))
      {
        case PARAM_NN:
          fitset = std::unique_ptr<FitPDFSet>(FitPDFSet::Generate<MultiLayerPerceptron,GAMinimizer>(settings, fitbasis.get())); // need to rewrite generate
          cout << Colour::FG_BLUE << "Parametrisation: Neural Network" << Colour::FG_DEFAULT << endl;
          break;

        case PARAM_SLNPP:
          fitset = std::unique_ptr<FitPDFSet>(FitPDFSet::Generate<SingleLayerPerceptronPreproc,GAMinimizer>(settings, fitbasis.get())); // need to rewrite generate
          cout << Colour::FG_BLUE << "Parametrisation: Single layer network (preprocessed)" << Colour::FG_DEFAULT << endl;
          break;

        case PARAM_SLN:
          fitset = std::unique_ptr<FitPDFSet>(FitPDFSet::Generate<SingleLayerPerceptron,GAMinimizer>(settings, fitbasis.get())); // need to rewrite generate
          cout << Colour::FG_BLUE << "Parametrisation: Single layer network" << Colour::FG_DEFAULT << endl;
          break;

        default:
          cout << Colour::FG_RED << "ERROR: Invalid Parametrisation" << Colour::FG_DEFAULT << endl;
          exit(-1);
          break;
      }
      fitset->ValidateStartingPDFs();

      // Stopping criterion
      std::unique_ptr<StoppingCriterion> stop;
      switch (NNPDFSettings::getStopType(settings.Get("stopping","stopmethod").as<string>()))
      {
        case STOP_NONE:
          stop = std::make_unique<StoppingCriterion>(settings);
          cout << Colour::FG_BLUE << "Stopping Criterion: Fixed Length Fit" << Colour::FG_DEFAULT << endl;
          break;

        case STOP_LB:
          stop = std::make_unique<LookBackCV>(settings);
          cout << Colour::FG_BLUE << "Stopping Criterion: Look-Back Cross-Validation" << Colour::FG_DEFAULT << endl;
          break;

        default:
          cout << Colour::FG_RED << "ERROR: Invalid Stopping Criterion" << Colour::FG_DEFAULT << endl;
          exit(-1);
          break;
      }

      cout << endl;

      int nData = 0;
      for (size_t i=0; i<training.size(); i++)
        nData+=training[i]->GetNData();

      cout << "Training upon "<<nData<<" datapoints"<<endl;

      Timer time;
      if (settings.Get("debug").as<bool>())
      {
        time.start();
        state = FIT_ITER;
      }

      // Initialise minimiser
      minim->Init(fitset.get(),training, pos);

      for (int i = 0; i < settings.Get("fitting","ngen").as<int>(); i++)
      {
        // Abort signal
        if (state == FIT_ABRT)
          break;

        minim->Iterate(fitset.get(), training, pos);

        if (stop->Stop(fitset.get(), training, validation, pos)) break;

        if (settings.Get("debug").as<bool>())
        {
          cout <<  "Generation "<<fitset->GetNIte()<<"  ";
          time.printTime(time.stop());
          time.start();
        }

        if (i % 100 == 0)
            LogChi2(fitset.get(), pos, training, validation);
      }

      state = FIT_END;

      int doftrn = 0;
      int dofval = 0;
      real erf_val = 0;
      real erf_trn = 0;

      // Compute training error function and free training sets
      for (size_t i = 0; i < training.size(); i++)
      {
        vector<real> theory(training[i]->GetNData());

        Convolute(fitset.get(),training[i],theory.data());
        NNPDF::ComputeChi2(training[i],1,theory.data(),&erf_trn);

        doftrn += training[i]->GetNData();

        delete training[i];
      }

      // Compute validation error function and free validation sets
      for (size_t i = 0; i < validation.size(); i++)
        if (validation[i])
        {
          vector<real> theory(validation[i]->GetNData());

          Convolute(fitset.get(),validation[i],theory.data());
          NNPDF::ComputeChi2(validation[i],1,theory.data(),&erf_val);

          dofval += validation[i]->GetNData();

          delete validation[i];
        }

      // Check for empty validation set
      if (dofval == 0)
      {
        erf_val = 0;
        dofval = 1;
      }

      training.clear();
      validation.clear();

      // Compute Final Chi2
      cout << Colour::FG_YELLOW << "Final Chi2 Test" << Colour::FG_DEFAULT << endl;
      int dof = 0;
      real chi2 = 0;
      auto T0Set = std::make_unique<LHAPDFSet>(settings.Get("datacuts","t0pdfset").as<string>(), PDFSet::erType::ER_MCT0);
      for (int i = 0; i < settings.GetNExp(); i++)
      {
        if (settings.GetExpName(i) == "REWEIGHT") // Don't fit RW experiment
          continue;

        vector<DataSet> datasets;
        for (size_t j = 0; j < settings.GetExpSets(i).size(); j++)
          {
            datasets.push_back(LoadDataSet(settings, settings.GetExpSets(i)[j], DATA_FILTERED));
            MakeT0Predictions(T0Set.get(), datasets[j]);
          }

        // Load Experiments
        auto exp = std::make_unique<Experiment>(datasets, settings.GetExpName(i));
        vector<real> theory(exp->GetNData());

        Convolute(fitset.get(),exp.get(),theory.data());
        NNPDF::ComputeChi2(exp.get(),1,theory.data(),&chi2);

        dof += exp->GetNData();
      }

      // Check Positivity Veto
      bool posVeto = false;
      // Compute Final Chi2
      cout << Colour::FG_BLUE << "\n- Final Positivity Test" << Colour::FG_DEFAULT << endl;
      for (size_t i = 0; i < pos.size(); i++)
      {
        // Load Experiments
        int res;
        pos[i].ComputeNUnacceptable(fitset.get(),&res);
        if (res != 0)
          {
            cout << Colour::FG_RED << "- Positivity Vetoed\n" << Colour::FG_DEFAULT << endl;
            posVeto = true;
            break;
          }
        else
            cout << Colour::FG_GREEN << "- Passed all points for " << settings.GetPosName(i) << Colour::FG_DEFAULT << endl;
      }
      pos.clear();

      // Export meta file
      fitset->ExportMeta(replica, erf_val/dofval, erf_trn/doftrn, chi2/dof, posVeto);

      // Export fit results to an initial scale grid
      std::string basefile = settings.GetResultsDirectory() + "/nnfit/";
      std::string gridfile = basefile + "replica_" + std::to_string(replica) + "/"
                           + settings.GetPDFName() +".exportgrid";
      const auto eg = ExportGrid(*fitset, 0, replica, fitset->GetQ20());
      eg.Write(gridfile);

      // Export evolved fit
      std::string infofile = basefile + settings.GetPDFName() + ".info";
      std::string replica_file = basefile + "replica_" + std::to_string(replica) + "/"
                               + settings.GetPDFName() + ".dat";

      const vector<ExportGrid> egrid = {eg};
      auto dglapg = EvolveGrid(egrid, settings.GetTheoryMap());
      dglapg.WriteInfoFile(infofile);
      const auto outstream = dglapg.WriteLHAFile();
      write_to_file(replica_file, outstream[0].str());

      // Export Logs
      LogManager::ExportLogs();

      cout << Colour::FG_GREEN << endl;
      cout << " -------------------------------------------------\n";
      cout <<   " - nnfit completed with success" << endl;
      cout <<   " - please go "<< settings.GetResultsDirectory() << "/nnfit \n";
      cout <<   " -------------------------------------------------\n";
      cout << Colour::FG_DEFAULT << endl;

    }

  MPI::Finalize();

  return 0;
}

// Load data and perform trainng validation split
void LoadAllDataAndSplit(NNPDFSettings const& settings,
                         vector<Experiment*> & training,
                         vector<Experiment*> & validation,
                         vector<PositivitySet> & pos)
{
  auto T0Set = std::make_unique<LHAPDFSet>(settings.Get("datacuts","t0pdfset").as<string>(), PDFSet::erType::ER_MCT0);
  for (int i = 0; i < settings.GetNExp(); i++)
    {
      if (settings.GetExpName(i) == "REWEIGHT") // Don't fit RW experiment
        continue;

      vector<DataSet> datasets;
      for (int j = 0; j < (int) settings.GetExpSets(i).size(); j++)
        {
          datasets.push_back(LoadDataSet(settings, settings.GetExpSets(i)[j], DATA_FILTERED));
          MakeT0Predictions(T0Set.get(),datasets[j]);
        }

      auto exp = std::make_unique<Experiment>(datasets, settings.GetExpName(i));

      // Apply MC shifts
      if (settings.Get("fitting","genrep").as<bool>())
        exp->MakeReplica();

      training.push_back(NULL);
      validation.push_back(NULL);

      TrainValidSplit(settings, exp.get(), training.back(), validation.back());
    }

  // Read Positivity Sets
  if (settings.GetNPos())
    cout << Colour::FG_YELLOW << " ----------------- Reading positivity sets ----------------- " << Colour::FG_DEFAULT << endl;

  // Positivity sets
  for (int i = 0; i < settings.GetNPos(); i++)
    {
      cout << Colour::FG_BLUE << "\n- Loading: " << Colour::FG_DEFAULT << settings.GetPosName(i) << endl;
      pos.push_back(LoadPositivitySet(settings,settings.GetPosName(i),settings.GetPosInfo(settings.GetPosName(i)).tLambda));
      pos[i].SetBounds(T0Set.get());
    }
  cout << endl;
}

// Chi2 per experiment logger
void LogChi2(const FitPDFSet* pdf,
             vector<PositivitySet> const& pos,
             vector<Experiment*> const& train,
             vector<Experiment*> const& valid)
{

  if (pdf->GetMembers() != 1)
  {
    cerr << Colour::FG_RED << "LogChi2 Error: number of PDFs in FitPDFSet is not = 1" << Colour::FG_DEFAULT << endl;
    exit(-1);
  }

  // Total chi^2 values
  real TrnChi2Tot = 0;
  real ValChi2Tot = 0;
  real PosChi2Tot = 0;

  // Data points
  int nDataTrn = 0;
  int nDataVal = 0;

  const int Nexp = train.size();

  stringstream logString;
  logString << "Generation "<<pdf->GetNIte()<<"  NExp: "<<Nexp<<endl;

  for (int i=0; i<Nexp; i++)
  {
    string ExpName;
    real   ExpVal = 0;
    real   ExpTrn = 0;

    // Training chi^2
    if (train[i] != NULL)
    {
      ExpName = train[i]->GetExpName();
      FastAddChi2(pdf, train[i], &ExpTrn);
      TrnChi2Tot += ExpTrn;
      ExpTrn /= train[i]->GetNData();
      nDataTrn += train[i]->GetNData();
    }

    // Validation chi^2
    if (valid[i] != NULL)
    {
      ExpName = train[i]->GetExpName();
      FastAddChi2(pdf, valid[i], &ExpVal);
      ValChi2Tot += ExpVal;
      ExpVal /= valid[i]->GetNData();
      nDataVal += valid[i]->GetNData();
    }

    // Write to log string
    logString << "\t"<<ExpName<<"  "<<ExpTrn<<"  "<<ExpVal <<endl;
  }

  // Positivity info
  for (size_t i=0; i<pos.size(); i++)
    pos[i].ComputeErf(pdf,&PosChi2Tot);

  logString << "Total Trn: " << TrnChi2Tot/nDataTrn<<" Val: "<< ValChi2Tot/nDataVal<<" Positivity:"<<PosChi2Tot<<endl;

  // Push log into logger class
  LogManager::AddLogEntry("Chi2Breakdown",logString.str());

  return;
}

void TrainValidSplit(NNPDFSettings const& settings,
                     Experiment* const& exp, Experiment* &tr, Experiment* &val)
{
  vector<DataSet> trainingSets;
  vector<DataSet> validationSets;

  int expValSize = 0; // size of validation experiment

  for (int s = 0; s < exp->GetNSet(); s++)
    {
      const DataSet& set = exp->GetSet(s);

      // Fraction of data in training and validation sets
      const double trFrac = settings.GetSetInfo(set.GetSetName()).tTrainingFraction;
      const int trMax =(trFrac*set.GetNData());

      // Creating Masks
      vector<int> mask;
      for (int i = 0; i < set.GetNData(); i++) mask.push_back(i);
      RandomGenerator::GetRNG()->ShuffleVector(mask);

      const vector<int> trMaskset(mask.begin(), mask.begin() + trMax);
      const vector<int> valMaskset(mask.begin() + trMax, mask.end());

      // Initializing new datasets
      trainingSets.push_back(DataSet(exp->GetSet(s), trMaskset));
      if ((int) valMaskset.size() != 0)
      {
        validationSets.push_back(DataSet(exp->GetSet(s), valMaskset));
        expValSize += valMaskset.size();
      }
    }

  cout << Colour::FG_BLUE << "- Building Training" << Colour::FG_DEFAULT << endl;
  tr = new Experiment(*exp, trainingSets);

  cout << Colour::FG_BLUE << "- Building Validation" << Colour::FG_DEFAULT << endl;
  if (expValSize != 0)
      val = new Experiment(*exp, validationSets);
}
