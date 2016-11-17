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
#include "svn.h"

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

// Signal catcher
static void catch_int(int signal) {
  if ( state == FIT_ITER )
  {
    cout << Colour::FG_RED << endl << endl << "----------------- SIGINT - Interrupting fit ----------------- "<<endl<<endl;
    state = FIT_ABRT;
  }
  else
    exit(-1);
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
      string filename;
      int replica = 0;
      if (argc > 2)
        {
          replica = atoi(argv[1]);
          filename.assign(argv[2]);
        }
      else
        {
          cout << Colour::FG_RED << "\nusage: nnfit [replica number] [configuration filename]\n" << endl;
          exit(-1);
        }

      if (replica <= 0)
        {
          cerr << Colour::FG_RED << "\nError replica must be > 0" << endl;
          exit(-1);
        }

      // Creates the configuration class
      NNPDFSettings settings(configPath() + filename);

      // Creating output folder
      CreateResultsFolder(settings, replica);

      // Configuration file verification
      stringstream targetFile;
      targetFile << "nnfit/replica_"<<replica<<"/nnfit.yml";

      settings.PrintConfiguration(targetFile.str());
      settings.VerifyConfiguration(targetFile.str());

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

      // Setting the correct replica based seed
      unsigned long int seed = 0;
      for (int i = 0; i < replica; i++)
        seed = RandomGenerator::GetRNG()->GetRandomInt();
      RandomGenerator::GetRNG()->SetSeed(seed);

      // Creating folders
      cout << "\n- MC generation of replica " << replica << ":" << endl;

      // Read data and perform training-validation split
      cout << Colour::FG_YELLOW << "-----------------  Performing training - validation split ----------------- " << Colour::FG_DEFAULT << endl;
      vector<Experiment*> Exp;
      vector<Experiment*> training;
      vector<Experiment*> validation;

      PDFSet* T0Set = new LHAPDFSet(settings.Get("datacuts","t0pdfset").as<string>(), PDFSet::ER_MCT0);
      for (int i = 0; i < settings.GetNExp(); i++)
        {
          if (settings.GetExpName(i) == "REWEIGHT") // Don't fit RW experiment
            continue;

          vector<DataSet> datasets;
          for (int j = 0; j < settings.GetExpSets(i).size(); j++)
            {
              datasets.push_back(LoadDataSet(settings, settings.GetExpSets(i)[j], DATA_FILTERED));
              MakeT0Predictions(T0Set,datasets[j]);
            }

          Experiment *exp = new Experiment(datasets, settings.GetExpName(i));

          if (settings.Get("closuretest","printpdf4gen").as<bool>())
            Exp.push_back(new Experiment(datasets, settings.GetExpName(i)));

          // Apply MC shifts
          if (settings.Get("fitting","genrep").as<bool>())
            exp->MakeReplica();

          training.push_back(NULL);
          validation.push_back(NULL);

          TrainValidSplit(settings, exp, training.back(), validation.back());

          delete exp;
        }

      // Read Positivity Sets
      if (settings.GetNPos())
        cout << Colour::FG_YELLOW << " ----------------- Reading positivity sets ----------------- " << Colour::FG_DEFAULT << endl;

      // Positivity sets
      vector<PositivitySet> pos;
      for (int i = 0; i < settings.GetNPos(); i++)
        {
          cout << Colour::FG_BLUE << "\n- Loading: " << Colour::FG_DEFAULT << settings.GetPosName(i) << endl;
          pos.push_back(LoadPositivitySet(settings,settings.GetPosName(i),settings.GetPosInfo(settings.GetPosName(i)).tLambda));
          pos[i].SetBounds(T0Set);
        }

      cout <<endl;

      // Fit Basis
      FitBasis* fitbasis = getFitBasis(settings, NNPDFSettings::getFitBasisType(settings.Get("fitting","fitbasis").as<string>()));

      // Minimiser control
      Minimizer* minim = NULL;
      switch (NNPDFSettings::getFitMethod(settings.Get("fitting","fitmethod").as<string>()))
      {
        case MIN_GA:
          minim = new GAMinimizer(settings);
          cout  << Colour::FG_BLUE << "Minimiser: Genetic Algorithm" << Colour::FG_DEFAULT << endl;
          break;

        case MIN_NGA:
          minim = new NGAMinimizer(settings);
          cout  << Colour::FG_BLUE << "Minimiser: Genetic Algorithm w/ nodal mutations" << Colour::FG_DEFAULT << endl;
          break;

        case MIN_NGAP:
          minim = new NGAPMinimizer(settings);
          cout  << Colour::FG_BLUE << "Minimiser: Genetic Algorithm w/ nodal mutations for alpha/beta preprocessing" << Colour::FG_DEFAULT << endl;
          break;

        case MIN_NGAFT:
          minim = new NGAFTMinimizer(settings);
          cout  << Colour::FG_BLUE << "Minimiser: Genetic Algorithm w/ fixed threshold term NN(x)-NN(1)" << Colour::FG_DEFAULT << endl;
          break;

        case MIN_CMAES:
          minim = new CMAESMinimizer(settings);
          cout  << Colour::FG_BLUE << "Minimiser: CMA-ES" << Colour::FG_DEFAULT << endl;
          break;  

        default:
          cout << Colour::FG_RED << "ERROR: Invalid Minimiser" << Colour::FG_DEFAULT <<endl;
          exit(-1);
          break;
      }

      // Parametrisation control
      FitPDFSet* fitset = NULL;
      switch (NNPDFSettings::getParamType(settings.Get("fitting","paramtype").as<string>()))
      {
        case PARAM_NN:
          fitset = FitPDFSet::Generate<MultiLayerPerceptron,GAMinimizer>(settings, fitbasis); // need to rewrite generate
          cout << Colour::FG_BLUE << "Parametrisation: Neural Network" << Colour::FG_DEFAULT << endl;
          break;

        case PARAM_CHEBYSHEV:
          fitset = FitPDFSet::Generate<ChebyshevPolynomial,GAMinimizer>(settings, fitbasis); // need to rewrite generate
          cout << Colour::FG_BLUE << "Parametrisation: Chebyshev Polynomial (Order 10)" << Colour::FG_DEFAULT << endl;
          break;

        case PARAM_QUADNN:
          fitset = FitPDFSet::Generate<QuadMultiLayerPerceptron,GAMinimizer>(settings, fitbasis);
          cout << Colour::FG_BLUE << "Parametrisation: Quadratic Neural Network" << Colour::FG_DEFAULT << endl;
          break;

        case PARAM_NNP:
          fitset = FitPDFSet::Generate<MultiLayerPerceptronPreproc,GAMinimizer>(settings, fitbasis); // need to rewrite generate
          cout << Colour::FG_BLUE << "Parametrisation: Neural Network Preprocessing" << Colour::FG_DEFAULT << endl;
          break;

        default:
          cout << Colour::FG_RED << "ERROR: Invalid Parametrisation" << Colour::FG_DEFAULT << endl;
          exit(-1);
          break;
      }
      fitset->ValidateStartingPDFs();

      // Stopping criterion
      StoppingCriterion *stop = NULL;
      switch (NNPDFSettings::getStopType(settings.Get("stopping","stopmethod").as<string>()))
      {
        case STOP_NONE:
          stop = new StoppingCriterion(settings);
          cout << Colour::FG_BLUE << "Stopping Criterion: Fixed Length Fit" << Colour::FG_DEFAULT << endl;
          break;

        case STOP_GRAD:
          stop = new SimpleGradientStop(settings);
          cout << Colour::FG_BLUE << "Stopping Criterion: Gradient Stopping"<< Colour::FG_DEFAULT << endl;
          break;

        case STOP_VAR:
          stop = new SimpleVarianceStop(settings);
          cout << Colour::FG_BLUE << "Stopping Criterion: Variance Stopping" << Colour::FG_DEFAULT << endl;
          break;

        case STOP_LB:
          stop = new LookBackCV(settings);
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
      minim->Init(fitset,training, pos);

      for (int i = 0; i < settings.Get("fitting","ngen").as<int>(); i++)
      {
        // Abort signal
        if (state == FIT_ABRT)
          break;

        minim->Iterate(fitset,training, pos);

        if (stop->Stop(fitset,training,validation,pos)) break;

        if (settings.Get("debug").as<bool>())
        {
          cout <<  "Generation "<<fitset->GetNIte()<<"  ";
          time.printTime(time.stop());
          time.start();
        }

        if (i % 10 == 0)
          {
            LogChi2(settings,fitset,pos,training,validation,Exp);
            LogPDF(settings,fitset,replica);
          }
      }

      state = FIT_END;

      int doftrn = 0;
      int dofval = 0;
      real erf_val = 0;
      real erf_trn = 0;

      // Compute training error function and free training sets
      for (size_t i = 0; i < training.size(); i++)
      {
        real* theory = new real[training[i]->GetNData()];

        Convolute(fitset,training[i],theory);
        NNPDF::ComputeChi2(training[i],1,theory,&erf_trn);

        doftrn += training[i]->GetNData();

        delete[] theory;
        delete training[i];
      }

      // Compute validation error function and free validation sets
      for (size_t i = 0; i < validation.size(); i++)
        if (validation[i])
        {
          real* theory = new real[validation[i]->GetNData()];

          Convolute(fitset,validation[i],theory);
          NNPDF::ComputeChi2(validation[i],1,theory,&erf_val);

          dofval += validation[i]->GetNData();

          delete[] theory;
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

      // Free experiment array when PRINTPDF4GEN is activated
      for (size_t i = 0; i < Exp.size(); i++)
        if (Exp[i]) delete Exp[i];
      Exp.clear();

      // Compute Final Chi2
      cout << Colour::FG_YELLOW << "Final Chi2 Test" << Colour::FG_DEFAULT << endl;
      int dof = 0;
      real chi2 = 0;
      for (int i = 0; i < settings.GetNExp(); i++)
      {
        if (settings.GetExpName(i) == "REWEIGHT") // Don't fit RW experiment
          continue;

        vector<DataSet> datasets;
        for (int j = 0; j < settings.GetExpSets(i).size(); j++)
          {
            datasets.push_back(LoadDataSet(settings, settings.GetExpSets(i)[j], DATA_FILTERED));
            MakeT0Predictions(T0Set, datasets[j]);
          }

        // Load Experiments
        Experiment *exp = new Experiment(datasets, settings.GetExpName(i));
        real* theory = new real[exp->GetNData()];

        Convolute(fitset,exp,theory);
        NNPDF::ComputeChi2(exp,1,theory,&chi2);

        dof += exp->GetNData();

        delete[] theory;
        delete exp;
      }

      // Check Positivity Veto
      bool posVeto = false;
      // Compute Final Chi2
      cout << Colour::FG_BLUE << "\n- Final Positivity Test" << Colour::FG_DEFAULT << endl;
      for (int i = 0; i < pos.size(); i++)
      {
        // Load Experiments
        int res;
        pos[i].ComputeNUnacceptable(fitset,&res);
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

      // Export LHgrid part
      fitset->ExportPDF(replica, erf_val/dofval, erf_trn/doftrn, chi2/dof, posVeto);

      // Export Logs
      delete minim;
      LogManager::ExportLogs();

      // Delete t0set
      delete T0Set;

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


// Chi2 per experiment logger
void LogChi2(NNPDFSettings const& settings,
             const FitPDFSet* pdf,
             vector<PositivitySet> const& pos,
             vector<Experiment*> const& train,
             vector<Experiment*> const& valid,
             vector<Experiment*> const& exp)
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

  // Now fill the log with chi2 to real data
  if (settings.Get("closuretest","printpdf4gen").as<bool>())
    {
      stringstream logStringrd;
      logStringrd << "Generation "<<pdf->GetNIte()<<"  NExp: "<<Nexp<<endl;

      // Total chi^2 values for real data
      real ExpChi2Tot = 0;
      int nDataExp = 0;

      const int Nexprd = exp.size();

      for (int i = 0; i < Nexprd; i++)
        {
          string ExpName;
          real Exp = 0;

          // compute the chi^2 to real data
          if (exp[i] != NULL)
            {
              ExpName = exp[i]->GetExpName();
              FastAddChi2(pdf,exp[i],&Exp);
              ExpChi2Tot += Exp;
              Exp /= exp[i]->GetNData();
              nDataExp += exp[i]->GetNData();
            }

          // Write to log string
          logStringrd << "\t" << ExpName << "  " << Exp << endl;
        }

      logStringrd << "Total Exp: " << ExpChi2Tot/nDataExp << endl;

      LogManager::AddLogEntry("Chi2RealData",logStringrd.str());
    }
  
  return;
}

void LogPDF(NNPDFSettings const& settings,
            FitPDFSet* pdf,
            int replica)
{
  if (settings.Get("closuretest","printpdf4gen").as<bool>())
    LogManager::AddLogEntry("PDFgenerations",pdf->ExportPDF(replica));
}

void TrainValidSplit(NNPDFSettings const& settings,
                     Experiment* const& exp, Experiment *&tr, Experiment *&val)
{
  vector<DataSet> trainingSets;
  vector<DataSet> validationSets;

  int expValSize = 0; // size of validation experiment
  RandomGenerator *rng = RandomGenerator::GetRNG(); // Random number generator

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
