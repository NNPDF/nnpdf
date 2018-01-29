// pdfs.cc
// PDF classes with RW methods

#include <sys/types.h>
#include <sys/stat.h>
#include <numeric>
#include <algorithm>
#include <cmath>

#include "pdfs.h"
#include "utils.h"

// ROOT
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TDecompBK.h"

// LHAPDF
#include "LHAPDF/LHAPDF.h"

using namespace std;

// --------------------------------------------------
//              Reweighting statics
// --------------------------------------------------

// Return unweighted replica numbers
static vector<int> unweight(
                            vector<double> const& w_old,
                            vector<int> const& repnums_old,
                            size_t const& nrep_unw
                            )
{
  double step=(1.0/(double)(nrep_unw+1.0));
  size_t nrep_old = repnums_old.size();
  
  // vector of new replica labels
  vector<int> repnums;
  
  // Vector of weight cumulants
  vector<double> wC;
  
  wC.push_back(w_old[0]/nrep_old);
  
  for (size_t i=1; i<nrep_old; i++)
    wC.push_back(wC[i-1]+w_old[i]/(double)nrep_old);
  
  if ((wC.back()-1)>1e-12)
    {  
      cerr << "Error with weight normalisation" <<endl;
      exit(1);
    }
  
  for (size_t i=1; i<nrep_unw+1; i++)
    {
      size_t j=0;
      double xpos = (double)i*step;
      
      while (xpos>wC[j])
	j++;
      
      repnums.push_back(repnums_old[j]);
    }		
  
  return repnums;
}

// Computation of PDF weights
static vector<double> computeweights(vector<double> const& chi2,vector<int> const& repnums, size_t ndata, bool norm)
{
  vector<double> logw,w;
  double exp_avg=0,wtot=0;
  size_t nrep = chi2.size();
  
  // Calculate Ln(wnn) (nn - not normalised) -> use repnums as index for unweighted PDFs
  for (size_t i=0;i<nrep;i++)
    logw.push_back( - chi2[repnums[i]]/(2.0) +( (((double) ndata)-1.0)/2.)*log(chi2[repnums[i]]));
  
  // Get maximum value of log(w)
  // Palpha doesnt cancel this term - as it doesnt normalise
  if (norm)
    exp_avg = *(max_element(logw.begin(), logw.end()));
  else
    exp_avg=0;
  
  // Calculate weights
  for (size_t i=0 ;i<nrep; i++)
    w.push_back(exp(logw[i] - exp_avg ));
  
  // Drop any weights smaller than 1e-12
  wtot=accumulate(w.begin(),w.end(),0.0); 
  
  for (size_t i=0; i<w.size(); i++)
    if ((w[i]*(nrep/wtot)) < 1e-12)
      w[i]=0;
  
  if (norm) // Normalise weights so Sum(weights)=N
    {
      wtot=accumulate(w.begin(),w.end(),0.0); 
      for (size_t i=0;i<w.size();i++)
	w[i]*=(nrep/wtot); 
    }
  
  return w;
}

// Computes total (over all data in the new set) chi2 between one observable replica and the experimental value.
// Takes an inverted CovMatrix
static double Compute_totChi2(vector<double> const& ex, vector<double> const& th, TMatrixDSym const& iCovMat)
{	
  if (ex.size() != th.size())
    {
      cerr <<"Error: (Compute_totChi2) Ex and Th vectors are not the same size"<<endl;
      exit(1);
    }
  
  double chi2tot=0;
  
  for (size_t i=0; i<ex.size(); i++)
    for (size_t j=0; j<ex.size(); j++)
      chi2tot+=(th[i]-ex[i])*(iCovMat[i][j])*(th[j]-ex[j]);
  
  // return chi2 per d.o.f    
  return chi2tot/ex.size();
}


// --------------------------------------------------
//                 PDF Constructors
// --------------------------------------------------

PDF::PDF(const PDFparams& par) 
  : nflav   (par.nflav),
    nrep    (par.nrep),
    nlha    (par.nrep),
    mult    (0),
    NAME    (par.NAME),
    xpdf    (new vector<double>[par.nrep]),
    weights (par.nrep,1),
    unweighted (false)
{
  
#ifdef VERBOSE
  cout << "PDF::PDF()" << endl;
#endif
  
  for (size_t i=0; i<nrep; i++)
    repnums.push_back(i);
  
}

PDF::PDF(const PDFparams& par, const size_t nset) 
  : nflav   (par.nflav),
    nrep    (par.nrep),
    nlha    (par.nrep),
    mult    (0),
    NAME    (par.NAME),
    xpdf    (new vector<double>[par.nrep]),
    weights (par.nrep,1),
    unweighted (false)
{
  
#ifdef VERBOSE
  cout << "PDF::PDF()" << endl;
#endif
  
  for (size_t i=0; i<nrep; i++)
    repnums.push_back(i);
  mult=nset;
  
  cout << "PDF initialized using set "
       << mult
       << " from a multiple set choice"
       << endl;
  
}
//////////////////
PDF::PDF(const PDF& pdf)
  : nflav   (pdf.nflav),
    nrep    (pdf.nrep),
    nlha    (pdf.nlha),
    mult    (pdf.mult),
    NAME    (pdf.NAME),
    weights (pdf.weights),
    repnums (pdf.repnums),
    unweighted (false)
{
  
#ifdef VERBOSE
  cout << "PDF::PDF()" << endl;
#endif
  
}

//////////////////
PDF::PDF(const PDF& pdf, const size_t Nuwrep) // Unweighting copy constructor
  : nflav   (pdf.nflav),
    nrep    (Nuwrep),
    nlha    (pdf.nlha),
    mult    (pdf.mult),
    NAME    (pdf.NAME),
    xpdf    (new vector<double>[Nuwrep]),
    weights (Nuwrep,1),
    repnums (unweight(pdf.weights,pdf.repnums,nrep)),
    unweighted (true)
{
#ifdef VERBOSE
  cout << "PDF::PDF() unweighted" << endl;
#endif
  
}


PDF::~PDF() {
  
#ifdef VERBOSE
  cout << "PDF::~PDF()" << endl;
#endif
  
  delete [] xpdf;
  
}

// --------------------------------------------------
//              PDF Public Methods
// --------------------------------------------------

void PDF::PDFGet(double &x, double &Q){
  
  switch (mult) {
  case 0:
    
    for (size_t n=0; n<nrep; n++)
      {
	LHAPDF::initPDF(repnums[n]+1);
	xpdf[n]=LHAPDF::xfx(x,Q);
      };
    
    break;
    
  default:
    
    for (size_t n=0; n<nrep; n++) 
      {
	LHAPDF::initPDFM(mult,repnums[n]+1);
	xpdf[n]=LHAPDF::xfxM(mult,x,Q);
      };  
    
    break;
  };
  
  xval=x;
  qval=Q;
  
  return;
}

void PDF::PDFInfo() 
{
  cout << "# PDF set name: " << NAME << endl;
  if (unweighted)
    cout << " ** UNWEIGHTED ** "<<endl;
  cout << "# PDF N_rep: " << nrep << endl;
  cout << "# PDF N_lha: " << nlha << endl;
  
}
// --------------------------------------------------
//              PDF Reweighting Methods
// --------------------------------------------------

void PDF::Reweight(rwparam& par)
{
  struct stat st;
  vector<double> chi2, exp, th;
  stringstream filename;
  
  switch (par.mode) {
  case 0:
    {
      cout << "RWMODE: Weight File"<<endl;
            
      if (stat(par.rwdata.c_str(),&st) == 0 )
	{
	  cout<< "Found weight data file: "<<par.rwdata.c_str()<<endl;
	  weights=ReadVector(par.rwdata.c_str(),2); 
                
	  if (weights.size()!=nlha)
	    {
	      cerr << "Error: weight file has "<<weights.size()<<" entries"<<endl;
	      cerr << "PDF consists of "<<nlha<<" replicas."<<endl;
	      exit(1);
	    }
	  
	  return;
	} else
	{
	  cerr << "Error: weight file at: "<<par.rwdata <<" cannot be found"<<endl;
	  exit(1);
	}
            
      break;
    }
    
  case 1:
    {
      cout << "RWMODE: Chi2 File"<<endl;
      
      // Check fot chi2 values file
      filename.str("");
      filename<<"../data/"<<par.rwdata<<"/"<<par.prior<<"/chi2/"<<par.rwdata<<"_"<<par.prior<<"-chi2.res";
      
      if (stat(filename.str().c_str(),&st) == 0 )
	{
	  cout<< "Found chi2 data file: "<<filename.str().c_str()<<endl;
	  cout <<"Reweighting with n_dat: "<<par.ndat<<endl;
          
	  vector<double> chi2vals=ReadVector(filename.str().c_str(),2);
          
	  if (chi2vals.size() != nlha)
	    {
	      cerr << "Error: chi2 file has "<<weights.size()<<" entries"<<endl;
	      cerr << "PDF consists of "<<nrep<<" replicas."<<endl;
	      exit(1);
	    }
	  
	  ComputeWeights(chi2vals,par.ndat);
          
	  return;
	} else
	{
	  cerr << "Error: chi2 file at: "<<filename.str() <<" cannot be found"<<endl;
	  exit(1);
	}
      
      break;
    }
    
  case 3:
    {
      cout << "RWMODE: Th Predictions"<<endl;
      
      // Check for experimental data
      filename.str("");
      filename<<"../data/"<<par.rwdata<<"/"<<par.rwdata<<"_exp.res";
      
      if (stat(filename.str().c_str(),&st) != 0 )
	{
	  cerr<< "Cannot find experimental data: "<<filename.str()<<endl;        
	  exit(0);
	}
      
      // Read experimental data  
      exp=ReadVector((filename.str()).c_str(),3);
      
      size_t edata =exp.size();
      if ( edata != par.ndat )
	cerr << "Parameter file NDATA differs from number of experimental points"<<endl;
      
      // Check for covariance matrix
      filename.str("");
      filename<<"../data/"<<par.rwdata<<"/"<<par.rwdata<<"_CovMat.res";
      
      if (stat(filename.str().c_str(),&st) != 0 )
	{
	  cerr<< "Cannot find Covariance Matrix: "<<filename.str()<<endl;        
	  exit(0);
	}
      
      //Read Covariance matrix  
      TMatrixDSym CovMat(ReadSymMatrix((filename.str()).c_str(), par.ndat));
      
      // Invert Covariance Matrix 
      TDecompBK dc(CovMat);
      TMatrixDSym iCovMat = dc.Invert();
      
      // Check for theory predictions
      filename.str("");
      filename<<"../data/"<<par.rwdata<<"/"<<par.prior<<"/obs/"<<par.rwdata<<"_"<<1001<<".res";
      
      if(stat(filename.str().c_str(),&st) != 0)
	{
	  cerr << "No suitable data found"<<endl;
	  exit(1);
	}
      
      vector<double> chi2vals;
      
      //Read theory preditions - note use nlha (total reps in prior pdf) for obs labelling.
      //Any UW related juggling is done simply in computeweights using numreps vector
      for (size_t i=1; i<(nlha+1); i++)
	{
	  th.clear();
	  filename.str("");
	  filename<<"../data/"<<par.rwdata<<"/"<<par.prior<<"/obs/"<<par.rwdata<<"_"<<(i+1000)<<".res";
	  th=ReadVector((filename.str()).c_str(),3);
          
	  if (th.size()!=par.ndat)
	    {
	      cerr << "ERR: Ndat in th: "<<filename.str()<< " does not match NDATA"<<endl;
	      exit(1);
	    }
	  
	  chi2vals.push_back(Compute_totChi2(exp,th,iCovMat));
	}
      
      // Feed chi2 vector to Computeweights
      ComputeWeights(chi2vals, par.ndat);
      
      break;
    }
    
  default:
    {
      cerr << "ERROR: Invalid RW Mode - Options are RWMODE: WGT/CHI2/DATA" <<endl;
      exit(0);
      
      break;
    }
  };
  
  return;
}

// Export LHAPDF grid
void PDF::Export(rwparam& rpar) const
{
  size_t npx=100, npq2=50;
  typedef double xdim[npx][npq2][13];
  
  // Init x and q2 grid points
  vector<double> xg, qg;
  double qmin=LHAPDF::getQ2min(0);
  double qmax=LHAPDF::getQ2max(0);
  double xmin=LHAPDF::getXmin(0);
  double xmax=LHAPDF::getXmax(0);
  
  // FORTRAN STANDARDS
  double XCH=0.1;
  
    // Set up x, q2 grids
  for (size_t i=1; i<(npx+1); i++)
    if (i<50)
      xg.push_back(xmin*pow(XCH/xmin,2.0*(((double) i ) -1)/(((double) npx ) -1)));
    else
      xg.push_back(XCH+(xmax-XCH)*(((double) i ) -51)/(((double) npx ) -51) );
  
  for (size_t i=1; i<(npq2+1); i++)
    qg.push_back(qmin*pow(qmax/qmin,(((double) i ) -1)/(((double) npq2 ) -1)));
  
  
  double oas = LHAPDF::getOrderAlphaS();
  double opdf = LHAPDF::getOrderPDF();
  
  cout <<endl<<"Writing out LHAPDF grid: "<<rpar.outfile<<endl;
  cout <<"Using LHAPDF version: "<< LHAPDF::getVersion()<<endl;


  stringstream ofilename;
  ofilename.str("");
  ofilename << rpar.outdir.c_str() << "/" << rpar.outfile.c_str();

  ofstream lhaout(ofilename.str().c_str());
  
  // Write out LHAPDF preamble
  lhaout.precision(8);
  lhaout << scientific;
  
  lhaout<<" \'Version\' \'"<<LHAPDF::getVersion()<<"\'"<<endl;
  lhaout <<" \'Description:\'"<<endl;
  
  lhaout << rpar.desc;
  
  lhaout<< " \'Alphas:\'"<<endl;
  if (oas==0.0)
    lhaout<< " \'Variable', \'lo\', \'EvolCode\'"<<endl;
  else if (oas==1.0)
    lhaout<< " \'Variable', \'nlo\', \'EvolCode\'"<<endl;
  else if (oas==2.0)
    lhaout<< " \'Variable', \'nnlo\', \'EvolCode\'"<<endl;
  else
    {
      cout <<"ERR: invalid asorder"<<endl;
      cout <<oas<<endl;
      exit(1);
    }
  
  lhaout << " 1 ,   91.2      ,   "<<LHAPDF::getQMass(4)<<"      ,   "<<LHAPDF::getQMass(5)<<"      ,   "<<LHAPDF::getQMass(6)<<endl;     
  
  lhaout<< " \'MinMax:\'"<<endl;
  lhaout <<" "<<nrep<<",  1"<<endl;
  lhaout <<" "<<LHAPDF::getXmin(0)<<" , "<<LHAPDF::getXmax(0)<<" , "<<LHAPDF::getQ2min(0)<<" , "<<LHAPDF::getQ2max(0)<<endl; 
  
  lhaout<< " \'QCDparams:\'"<<endl;
  lhaout <<" "<<nrep<<",  1"<<endl;
  lhaout << " "<<LHAPDF::getLam4(0)<<" ,  "<<LHAPDF::getLam5(0)<<endl;
  
  // alphas values for each member
  const double mz = 91.2;
  lhaout<< " \'Parameterlist:\'"<<endl;
  lhaout<< " \'list', "<<nrep<<" , "<<" 1"<<endl;
  for (size_t i=0; i<(nrep+1); i++)
    lhaout << " "<<LHAPDF::alphasPDF(mz)<<endl;
  
  // Order of evolution
  lhaout << " \'Evolution:\'"<<endl;
  
  if (opdf==0.0)
    lhaout<< " \'lo\', "<<LHAPDF::getQ2min(0)<<" , "<<1<<endl;
  else if (opdf==1.0)
    lhaout<< " \'nlo\', "<<LHAPDF::getQ2min(0)<<" , "<<1<<endl;
  else if (opdf==2.0)
    lhaout<< " \'nnlo\', "<<LHAPDF::getQ2min(0)<<" , "<<1<<endl;
  else
    {
      cerr <<"ERR: invalid pdf order"<<endl;
      cerr <<opdf<<endl;
      exit(1);
    }
  
  lhaout <<  " \'NNPDF20int\'"<<endl; 
  lhaout <<  " "<<nrep<<", "<<1<<endl;
  
  
  // Write out x, q2 grid values
  lhaout <<fixed<< " "<<npx <<endl;
  
  // x in scientific
  lhaout.precision(18);
  lhaout << scientific;
  
  for (size_t i=0; i<npx; i++)
    lhaout<<  " "<<xg[i]<<endl;
  
  lhaout <<" "<<npq2 <<endl;
  
  // Q in fixed
  lhaout.precision(18);
  lhaout << scientific;
  lhaout <<" "<<qg[0]<<endl;
  for (size_t i=0; i<npq2; i++)
    lhaout<<" "<<  qg[i]<<endl;
  
  lhaout <<" "<<nrep<<endl;
  
  // rest in fixed 
  lhaout.precision(8);
  lhaout << fixed;
  
  // Determine values of xfx, place in xfxval array, compute average over replicas for 0th member PDF
  vector<double> xfxavg;
  xdim *xfxval = new xdim[nrep];
  
  for (size_t x=0; x<npx; x++)
    {
      for (size_t q=0; q<npq2; q++)
        {
	  double xval=xg[x], qval=qg[q];
	  for (int i=-6; i<7; i++) 
            {
	      xfxavg.clear();
	      for (size_t n=0; n<nrep; n++)
                {
		  LHAPDF::initPDF(repnums[n]+1); 
		  double val=LHAPDF::xfx(xval,sqrt(qval),i);
		  if (abs(val)<1e-50)
		    val=0;
		  xfxavg.push_back(val);
		  xfxval[n][x][q][(i+6)]=val;
                }
	      
	      double avgval = favg(xfxavg);
	      if (avgval < 1e-16)
		avgval = 0;
	      
	      lhaout <<" "<<avgval<<"  ";                
            }
	  lhaout <<endl; 
        }
      cout << "Generating grid: "<<x+1<<"/"<<npx;
      cout<<"\n\033[F\033[J";
    }
  
  cout << "USING NREP: "<<nrep<<endl;
  // Write out the contents of the xfxvals array to the LHgrid
  for (size_t n=0; n<nrep; n++)
    {
      LHAPDF::initPDF(repnums[n]+1); 
		
      for (size_t x=0; x<npx; x++)
        {
	  for (size_t q=0; q<npq2; q++)
            {
	      for (int i=0; i<13; i++)
		lhaout << " "<<xfxval[n][x][q][i]<<"  ";                
	      lhaout <<endl; 
            }
        }
      cout << "Writing replica: "<<n+1<<"/"<<nrep;
      cout<<"\n\033[F\033[J";
      
    }
  
  lhaout << "'End:'"<<endl; 
  
  cout <<"LHAPDF Writeout successful" <<endl<<endl;
  delete [] xfxval;
  
  lhaout.close();
  
  return;
}

// Note vector chi2 must be chi2/d.o.f
void PDF::ComputeWeights(vector<double> _chi2, size_t _ndatrw)
{
  if (_chi2.size()!=nlha)
    {
      cerr <<"ERR: Size of chi2 vector is not the same as LHAGrid Nrep"<<endl;
      exit(1);
    }
  
  if (_ndatrw==0)
    {
      cerr <<"ERR: ndata is zero"<<endl;
      exit(1);
    }
  
  if ( repnums.size() != nlha )
    {
      cerr << "ERR: (PDF::ComputeWeights) repnums incorrectly initialised"<<endl;
      exit(1);
    }
  
  // set PDF ndat
  ndatrw=_ndatrw;
  
  // Clear prior vectors
  weights.clear();
  chi2.clear();
  
  //reintroduce d.o.f
  for (size_t i=0; i<_chi2.size(); i++)
    _chi2[i]=_chi2[i]*((double) _ndatrw);
  
  chi2=_chi2;
  
  weights=computeweights(chi2,repnums,ndatrw,true);
  return;
}

void PDF::CheckWeights() const {
  
  if (weights.size()==0)
    {
      cerr <<"ERR: weights vector is empty"<<endl;
      exit(1);
    }
  
  double wmax = *( max_element( weights.begin(), weights.end() ) ); 
  double wmin = *( min_element( weights.begin(), weights.end() ) ); 
  cout <<endl<<"******************************* "<<endl;
  cout <<"Nweights:   "<<weights.size()<<endl;
  cout <<"Sum of Weights:  "<<accumulate(weights.begin(),weights.end(),0.0)<<endl;
  cout <<"Average weight: "<<favg(weights)<<endl;;
  cout <<"Max Weight: "<<wmax<<endl;
  cout <<"Min Weight: "<<wmin<<endl;	
  cout <<"Shannon:  " <<Shannon()<<endl;
  cout <<"N_rep,eff:  " <<Neff()<<endl;
  cout <<"******************************* "<<endl<<endl;  
}

double PDF::Palpha(double alpha) const
{
  double p=0;
  
  if (chi2.size()==0)
    {
      cerr << "ERR: (PDF::palpha) chi2 values not initialised"<<endl;
      exit(1);
    }
  
  vector<double> w_alpha;
  vector<double> chi2_alpha=chi2;
  
  for (size_t i=0; i<chi2.size(); i++)
    chi2_alpha[i]=chi2_alpha[i]/(alpha*alpha);
  
  w_alpha=computeweights(chi2_alpha,repnums,ndatrw,false);
  
  for (size_t i=0; i<w_alpha.size(); i++)
    p+=w_alpha[i];
  
  return p/alpha;
}

// Computes the Shannon Entropy for a set of weights w
double PDF::Shannon() const
{    
  double N(weights.size());
  
  if (N==0)
    {
      cerr << "Error: ( PDF::Shannon() ) Weights not initialised"<<endl;
      exit(1);
    }
  
  double lNeff=0,Neff=0;
  
  for (size_t i=0; i<N; i++)
    if (weights[i]>0)
      lNeff-=weights[i]*log(weights[i]/N);
  
  Neff=exp(lNeff/N);
  
  return Neff;
}

// Computes Neff (eqn 42 RW1 paper)
double PDF::Neff() const
{
  double N(weights.size());
  
  if (N==0)
    {
      cerr << "Error: ( PDF::Shannon() ) Weights not initialised"<<endl;
      exit(1);
    }
  
  double d=0;
  
  for (size_t i=0; i<N; i++)
    d+=weights[i]*weights[i];
  
  return nrep*nrep/d;
}
