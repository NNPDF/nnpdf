// NNPDF14
// Partitioned distance calculation
// n.p.hartland@ed.ac.uk

#include "LHAPDF/LHAPDF.h"

#include "pdfs.h"
#include "obs.h"
#include "pdffuns.h"

#include <iterator>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <algorithm>
#include <list>

using namespace std;

int main(int argc, char **argv)
{  
	//	********************* PARAMETERS ***********************************
  
  // Read configuration filename from arguments
  string SETONE;
  string SETTWO;
  
  if (argc == 3 ) 
  {
    SETONE.assign(argv[1]);
    SETTWO.assign(argv[2]);
  } else {
    cout << "Usage:"<<endl;
    cout << "distances <SETNAME1> <SETNAME2>"<<endl;
    exit(-1);
  }

	// PDF distance output
	ofstream eout("./evdistances.txt");
	ofstream fout("./fldistances.txt");

	// Prior PDF parameters
	const int    SUBSET  = 0;
	const int    NPOINTS = 100;
	const double xmin    = 1e-5;
	const double xmax    = 0.9;
	double Q = sqrt(2.0);
	    
	LHAPDF::initPDFSet(SETONE, LHAPDF::LHGRID, SUBSET);
	
	PDFparams par1;
	par1.nrep = LHAPDF::numberPDF();
	par1.nflav = 13;
	par1.nx=0;
	par1.NAME=SETONE;
	
	PDF one(par1);
	  	
	vector<double> xvals;
	vector< vector<Obs*>* > evlist1;	
  vector< vector<Obs*>* > fllist1;	

	// *********************** SET ONE  *****************************
	double XCH=0.1,x=0;
	
	for (int npx=0; npx<(NPOINTS); npx++) 
	{
    // x values
    if (npx<NPOINTS/2)
      x=(xmin*pow(XCH/xmin,2*(((double) npx ))/((double) NPOINTS )));
    else
      x=(XCH+(xmax-XCH)*(((double) npx+1 ) -(NPOINTS/2+1))/(((double) NPOINTS ) -(NPOINTS/2+1)) );
		
		xvals.push_back(x);	
		
		vector<Obs*>* pelist = obsvect (one);
		vector<Obs*>* pflist = pdfvect (one);

		one.PDFGet(x,Q);

		// Compute pdfs - evln basis
		for (size_t i=0; i<pelist->size(); i++)
			(*pelist)[i]->ComputeObs();
    
    // Compute pdfs - flvr basis
		for (size_t i=0; i<pflist->size(); i++)
			(*pflist)[i]->ComputeObs();
		
		// Push them back
		evlist1.push_back(pelist);
    fllist1.push_back(pflist);

  }

  // ********************** SET TWO ******************************  
	LHAPDF::initPDFSet(SETTWO, LHAPDF::LHGRID, SUBSET);
	
	PDFparams par;
	par.nrep = LHAPDF::numberPDF();
	par.nflav = 13;
	par.nx=0;
	par.NAME=SETTWO;
	
	PDF two(par);
  vector<Obs*>* evlist2 = obsvect(two);
  vector<Obs*>* fllist2 = pdfvect(two);
  
	for (int npx=0; npx<(NPOINTS); npx++) 
	{
		x=xvals[npx];
    two.PDFGet(x,Q);
    
		eout    << x << "  " << Q;
    fout    << x << "  " << Q;

		// Compute all observables at this x, Q
		for (size_t i=0; i<evlist2->size(); i++)
		{
			(*evlist2)[i]->   ComputeObs();

			// For partition uncomment 
			//eout    <<"  "<<sqrt(bootstrap(cvdistance,favg,*(*evlist2)[i],*(*evlist1[npx])[i]) );
			//eout    <<"  "<<sqrt(bootstrap(vardistance,favg,*(*evlist2)[i],*(*evlist1[npx])[i]) );
		
			eout    <<"  "<<sqrt(cvdistance(*(*evlist2)[i],*(*evlist1[npx])[i]) );
			eout    <<"  "<<sqrt(vardistance(*(*evlist2)[i],*(*evlist1[npx])[i]) );
		}
    
    // Compute all observables at this x, Q
		for (size_t i=0; i<fllist2->size(); i++)
		{
			(*fllist2)[i]->   ComputeObs();
			
			// For partition uncomment
			//fout    <<"  "<<sqrt(bootstrap(cvdistance,favg,*(*fllist2)[i],*(*fllist1[npx])[i]) );
			//fout    <<"  "<<sqrt(bootstrap(vardistance,favg,*(*fllist2)[i],*(*fllist1[npx])[i]) );
			
			fout    <<"  "<<sqrt(cvdistance(*(*fllist2)[i],*(*fllist1[npx])[i]) );
			fout    <<"  "<<sqrt(vardistance(*(*fllist2)[i],*(*fllist1[npx])[i]) );
		}
		
		eout  << endl;
    fout  << endl;

		cout  << "Writing x: "<<npx<<"/"<<NPOINTS;
		cout  <<"\n\033[F\033[J";
	}
	
	exit(0);
}

