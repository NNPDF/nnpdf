// $Id
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cmath>

using namespace std;

int main()
{  

  cout<<"\n -------------------------------- \n"<<endl;
  cout<<" Les Houches heavy quark benchmarks "<<endl;
  cout<<"\n -------------------------------- \n"<<endl;

  int const np=120;
  //  int const np=15;
  double fk[np]={0};
  double bench[np]={0};
  double x[np];
  double Q2[np];
  double y[np];
  int idum;

  ifstream out1;
  out1.open("fkcheck.res");
  if(out1.fail()){
    cout<<"Unable to open unit"<<endl;
    exit(-10);
  }
  for(int ip=0; ip<np; ip++){
    out1>>idum>>fk[ip];
  }
  out1.close();

  double adum=0;
  ifstream out2;
  out2.open("bench.dat");
  if(out2.fail()){
    cout<<"Unable to open unit"<<endl;
    exit(-10);
  }
  cout<<" x[ip]    bench[ip]   fk[ip]   diff(%)"<<endl;
  for(int ip=0; ip<np; ip++){
    out2>>idum>>x[ip]>>Q2[ip]>>y[ip]>>bench[ip];
    //cout<<adum<<endl;
    double diff=1e2*fabs((bench[ip] - fk[ip])/fk[ip]); 
    //    cout<<x[ip]<<" "<<bench[ip]<<" "<<fk[ip]<<" "<<diff<<endl;
    cout<<ip<<" "<<diff<<endl;
    //    cout<<ip<<" "<<x[ip]<<" "<<Q2[ip]<<" "<<diff<<endl;
  }
  out2.close();

  
  return 0;
}
