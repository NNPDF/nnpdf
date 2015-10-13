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

  int const np=20;
  double fonlla_fk[np]={0};
  double fonlla_bench[np]={0};
  double x[np];
  int idum;

  ifstream out1;
  out1.open("fkcheck-fonllb.res");
  if(out1.fail()){
    cout<<"Unable to open unit"<<endl;
    exit(-10);
  }
  for(int ip=0; ip<np; ip++){
    out1>>idum>>fonlla_fk[ip];
  }
  out1.close();

  double adum=0;
  ifstream out2;
  out2.open("bench-fonllb.dat");
  if(out2.fail()){
    cout<<"Unable to open unit"<<endl;
    exit(-10);
  }
  cout<<" x[ip]    fonlla_bench[ip]   fonlla_fk[ip]   diff(%)"<<endl;
  for(int ip=0; ip<np; ip++){
    out2>>x[ip]>>adum>>fonlla_bench[ip]>>adum;
    //cout<<adum<<endl;
    double diff=1e2*fabs((fonlla_bench[ip] - fonlla_fk[ip] )/fonlla_fk[ip]); 
    cout<<x[ip]<<" "<<fonlla_bench[ip]<<" "<<fonlla_fk[ip]<<" "<<diff<<endl;

  }
  out2.close();

  
  return 0;
}
