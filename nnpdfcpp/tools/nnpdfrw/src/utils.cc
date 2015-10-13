// utils.cc
// PDF utilities

#include <algorithm>
#include <iterator>
#include <iostream>
#include <cmath>

#include "pdfs.h"
#include "utils.h"

// ROOT
#include "TH1.h"
#include "TMatrixDSym.h"

// LHAPDF
#include "LHAPDF/LHAPDF.h"


using namespace std;

void parse_param_input(string infile, rwparam& rwdata) 
{
    // Set up defaut values
    rwdata.ndat=0;
    rwdata.mode=99;
    
    rwdata.outdir="../res/rwreport/";
    rwdata.outdesc="rw";
    rwdata.desc="' NO DESCRIPTION '";
    rwdata.plotform="eps";
    
    rwdata.lhgrid=false;
    rwdata.size=100;
    rwdata.outfile="OUTPDF.LHgrid";
	
    ifstream in (infile.c_str());
	
    if (in.is_open()) {
        string line,name;
		
        while (getline (in, line)) 
        {
            vector<string> res = split(line);
			
            if (res.size()==0)
                continue;
			
            if (res[0]=="PRIOR:")
                if (res.size()==2)
                    rwdata.prior=res[1];
            
            if (res[0]=="RWDATA:")
                if (res.size()==2)
                    rwdata.rwdata=res[1];
            
            if (res[0]=="RWMODE:")
                if (res.size()==2)
                {
                    if (res[1]=="WGT")
                    {  rwdata.mode=0; }
                    else if (res[1]=="CHI2")
                    {  rwdata.mode=1; }
                    else if (res[1]=="DATA")
                    {  rwdata.mode=3; }                    
                }
            
            
            if (res[0]=="NDATA:")
                if (res.size()==2)
                {
                    stringstream cvrt(res[1]);
                    cvrt >> rwdata.ndat;
                }
            
            // OUTPUT
            if (res[0]=="OUTDIR:")
                if (res.size()==2)
                    rwdata.outdir=res[1];
            
            if (res[0]=="OUTDESC:")
                if (res.size()==2)
                    rwdata.outdesc=res[1];
            
            if (res[0]=="PLOTFORMAT:")
                if (res.size()==2)
                    rwdata.plotform=res[1];
			
            // LHGRID Settings
            if (res[0]=="LHGRID:")
                if (res.size()==2)
                    if (res[1]=="Y")
                        rwdata.lhgrid=true;
            
            if (res[0]=="LHGRIDFILE:")
                if (res.size()==2)
                    rwdata.outfile=res[1];
            
            if (res[0]=="LHGRIDNREP:")
                if (res.size()==2)
                {
                    stringstream cvrt(res[1]);
                    cvrt >> rwdata.size;
                }
            
            if (res[0]=="DESCRIPTION:")
            {
                stringstream desc;
                string tmpres;
                size_t n=0;
                while (n<20)
                {
                    if (!getline(in,tmpres))
                    {
                        cout << "ERR: PARAMETER FILE INCORRECTLY FORMATTED"<<endl;
                        exit(1);
                    }
                    
                    if (tmpres=="ENDDESC")
                        break;
                    
                    desc <<" \'"<<tmpres<<"\'"<<endl;
                    n++;
					
                }
                rwdata.desc=desc.str();
            }
        }
    } 
}

TMatrixDSym ReadSymMatrix(const char* filename, size_t n)
{
    TMatrixDSym A(n);
	
    ifstream file(filename);
	
    if (!file.good())
    {
        cout << "ERR: CovMat input file is not good: "<<filename<<endl;
        exit(1);
    }
	
    string line;
    size_t row=0;
    while (getline(file,line))
    {
		
        if (row==n)
        {
            cout << "ERR: Data in CovMat file exceeds length of data vector"<<endl;
            exit(0);
        }
		
        istringstream linestream(line);
		
        double val=0;
        size_t col=0;
        while (linestream>>val)
        {
            if (col==n)
            {
                cout << "ERR: Data in CovMat file exceeds length of data vector"<<endl;
                exit(0);
            }
            A[row][col]=val;
            col++;
        }
        if (col<n)
        {
            cout << "ERR: Not enough data in CovMat"<<endl;
            exit(0);
        }
        row++;		
    }
	
    if (row<n)
    {
        cout << "ERR: Not enough data in CovMat"<<endl;
        exit(0);
    }
	
    return A;
	
}

vector<double> ReadVector(const char* filename, size_t col)
{
    vector<double> out;
    ifstream file(filename);
	
    if (!file.good())
    {
        cout << "ERR: Data file is not good: "<<filename<<endl;
        exit(1);
    }
	
    string line;
    double dummy,val;
    while (getline(file,line))
    {
        stringstream linestream(line);
        for (size_t i=0; i<(col-1); i++)
            if (!(linestream >> dummy))
            {
                cout << "ERR: Not enough columns in data file: " <<filename<<endl;
                exit(1);
            }
		
        if (linestream >> val)
            out.push_back(val);
        else {
            cout << "ERR: Not enough columns in data file: "<<filename<<endl;
            exit(1);
        }
    }
	
    return out;
}

vector<string> split(string& input)
{
    stringstream strstr(input);
    istream_iterator<string> it(strstr);
    istream_iterator<string> end;
    vector<string> results(it, end);
    return results;
}

vector<double> dsplit(string& input)
{
    stringstream strstr(input);
    istream_iterator<double> it(strstr);
    istream_iterator<double> end;
    vector<double> results(it, end);
    return results;
}

double integrate(double data[], size_t npoints, double h)
{
    double integral=0;
    
    integral+=data[0]+data[npoints-1];
    
    for ( size_t j=1; j<(npoints)/2 ; j++ )
        integral+=2*data[2*j -1];
    
    for (size_t j=1; j<(npoints)/2 + 1; j++)
        integral+=4*data[2*j - 2];
    
    return integral*h/3.0;
}

double favg(vector<double> const& data){
    
    double media=0.0;
    
    for (std::size_t i=0; i<data.size(); i++)
        media+=data[i];
    
    return media/((double)data.size());
    
};

void BinLogX(TH1*h) 
{
	
    TAxis *axis = h->GetXaxis(); 
    int bins = axis->GetNbins();
	
    Axis_t from = axis->GetXmin();
    Axis_t to = axis->GetXmax();
    Axis_t width = (to - from) / bins;
    Axis_t *new_bins = new Axis_t[bins + 1];
	
    for (int i = 0; i <= bins; i++) {
        new_bins[i] = pow(10, from + i * width);
    } 
    axis->Set(bins, new_bins); 
    delete new_bins; 
}




