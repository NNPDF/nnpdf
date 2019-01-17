// $Id$
//
// NNPDF++ 2012-2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#pragma once

#include <NNPDF/common.h>
using NNPDF::real;

#include <NNPDF/nnmpi.h>
using NNPDF::MPI;

#include <iostream>
#include <string>

using std::string;
using std::cout;
using std::cerr;
using std::cin;
using std::endl;
using std::ios;
using std::stringstream;

typedef real (*gpdf)(real*);

// ********* Physics **********************

enum {TBAR,BBAR,CBAR,SBAR,UBAR,DBAR,GLUON,D,U,S,C,B,T,PHT};

enum evlnBasis {  EVLN_GAM, EVLN_SNG, EVLN_GLU, EVLN_VAL, EVLN_V3, EVLN_V8, EVLN_V15, EVLN_V24, EVLN_V35,
                  EVLN_T3, EVLN_T8, EVLN_T15, EVLN_T24, EVLN_T35};

// PDF Sum rules
enum sumRule {SUM_MSR, SUM_UVL, SUM_DVL, SUM_SVL, SUM_CVL, SUM_USM, SUM_DSM, SUM_SSM, SUM_CSM};

// ************ Code Organisation Enums ***********

enum minType    {MIN_UNDEF, MIN_GA, MIN_NGA, MIN_NGAFT, MIN_CMAES};
enum stopType   {STOP_UNDEF, STOP_NONE, STOP_LB};
enum paramType  {PARAM_UNDEF, PARAM_NN, PARAM_SLN, PARAM_SLNPP };
enum basisType  {BASIS_UNDEF, BASIS_NN23, BASIS_NN23QED,
                 BASIS_EVOL, BASIS_EVOLQED,BASIS_EVOLS, BASIS_EVOLSQED,
                 BASIS_NN30, BASIS_NN30QED, BASIS_FLVR, BASIS_FLVRQED,
                 BASIS_NN30IC, BASIS_EVOLIC, BASIS_NN31IC, BASIS_LUX, BASIS_NN31ICQED,
                 BASIS_NSR, BASIS_DISEVOL};

enum covType {COV_EXP = false, COV_T0 = true};
enum filterType {DATA_UNFILTERED = false,DATA_FILTERED = true};

namespace Colour {
    enum Code {
        FG_RED      = 31,
        FG_GREEN    = 32,
        FG_YELLOW   = 33,
        FG_BLUE     = 34,
        FG_DEFAULT  = 39,
        BG_RED      = 41,
        BG_GREEN    = 42,
        BG_YELLOW   = 43,
        BG_BLUE     = 44,
        BG_DEFAULT  = 49
    };
    static inline std::ostream& operator<<(std::ostream& os, Code code) {
        return os << "\033[" << static_cast<int>(code) << "m";
    }
}
