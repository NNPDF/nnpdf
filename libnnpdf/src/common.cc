// $Id: common.cc 3319 2015-09-25 12:54:27Z s1044006 $
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#include "NNPDF/common.h"
#include "config.h"

namespace NNPDF
{
	std::string getVersion()
	{
		return std::string(VERSION);
	}
}