// $Id: common.cc 3319 2015-09-25 12:54:27Z s1044006 $
//
// NNPDF++ 2013
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#include "NNPDF/common.h"
#include "config.h"

namespace{

class NullBuffer : public std::streambuf
{
public:
  int overflow(int c) { return c; }
};

std::ostream& get_null_stream()
{
 static NullBuffer nullbuffer;
 static std::ostream nullstream(&nullbuffer);
 return nullstream;
}

int Verbosity = 1;

}

namespace NNPDF
{
	std::string getVersion()
	{
		return std::string(VERSION);
	}

   void SetVerbosity(int level)
   {
     Verbosity = level;
   }

   int GetVerbosity()
   {
     return Verbosity;
   }

   std::ostream& get_logger()
   {
      if (Verbosity){
        return std::cout;
      }
      return get_null_stream();
    }
}
