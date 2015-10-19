// $Id: timer.h 2079 2014-11-12 11:42:33Z s1044006 $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

namespace NNPDF
{

  /**
   * \class Timer
   * \brief Computes the calculation time.
   */
  class Timer {
   private:

    timeval startTime;

   public:

    void start(){
      gettimeofday(&startTime, NULL);
    }

    double stop(){
      timeval endTime;
      long seconds, useconds;
      double duration;

      gettimeofday(&endTime, NULL);

      seconds  = endTime.tv_sec  - startTime.tv_sec;
      useconds = endTime.tv_usec - startTime.tv_usec;

      duration = seconds + useconds/1E6f;

      return duration;
    }

    static void printTime(double duration){
      printf("%5.6f seconds\n", duration);
    }
    
  };

}