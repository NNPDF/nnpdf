// $Id$
//
// NNPDF++ 2015
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it

#pragma once

#include "common.h"

namespace NNPDF
{
  /**
   * @brief The MPI class: handler for openMPI world
   */
  class MPI
  {
  public:
    static void Init();
    static void Finalize();
    static int  TaskID() { return getInstance()->ftaskid; }
    static int  NumTasks() { return getInstance()->fnumbertasks; }
    static void DecomposeDomain(int domain_size, int taskid,
                                int* subdomain_start,
                                int* subdomain_size);
    static void SendJobs(int, int, int, real*, real*, real*);
    static void RecvJobs(int, int, real *);
    static void ResetWorld();

  private:
    MPI();
    ~MPI();
    void SlaveConvoluteLoop();

    static MPI* getInstance()
    {
      if (!mpiInstance)
        mpiInstance = new MPI();
      return mpiInstance;
    }

    // Members
    int fnumbertasks;
    int ftaskid;
    int fworld;
    static MPI* mpiInstance;
  };
}
