*     -*-fortran-*-
*
*     Parameters of the LHAPDF evolution grid     
*
      integer nxLHA,nxmLHA,nq2LHA
      double precision q2minLHA,q2maxLHA
      double precision xminLHA,xmLHA,xmaxLHA,Lambda2
*
      parameter(nxLHA    = 5)
      parameter(nxmLHA   = 0)
      parameter(nq2LHA   = 1)
c      parameter(nxLHA    = 6)
c      parameter(nxmLHA   = 3)
c      parameter(nq2LHA   = 2)
      parameter(q2minLHA = 2d0)
      parameter(q2maxLHA = 1d8)
      parameter(xminLHA  = 1d-9)
      parameter(xmLHA    = 1d-1)
      parameter(xmaxLHA  = 1d0)
      parameter(Lambda2  = 0.0625d0)
