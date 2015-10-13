!     Example how to use the fortram module
!     The NNPDF++ Collaboration - 2012
!     Author: Stefano Carrazza, stefano.carrazza@mi.infn.it

      program main

      use EvolQED_module
      type(EvolQED_type) evol

      integer nf, N
      logical run
      double precision Q2I, Q2F
      double precision SOLNS(10)
      double precision SOLSG(3,3)
      
      nf = 3  ! 3 flavours
      run = .true. ! activate alpha running QED
      Q2I = 2.0
      Q2F = 10000.0

      write(6,*) "Nf = ", nf
      write(6,*) "Running alphaQED: ", run
      write(6,*) ""
      
      call new(evol, nf, run)

      do N=2,2
         write(6,*) "Mellin variable N = ", N
         write(6,*) "Q2I = ", Q2I, "GeV^2"
         write(6,*) "Q2F = ", Q2F, "GeV^2"

         ! EvolFactNQED(N, Q2I, Q2F, SOLNS, SOLSG)
      enddo

      call delete(evol)
      
      end program main
