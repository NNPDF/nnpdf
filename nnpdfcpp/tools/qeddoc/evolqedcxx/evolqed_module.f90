! Module to access the C++ evolqed class
! from a fortran program.
!
! Availables subroutines:
!   new(EvolQED_type, integer, integer) ! Initialize
!   delete(EvolQED_type)                ! Destroys obj.
!
! Available functions:
!   GetNF(EvolQED_type)      ! Get the current nf.
!   GetQ20Ref(EvolQED_type)  ! Get the current Q20 ref.
!
!
! Usage, in the program main call:
!
! use EvolQED_module
! type(EvolQED_type) evol
!
! call new(evol, nf, run)
! write(6,!) GetNF(evol)
! call delete(evol)
!

module EvolQED_module  
  use, intrinsic :: ISO_C_Binding, only: C_int, C_ptr, C_NULL_ptr, C_double, C_bool
  implicit none
  
  private
  type EvolQED_type
     private
     type(C_ptr) :: object = C_NULL_ptr
  end type EvolQED_type
  
  interface 
     ! Constructor
     function C_EvolQED__new (nf, run) result(this) bind(C, name="EvolQED__new")
       import
       type(C_ptr) :: this
       integer(C_int), value :: nf
       logical(C_bool), value :: run
     end function C_EvolQED__new

     ! Destructor
     subroutine C_EvolQED__delete(this) bind(C, name="EvolQED__delete")
       import
       type(C_ptr), value :: this
     end subroutine C_EvolQED__delete

!       ! The most important function
!     subroutine C_EvolQED__EvolFactNQED(this, N, Q2I, Q2F, EFNNS, EFNSG) bind(C, name="EvolQED__EvolFactNQED")
!       import
!       type(C_ptr), value :: this
!       real(C_double), value :: N, Q2I, Q2F
!       type(C_ptr), value :: EFNNS, EFNSG
!     end subroutine C_EvolQED__EvolFactNQED

     ! Set Methods
     subroutine C_EvolQED__SetNF(this, nf) bind(C, name="EvolQED__SetNF")
       import
       type(C_ptr), value :: this
       integer(C_int), value :: nf
     end subroutine C_EvolQED__SetNF

     subroutine C_EvolQED__ActivateRunning(this, run) bind(C, name="EvolQED__ActivateRunning")
       import
       type(C_ptr), value :: this
       logical(C_bool), value :: run
     end subroutine C_EvolQED__ActivateRunning

     subroutine C_EvolQED__SetRefCoupling(this, Q20, AlphaQ20) bind(C, name="EvolQED__SetRefCoupling")
       import
       type(C_ptr), value :: this
       real(C_double), value :: Q20, AlphaQ20
     end subroutine C_EvolQED__SetRefCoupling

     ! GetNF
     function C_EvolQED__GetNF(this) result(nf) bind(C, name="EvolQED__GetNF")
       import
       integer(C_int) :: nf
       type(C_ptr), value :: this
     end function C_EvolQED__GetNF

     function C_EvolQED__GetQ20Ref(this) result(q2ref) bind(C, name="EvolQED__GetQ20Ref")
       import
       real(C_double) :: q2ref
       type(C_ptr), value :: this
     end function C_EvolQED__GetQ20Ref

     function C_EvolQED__GetAlphaQ20(this) result(alpha) bind(C, name="EvolQED__GetAlphaQ20")
       import
       real(C_double) :: alpha
       type(C_ptr), value :: this
     end function C_EvolQED__GetAlphaQ20

  end interface
  
  interface new
     module procedure EvolQED__new
  end interface new

  interface delete
     module procedure EvolQED__delete
  end interface delete

!  interface EvolFactNQED
!     module procedure EvolQED__EvolFactNQED
!  end interface EvolFactNQED

  interface SetNF
     module procedure EvolQED__SetNF
  end interface SetNF

  interface ActivateRunning
     module procedure EvolQED__ActivateRunning
  end interface ActivateRunning

  interface SetRefCoupling
     module procedure EvolQED__SetRefCoupling
  end interface SetRefCoupling

  interface GetNF
     module procedure EvolQED__GetNF
  end interface GetNF

  interface GetQ20Ref
     module procedure EvolQED__GetQ20Ref
  end interface GetQ20Ref

  interface GetAlphaQ20
     module procedure EvolQED__GetAlphaQ20
  end interface GetAlphaQ20

  public :: new, delete, EvolQED_type, SetNF, GetNF, GetQ20Ref, GetAlphaQ20, ActivateRunning!, EvolFactNQED

contains
  subroutine EvolQED__new(this, nf, run)
    type(EvolQED_type), intent(out) :: this
    integer :: nf
    logical :: run
    this%object = C_EvolQED__new(int(nf,C_int),logical(run,C_bool))
  end subroutine EvolQED__new

  subroutine EvolQED__delete(this)
    type(EvolQED_type), intent(inout) :: this
    call C_EvolQED__delete(this%object)
    this%object = C_NULL_ptr
  end subroutine EvolQED__delete

!  subroutine EvolQED__EvolFactNQED(this, N, Q2I, Q2F, EFNNS, EFNSG)       
!    type(EvolQED_type) :: this
!    type(C_ptr), value :: EFNNS, EFNSG
!    real(C_double), value :: N, Q2I, Q2F
!    call C_EvolQED__EvolFactNQED(this%object,real(N,C_double),real(Q2I,C_double),&
!         & real(Q2F,C_double), real(EFNNS,C_double), real(EFNSG, C_double))
!  end subroutine EvolQED__EvolFactNQED

  subroutine EvolQED__SetNF(this, nf)
    type(EvolQED_type) :: this
    integer :: nf
    call C_EvolQED__SetNF(this%object,int(nf,C_int))
  end subroutine EvolQED__SetNF

  subroutine EvolQED__ActivateRunning(this, run)
    type(EvolQED_type) :: this
    logical :: run
    call C_EvolQED__ActivateRunning(this%object,logical(run,C_bool))
  end subroutine EvolQED__ActivateRunning

  subroutine EvolQED__SetRefCoupling(this, Q20, AlphaQ20)
    type(EvolQED_type) :: this
    real*8 Q20, AlphaQ20
    call C_EvolQED__SetRefCoupling(this%object,real(Q20,C_double),real(AlphaQ20,C_double))
  end subroutine EvolQED__SetRefCoupling

  function EvolQED__GetNF(this) result(nf)
    type(EvolQED_type), intent(in) :: this
    integer nf
    nf = C_EvolQED__GetNF(this%object)
  end function EvolQED__GetNF

  function EvolQED__GetQ20Ref(this) result(q2ref)
    type(EvolQED_type), intent(in) :: this
    real*8 q2ref
    q2ref = C_EvolQED__GetQ20Ref(this%object)
  end function EvolQED__GetQ20Ref

  function EvolQED__GetAlphaQ20(this) result(alpha)
    type(EvolQED_type), intent(in) :: this
    real*8 alpha
    alpha = C_EvolQED__GetAlphaQ20(this%object)
  end function EvolQED__GetAlphaQ20

end module EvolQED_module
