!> summary:  Variable Additive Runge-Kutta methods
!> author:   Joerg Stiller
!> date:     2019/03/21
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Variable Additive Runge-Kutta methods
!>
!> This module provides pairs of additive Runge-Kutta methods consisting of an
!> implicit (EDIRK) part
!>
!>     0      |
!>     c(2)   | a_im(2,1)  a_im(2,3)
!>     :      | :
!>     c(s)   | b_im(1)       ..    b_im(s-1)   b_im(s)
!>     -------+----------------------------------------
!>     1      | b_im(1)       ..    b_im(s-1)   b_im(s)
!>
!>
!> where `c(s) = 1`, and an explicit (ERK) part with `s-1` stages
!>
!>     0      |
!>     c(2)   | a_ex(2,1)
!>     :      | :
!>     c(s-1) | a_ex(s-1,1)   ..    a_ex(s-1,s-2)
!>     -------+-------------------------------------------
!>     1      | b_ex(1)       ..    b_ex(s-2)      b_ex(s)
!>
!> The sub-schemes are synchronized, i.e. their nodes c(i) coincide.
!> In contrast to common RK schemes they accomodate variable intermediate points
!> `c(1:s-1)`. This allows to adjust them to arbitrary point distributions as,
!> e.g., emerging in spectral deferred correction (SDC) methods.
!>
!> #### Application
!>
!> TBD
!>
!===============================================================================

module VA_Runge_Kutta_Method
  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT
  use Kind_Parameters,               only: RNP, RHP
  use Constants,                     only: ZERO, HALF, ONE
  use Execution_Control
  implicit none
  private

  public :: VARK_Method

  !-----------------------------------------------------------------------------
  !> Type for keeping the Butcher tableau of an IMEX Runge-Kutta method

  type VARK_Method
    character(len=80)      :: name  = ' ' !< name of RK method
    integer                :: s     = 0   !< number of stages
    integer                :: order = 0   !< convergence order
    real(RNP), allocatable :: a_im(:,:)   !< EDIRK matrix
    real(RNP), allocatable :: a_ex(:,:)   !< ERK matrix
    real(RNP), allocatable :: b_ex(:)     !< ERK weights
    real(RNP), allocatable :: c(:)        !< RK nodes
  contains
    procedure :: Init_VARK_Method
    procedure :: Write => Write_VARK_Method
  end type VARK_Method

  ! constructor
  interface VARK_Method
    module procedure New_VARK_Method
  end interface

contains

!-------------------------------------------------------------------------------
!> VARK_Method constructor

type(VARK_Method) function New_VARK_Method(t, sh, r) result(this)
  real(RNP),           intent(in) :: t(0:) !< principal points
  integer,             intent(in) :: sh    !< number of helper stages
  real(RNP), optional, intent(in) :: r     !< last coefficient of R(z)

  call Init_VARK_Method(this, t, sh, r)

end function New_VARK_Method

!-------------------------------------------------------------------------------
!> Initialization of VARK_Method

subroutine Init_VARK_Method(this, t, sh, r)
  class(VARK_Method),  intent(inout) :: this
  real(RNP),           intent(in)    :: t(0:) !< principal points
  integer,             intent(in)    :: sh    !< number of helper stages
  real(RNP), optional, intent(in)    :: r     !< last coefficient of R(z)

  integer :: sp, se, si

  sp = size(t)
  si = sp + sh
  se = si - 1

  this % s = si

  if (allocated(this % a_im)) deallocate( this % a_im )
  if (allocated(this % a_ex)) deallocate( this % a_ex )
  if (allocated(this % b_ex)) deallocate( this % b_ex )
  if (allocated(this % c   )) deallocate( this % c    )

  allocate(this % a_im (si, si) )
  allocate(this % a_ex (se, se) )
  allocate(this % b_ex (se)     )
  allocate(this % c    (si)     )

  select case(sp)
  case(4)
    select case(sh)
    case(0)
      this % name  = 'EDIRK_422 + ERK_321'
      this % order = 2
      call EDIRK_422  (t, this % a_im, this % c)
      call ERK_321    (t, this % a_ex, this % b_ex, this % c(1:se), r3 = r)
    case(1)
      this % name  = 'EDIRK_522 + ERK_432'
      this % order = 2
     !call EDIRK_522e (t, this % a_im, this % c)  ! explicit helper
      call EDIRK_522i (t, this % a_im, this % c)  ! implicit helper
      call ERK_432    (t, this % a_ex, this % b_ex, this % c(1:se), r4 = r)
    end select
  end select

end subroutine Init_VARK_Method

!-------------------------------------------------------------------------------
!> EDIRK with 4 stages, order 2 and stage-order 2

subroutine EDIRK_422(t, a, c)
  real(RNP), intent(in)  :: t(0:3) !< principal points
  real(RNP), intent(out) :: a(4,4) !< ERK matrix
  real(RNP), intent(out) :: c(4)   !< ERK nodes

  real(RNP) :: q, w

  a = ZERO
  c = t

  ! CN for stage 2
  a(2,1) = HALF * c(2)
  a(2,2) = c(2) - a(2,1)

  ! BDF2 for stage 3
  w = (c(3) - c(2)) / (c(2) - c(1))
  q = ONE / (1 + 2*w)
  a(3,3) = (c(3) - c(2)) * (1 + w) * q
  a(3,2) = (1 + w)**2 * q * c(2)/2
  a(3,1) = c(3) - a(3,2) - a(3,3)

  ! BDF2 for stage 4
  w = (c(4) - c(3)) / (c(3) - c(2))
  q = ONE / (1 + 2*w)
  a(4,4) = (c(4) - c(3)) * (1 + w) * q
  a(4,3) = (1 + w)**2 * q * a(3,3)
  a(4,2) = (1 + w)**2 * q * a(3,2) - w**2 * q * c(2)/2
  a(4,1) = 1 - (a(4,2) + a(4,3) + a(4,4))

end subroutine EDIRK_422

!-------------------------------------------------------------------------------
!> EDIRK with 5 stages, order 2 and main stage-order 2, implicit helper stage

subroutine EDIRK_522i(t, a, c)
  real(RNP), intent(in)  :: t(0:3) !< principal points
  real(RNP), intent(out) :: a(5,5) !< ERK matrix
  real(RNP), intent(out) :: c(5)   !< ERK nodes

  real(RNP) :: a_(4,4), c_(4)

  a = ZERO
  c = [ t(0:1), t(1:3) ]

  ! coefficients of the underlying base scheme
  call EDIRK_422(t, a_, c_)

  ! backward Euler for stage 2 (helper)
  a(2,2) = c(2)

  ! adopting higher stages from base scheme
  a(3:5,1) = a_(2:4,1)
  a(3:5,3) = a_(2:4,2)
  a(4:5,4) = a_(3:4,3)
  a(  5,5) = a_(  4,4)

end subroutine EDIRK_522i

!-------------------------------------------------------------------------------
!> EDIRK with 5 stages, order 2 and main stage-order 2, explicit helper stage

subroutine EDIRK_522e(t, a, c)
  real(RNP), intent(in)  :: t(0:3) !< principal points
  real(RNP), intent(out) :: a(5,5) !< ERK matrix
  real(RNP), intent(out) :: c(5)   !< ERK nodes

  real(RNP) :: a_(4,4), c_(4)

  a = ZERO
  c = [ t(0:1), t(1:3) ]

  ! coefficients of the underlying base scheme
  call EDIRK_422(t, a_, c_)

  ! forward Euler for stage 2 (helper)
  a(2,1) = c(2)

  ! adopting higher stages from base scheme
  a(3:5,1) = a_(2:4,1)
  a(3:5,3) = a_(2:4,2)
  a(4:5,4) = a_(3:4,3)
  a(  5,5) = a_(  4,4)

end subroutine EDIRK_522e

!-------------------------------------------------------------------------------
!> ERK with 3 stages, order 2 and main stage-order 1 (1,2)

subroutine ERK_321(t, a, b, c, r3)
  real(RNP), intent(in)  :: t(0:3) !< principal points
  real(RNP), intent(out) :: a(3,3) !< ERK matrix
  real(RNP), intent(out) :: b(3)   !< ERK weights
  real(RNP), intent(out) :: c(3)   !< ERK nodes
  real(RNP), intent(in)  :: r3     !< coefficient 3 of stability function
  optional :: r3

  real(RNP) :: r3_

  if (present(r3)) then
    r3_ = r3
  else
   !r3_ = ONE / 6    ! better accuracy
    r3_ = ONE / 15   ! better stability
  end if

  a = ZERO
  b = ZERO
  c = t(0:2)

  a(2,1) = c(2)

  a(3,2) = HALF * c(3)*c(3) / c(2)
  a(3,1) = c(3) - a(3,2)

  b(3) = r3_ / (a(3,2)*c(2))
  b(2) = (HALF - b(3)*c(3)) / c(2)
  b(1) = 1 - (b(2) + b(3))

end subroutine ERK_321

!-------------------------------------------------------------------------------
!> ERK with 4 stages, order 3 and main stage-order 2

subroutine ERK_432(t, a, b, c, r4)
  real(RNP), intent(in)  :: t(0:3) !< interior points
  real(RNP), intent(out) :: a(4,4) !< ERK matrix
  real(RNP), intent(out) :: b(4)   !< ERK weights
  real(RNP), intent(out) :: c(4)   !< ERK nodes
  real(RNP), intent(in)  :: r4     !< coefficient 4 of stability function
  optional :: r4

  real(RNP) :: r4_

  if (present(r4)) then
    r4_ = r4
  else
   !r4_ = ONE / 24   ! better accuracy ?
    r4_ = ONE / 50   ! better stability
  end if

  c = [ t(0:1), t(1:2) ]

  b(4) = (2 - 3*c(3)) / (6 * c(4) * (c(4) - c(3)))
  b(3) = (2 - 3*c(4)) / (6 * c(3) * (c(3) - c(4)))
  b(2) = 0
  b(1) = 1 - (b(2) + b(3) + b(4))

  a = ZERO

  a(2,1) = c(2)

  a(3,2) = HALF * c(3)*c(3) / c(2)
  a(3,1) = c(3) - a(3,2)

  a(4,3) = r4_ / (b(4) * a(3,2) * c(2))
  a(4,2) = (HALF * c(4)*c(4) - a(4,3) * c(3)) / c(2)
  a(4,1) = c(4) - (a(4,2) + a(4,3))

end subroutine ERK_432

!===============================================================================

subroutine Write_VARK_Method(this, unit)
  class(VARK_Method), intent(inout) :: this
  integer,  optional, intent(in)    :: unit  !< output unit

  character(len=*), parameter :: fmt_ca = '(F13.10," |",99F14.10)'
  character(len=*), parameter :: fmt_b =  '(13X,   " |",99F14.10)'
  integer :: i, io

  if (this % s < 1) return

  if (present(unit)) then
    io = unit
  else
    io = OUTPUT_UNIT
  end if

  write(io,'(A,/)')   'Variable additive Runge-Kutta method'
  write(io,'(2A,/)')  'name: ', trim(this % name)
  write(io,'(A,I0)')  'implicit stages = ', size(this % a_im, 1)
  write(io,'(A,I0)')  'explicit stages = ', size(this % a_ex, 1)
  write(io,'(A,I0)')  'order           = ', this % order

  write(io,'(/,A,/)') 'implicit part'
  do i = 1, this%s
    write(io,fmt_ca) this % c(i), this % a_im(i,1:i)
  end do
  write(io,'(A)') repeat('-', 16 + 14*this%s)

  write(io,'(/,A,/)') 'explicit part'
  do i = 1, this%s-1
    write(io,fmt_ca) this % c(i), this % a_ex(i,1:i-1)
  end do
  write(io,'(A)') repeat('-', 1 + 14*this%s)
  write(io,fmt_b) this % b_ex
  write(io,*)

end subroutine Write_VARK_Method

!===============================================================================

end module VA_Runge_Kutta_Method
