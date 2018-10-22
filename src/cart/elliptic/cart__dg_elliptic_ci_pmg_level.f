!> summary:  Polynomial multigrid level for use with DG elliptic solvers
!> author:   Joerg Stiller
!> date:     2018/02/01
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Polynomial multigrid level for use with DG elliptic solvers
!===============================================================================

module CART__DG_Elliptic_CI_PMG_Level

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO
  use Array_Assignments, only: AssignScalar

  use CART__DG_Element_Operators
  use CART__DG_Elliptic_CI_Schwarz

  implicit none
  private

  public :: PMG_Level

  !-----------------------------------------------------------------------------
  !> Polynomial multigrid level

  type PMG_Level

    integer :: po = -1  !< polynomial order
    integer :: ne = -1  !< number of local elements
    logical :: bottom   !< switch for bottom level
    logical :: top      !< switch for top level
    integer :: ns1 = 1  !< number of pre-smoothing steps
    integer :: ns2 = 1  !< number of post-smoothing steps

    type(DG_ElementOperators3D) :: eop      !< 1D DG element operators
    type(SchwarzOperator)       :: schwarz  !< Schwarz operator

    logical                :: has_var_nu    !< switch to variable diffusivity
    real(RNP)              :: lambda        !< Helmholtz parameter
    real(RNP)              :: nu_c          !< constant diffusivity
    real(RNP), allocatable :: nu_v(:,:,:,:) !< variable diffusivity
    real(RNP), allocatable :: nu_e(:)       !< element-mean diffusivity
    character, allocatable :: bc(:)         !< boundary conditions {D,N}

    real(RNP), allocatable :: u(:,:,:,:)    !< solution
    real(RNP), allocatable :: f(:,:,:,:)    !< RHS
    real(RNP), allocatable :: v(:,:,:,:)    !< residual / correction

  contains

    procedure :: New => New_PMG_Level

    generic   :: SetProblem => SetProblem_C
    procedure, private :: SetProblem_C

    procedure :: FreeWorkspace
    final     :: Delete_PMG_Level

  end type PMG_Level

contains

!===============================================================================
! type-bound procedures

!-------------------------------------------------------------------------------
!> Initialization of a PMG_Level object

subroutine New_PMG_Level( this                      &
                        , po, bottom, top           &
                        , ne, dx, penalty           &
                        , ns1, ns2                  &
                        , delta, no_min, weighting  &
                        )

  ! arguments ..................................................................

  class(PMG_Level), intent(inout) :: this

  integer,           intent(in) :: po        !< polynomial order
  logical,           intent(in) :: bottom    !< switch for bottom level
  logical,           intent(in) :: top       !< switch for top level
  integer,           intent(in) :: ne        !< number of local elements
  real(RNP),         intent(in) :: dx(3)     !< element extensions
  real(RNP),         intent(in) :: penalty   !< penalty parameter (> 1)
  integer,           intent(in) :: ns1       !< number of pre-smoothing steps
  integer,           intent(in) :: ns2       !< number of post-smoothing steps
  real(RNP),         intent(in) :: delta(3)  !< relative element overlap
  integer, optional, intent(in) :: no_min    !< min overlap in points [0]
  integer, optional, intent(in) :: weighting !< weighting method [5]

  ! clean-up ...................................................................

  call Delete_PMG_Level(this)

  ! initialization of components ...............................................

  this % bottom = bottom
  this % top    = top
  this % po     = po
  this % ne     = ne
  this % ns1    = ns1
  this % ns2    = ns2

  call this % eop     % New(po, dx, penalty)
  call this % schwarz % New(this%eop, delta, no_min, weighting)

end subroutine New_PMG_Level

!-------------------------------------------------------------------------------
!> (Re)initialize problem with constant diffusivity

subroutine SetProblem_C(this, lambda, nu, bc)
  class(PMG_Level), intent(inout) :: this
  real(RNP), intent(in) :: lambda    !< Helmholtz parameter
  real(RNP), intent(in) :: nu        !< constant diffusivity
  character, intent(in) :: bc(:)     !< boundary conditions {P,D,N}

  integer :: np

  associate(po => this%po, ne => this%ne)

    this % has_var_nu = .false.
    this % lambda     =  lambda
    this % nu_c       =  nu
    this % bc         =  bc

    np = po + 1

    if (allocated(this % nu_e)) deallocate(this % nu_e)
    if (allocated(this % nu_v)) deallocate(this % nu_v)

    if (allocated(this % u)) then
      if (any(shape(this % u) /= [np, np, np, ne])) then
        deallocate(this % u)
        deallocate(this % f)
        deallocate(this % v)
      end if
    end if

    if (.not. allocated(this % u   )) allocate(this % u(0:po, 0:po, 0:po, ne))
    if (.not. allocated(this % f   )) allocate(this % f(0:po, 0:po, 0:po, ne))
    if (.not. allocated(this % v   )) allocate(this % v(0:po, 0:po, 0:po, ne))

    call AssignScalar(this % u, ZERO)
    call AssignScalar(this % f, ZERO)
    call AssignScalar(this % v, ZERO)

  end associate

end subroutine SetProblem_C

!-------------------------------------------------------------------------------
!> Deallocate work arrays

subroutine FreeWorkspace(this)
  class(PMG_Level), intent(inout) :: this

  if (allocated( this % nu_v )) deallocate(this % nu_v)
  if (allocated( this % nu_e )) deallocate(this % nu_e)
  if (allocated( this % bc   )) deallocate(this % bc  )
  if (allocated( this % u    )) deallocate(this % u   )
  if (allocated( this % f    )) deallocate(this % f   )
  if (allocated( this % v    )) deallocate(this % v   )

end subroutine FreeWorkspace

!-------------------------------------------------------------------------------
!> Finalization of a PMG_Level object

subroutine Delete_PMG_Level(this)
  type(PMG_Level), intent(inout) :: this

  call this % FreeWorkspace()

end subroutine Delete_PMG_Level

!===============================================================================

end module CART__DG_Elliptic_CI_PMG_Level
