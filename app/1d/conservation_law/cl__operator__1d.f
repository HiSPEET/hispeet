!> summary:  DG-SEM operator for 1D conservation problems
!> author:   Joerg Stiller
!> date:     2023/04/27
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CL__Operator__1D
  use Kind_Parameters, only: RNP
  use Standard_Operators__1D
  use Embedded_Interpolation__1D
  use DG__Element_Operators__1D
  use DG__Elliptic_Operator__1D
  use DG__Utilities__1D

  implicit none
  private

  public :: CL_Operator_1D
  public :: CL_Operator_Options_1D

  !-----------------------------------------------------------------------------
  !> DG-SEM mesh and operators for 1D conservation problems

  type :: CL_Operator_1D

    integer   :: ne !< number of elements
    integer   :: nc !< number of components
    real(RNP) :: dx !< element length

    type(DG_ElementOperators_1D)   :: eop    !< ops for solution u
    type(StandardOperators_1D)     :: qop    !< ops for quadrature of convection
    type(EmbeddedInterpolation_1D) :: iop_uq !< interpolation from u to q points
    type(DG_EllipticOperator_1D)   :: elliptic_op !< elliptic solvers+smoothers

    real(RNP), allocatable :: x(:,:)  !< mesh points
    real(RNP), allocatable :: Me(:)   !< element mass matrix

    ! element attributes
    integer, allocatable :: activity(:)   !< element activity
    !! - ` 1`  active
    !! - ` 0`  frozen
    !! - `-1`  undefined
    integer, allocatable :: refinement(:) !< element refinement
    !! - ` 1`  regular (active children)
    !! - ` 0`  closure (frozen children)
    !! - `-1`  none
    integer, allocatable :: mark(:)       !< element mark


  contains

    procedure :: Init_CL_Operator_1D

  end type CL_Operator_1D

  ! constructor
  interface CL_Operator_1D
    module procedure New_CL_Operator_1D
  end interface

  !-----------------------------------------------------------------------------
  !> Options for CL_Operator_1D initialization
  !>
  !> The parameters `nc`, `xb1`, `xb2` must be chosen in accordance with the
  !> related problem. It is recommended to extract them from the corresponding
  !> instance of CL_Problem_1D.

  type CL_Operator_Options_1D
    integer                          :: ne = 1  !< number of elements
    integer                          :: nc      !< number of components
    real(RNP)                        :: xb1     !< position of left boundary
    real(RNP)                        :: xb2     !< position of right boundary
    type(DG_ElementOptions_1D)       :: eop     !< DG operator options for v
    type(StandardOperatorOptions_1D) :: qop     !< quadrature for convection
    type(DG_SchwarzOptions_1D)       :: schwarz !< Schwarz preconditioner
  end type CL_Operator_Options_1D

contains

  !-----------------------------------------------------------------------------
  !> Returns a new CL_Operator_1D object

  type(CL_Operator_1D) function New_CL_Operator_1D(opt) result(this)
    class(CL_Operator_Options_1D), intent(in) :: opt

    call Init_CL_Operator_1D(this, opt)

  end function New_CL_Operator_1D

  !-----------------------------------------------------------------------------
  !> Initialization of CL_Operator_1D

  subroutine Init_CL_Operator_1D(this, opt)
    class(CL_Operator_1D),         intent(inout) :: this
    class(CL_Operator_Options_1D), intent(in)    :: opt

    ! free allocated components ................................................

    if (allocated(this % x         ))  deallocate(this % x         )
    if (allocated(this % Me        ))  deallocate(this % Me        )
    if (allocated(this % activity  ))  deallocate(this % activity  )
    if (allocated(this % refinement))  deallocate(this % refinement)
    if (allocated(this % mark      ))  deallocate(this % mark      )

    ! basic initialization .....................................................

    this % ne = opt % ne
    this % nc = opt % nc
    this % dx = (opt % xb2 - opt % xb1) / this % ne

    ! element operators ........................................................

    this % eop    = DG_ElementOperators_1D(opt % eop)
    this % qop    = StandardOperators_1D(opt % qop)
    this % iop_uq = EmbeddedInterpolation_1D(this % eop, this % qop % x)

    this % elliptic_op = DG_EllipticOperator_1D(opt % eop, opt % schwarz)

    ! mesh and derived operators ...............................................

    associate(ne => this % ne, po => this % eop % po)

      ! mesh points
      allocate(this % x(0:po,ne))
      call DG_GetMeshPoints_1D(this%eop, opt%xb1, opt%xb2, this%dx, this%x)

      ! element mass matrix
      allocate(this % Me(0:po), source = this%dx/2 * this%eop%w)

      ! attributes
      allocate(this % activity  (ne), source =  1) ! default: active
      allocate(this % refinement(ne), source = -1) ! default: no refinement
      allocate(this % mark      (ne), source = -1) ! default: no mark

    end associate

  end subroutine Init_CL_Operator_1D

  !=============================================================================

end module CL__Operator__1D
