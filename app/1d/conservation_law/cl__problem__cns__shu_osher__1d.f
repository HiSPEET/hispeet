!> summary:  Shu-Osher test case for the compressible Navier-Stokes equations
!> author:   Benedikt Wex, Joerg Stiller
!> date:     2023/12/14
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CL__Problem__CNS__Shu_Osher__1D

  use Kind_Parameters, only: RNP
  use Constants
  use Execution_Control
  use CL__Problem__CNS__1D
  use CL__Operator__1D

  implicit none
  private

  public :: CL_Problem_CNS_ShuOsher_1D
  public :: CL_Problem_CNS_ShuOsher_Options_1D

  !-----------------------------------------------------------------------------
  !> Type defining the shu osher problem

  type, extends(CL_Problem_CNS_1D) :: CL_Problem_CNS_ShuOsher_1D

    real(RNP) :: x_m    !< initial shock position
    real(RNP) :: p_l    !< pressure left of shock
    real(RNP) :: p_r    !< pressure right of shock
    real(RNP) :: rho_l  !< density right of shock
    real(RNP) :: rho_r  !< density left of shock
    real(RNP) :: rho_w  !< density wave ammplitude
    real(RNP) :: v_l    !< velocity left of shock
    real(RNP) :: v_r    !< velocity right of shock
    real(RNP) :: kappa  !< density wave number

  contains

    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetBoundaryValues

  end type CL_Problem_CNS_ShuOsher_1D

  !-----------------------------------------------------------------------------
  !> shu osher problem options

  type, extends(CL_Problem_CNS_Options_1D) :: CL_Problem_CNS_ShuOsher_Options_1D

    real(RNP) :: x_m   = -4         !< initial shock position
    real(RNP) :: p_l   = 31 * THIRD !< pressure left of shock
    real(RNP) :: p_r   =  1         !< pressure right of shock
    real(RNP) :: rho_l =  3.857143  !< density right of shock
    real(RNP) :: rho_r =  1         !< density left of shock
    real(RNP) :: v_l   =  2.629369  !< velocity left of shock
    real(RNP) :: v_r   =  0         !< velocity right of shock
    real(RNP) :: kappa =  5         !< density wave number
    real(RNP) :: rho_w =  0.2       !< density wave ammplitude

  end type CL_Problem_CNS_ShuOsher_Options_1D

contains

  !---------------------------------------------------------------------------
  !> Initialization of the flow problem

  subroutine SetProblem(this, file)
    class(CL_Problem_CNS_ShuOsher_1D), intent(inout) :: this
    character(len=*), optional, intent(in) :: file !< (*.prm)

    type(CL_Problem_CNS_ShuOsher_Options_1D) :: opt
    namelist/cns_shu_osher_prm/ opt

    logical :: exists, opened
    integer :: prm

    ! preset options ...........................................................

    opt % xb1 = -5
    opt % xb2 =  3 * PI / 2
    opt % bc  = 'S'

    ! check for input file .....................................................

    if (present(file)) then

      inquire(file=trim(file)//'.prm', exist=exists, opened=opened, number=prm)

      if (exists .and. .not. opened) then
        open(newunit=prm, file=trim(file)//'.prm', action='READ')
      else if (opened) then
        rewind(prm)
      else
        call Error( 'SetProblem'                                  &
                  , 'Input file "'//trim(file)//'.prm" not found' &
                  , 'CL__Problem__CNS__Shu_Osher__1D'             )
      end if

      ! read parameters ........................................................

      read(prm, nml=cns_shu_osher_prm)
      if (.not. opened) close(prm)

    end if

    ! set parameters ...........................................................

    ! initialize base type
    call this % Init_CL_Problem_CNS_1D(opt)

    ! specific parameters

    this % x_m   = opt % x_m
    this % p_l   = opt % p_l
    this % p_r   = opt % p_r
    this % rho_l = opt % rho_l
    this % rho_r = opt % rho_r
    this % v_l   = opt % v_l
    this % v_r   = opt % v_r
    this % rho_w = opt % rho_w
    this % kappa = opt % kappa

  end subroutine SetProblem

  !-----------------------------------------------------------------------------
  !> Provides the initial values u(x,0)

  subroutine GetInitialValues(this, cl_operator, u)
    class(CL_Problem_CNS_ShuOsher_1D), intent(in) :: this
    class(CL_Operator_1D), intent(in)  :: cl_operator
    real(RNP), contiguous, intent(out) :: u(0:,:,:) !< u⁰(0:po,1:ne,1:nc)

    real(RNP), dimension(3) :: u_l, u_r
    real(RNP) :: T_l, T_r
    integer   :: e, i

    associate( po => cl_operator % eop % po &
             , ne => cl_operator % ne       &
             , x  => cl_operator % x        )

      T_l = this%p_l / (this%r_gas * this%rho_l)
      call this % PrimitiveToConservative(this%v_l, T_l, this%p_l, u_l)

      T_r = this%p_r / (this%r_gas * this%rho_r)
      call this % PrimitiveToConservative(this%v_r, T_r, this%p_r, u_r)

      do e = 1, ne
      do i = 0, po
        if (x(i,e) < this % x_m) then
          u(i,e,:) = u_l
        else
          u(i,e,:) = u_r
          u(i,e,1) = u(i,e,1) + this%rho_w * sin(this%kappa * x(i,e))
        end if
      end do
      end do

    end associate

  end subroutine GetInitialValues

  !-----------------------------------------------------------------------------
  !> Provides boundary values as Values from inside of the domain

  subroutine GetBoundaryValues(this, t, bv)
    class(CL_Problem_CNS_ShuOsher_1D), intent(in)  :: this
    real(RNP), intent(in)  :: t       !< time
    real(RNP), intent(out) :: bv(:,:) !< boundary values

    real(RNP) :: T_b

    if (any(this % bc /= 'S')) then
      call Error( 'GetBoundaryValues'                 &
                , 'This problem requires constant BC' &
                , 'CL__Problem__CNS__Shu_Osher__1D'   )
    else

      T_b = this%p_l / (this%r_gas * this%rho_l)
      call this % PrimitiveToConservative(this%v_l, T_b, this%p_l, bv(:,1))

      T_b = this%p_r / (this%r_gas * this%rho_r)
      call this % PrimitiveToConservative(this%v_r, T_b, this%p_r, bv(:,2))

      bv(1,2) = bv(1,2) + this%rho_w * sin(this%kappa * this%xb2)

    end if

    ! avoid compiler warning
    if (t > 0) return

  end subroutine GetBoundaryValues

  !=============================================================================

end module CL__Problem__CNS__Shu_Osher__1D
