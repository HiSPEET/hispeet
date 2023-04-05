!> summary:  Interface to scalar conservation problem
!> author:   Joerg Stiller
!> date:     2019/11/23
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================
module CL__Problem__Scalar__1D

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, HALF
  use Execution_Control
  use CL__Problem__1D

  implicit none
  private

  public :: CL_Problem_Scalar_1D

  !-----------------------------------------------------------------------------
  !> Type for defining and handling 1D scalar problems

  type, abstract, extends(CL_Problem_1D) :: CL_Problem_Scalar_1D

    real(RNP)              :: nu_c = 0  !< constant regular  viscosity
    real(RNP), allocatable :: nu_v(:,:) !< variable diffusivity

  contains

    procedure :: DiffusiveFlux
    procedure :: RHS_Convection
    procedure :: RHS_Diffusion
    procedure :: BoundaryValue
    procedure :: BoundaryNormalFlux

    procedure(NumericalConvectiveFlux), deferred :: NumericalConvectiveFlux

  end type CL_Problem_Scalar_1D

  !=============================================================================

  abstract interface

    !---------------------------------------------------------------------------
    !> Numerical flux of scalar conservation problem

    elemental function NumericalConvectiveFlux(problem, ul, ur) result(hc)
      import
      class(CL_Problem_Scalar_1D), intent(in) :: problem
      real(RNP),                   intent(in) :: ul   !< left value
      real(RNP),                   intent(in) :: ur   !< right value
      real(RNP)                               :: hc
    end function NumericalConvectiveFlux

  end interface

contains

  !-----------------------------------------------------------------------------
  !> Diffusive flux of scalar conservation problem

  function DiffusiveFlux(problem, u) result(fd)
    class(CL_Problem_Scalar_1D), intent(in) :: problem
    real(RNP), contiguous,       intent(in) :: u(0:,:,:) !< u(x,t)
    real(RNP) :: fd(0:ubound(u,1), size(u,2), problem%nc)

    fd = 0

    ! constant viscosity .......................................................

    if (problem % nu_c > 0) then
      associate( dx => problem % dx      &
               , nu => problem % nu_c    &
               , Ds => problem % eop % D )

        fd(:,:,1) = nu * 2/dx * matmul(Ds, u(:,:,1))

      end associate
    end if

    ! variable viscosity .......................................................

    if (allocated(problem % nu_v)) then

      ! TBD

    end if

  end function DiffusiveFlux

  !-----------------------------------------------------------------------------
  !> Convective contribution to RHS of DG-SEM formulation

  function RHS_Convection(problem, t, u) result(rc)
    class(CL_Problem_Scalar_1D), intent(in) :: problem
    real(RNP),                   intent(in) :: t         !< time
    real(RNP), contiguous,       intent(in) :: u(0:,:,:) !< u(x,t)
    real(RNP) :: rc(0:ubound(u,1), size(u,2), problem%nc)

    real(RNP), allocatable :: MDt(:,:)  ! Ms (Ds)ᵀ
    real(RNP), allocatable :: fc(:,:)   ! convective fluxes F_c
    real(RNP), allocatable :: ul(:)     ! left  traces u⁻
    real(RNP), allocatable :: ur(:)     ! right traces u⁺
    real(RNP), allocatable :: hc(:)     ! numerical convective fluxes H_c

    integer :: i, j

    associate( po  =>  problem % eop % po  &
             , ne  =>  problem % ne        &
             , eop =>  problem % eop       )

      allocate(MDt(0:po, 0:po), fc(0:po, ne), ul(0:ne), ur(0:ne), hc(0:ne))

      ! interior fluxes ........................................................

      do i = 0, po
      do j = 0, po
        MDt(i,j) = eop%w(j) * eop%D(j,i)
      end do
      end do

      fc = reshape( problem % ConvectiveFlux(u), shape(fc) )
      rc = reshape( matmul(MDt, fc)            , shape(rc) )

      ! numerical fluxes .......................................................

      call GetElementBoundaryTraces(problem, t, u, ul, ur)

      ! local Lax-Friedrichs (LLF) flux
      hc = problem % NumericalConvectiveFlux(ul, ur)

      ! apply numerical fluxes
      rc( 0,:,1) = rc( 0,:,1) + hc(0:ne-1)  ! add to left  side
      rc(po,:,1) = rc(po,:,1) - hc(1:ne  )  ! add to right side

    end associate

  end function RHS_Convection

  !-----------------------------------------------------------------------------
  !> Diffusive contribution to RHS of DG-SEM formulation

  function RHS_Diffusion(problem, t, u) result(rd)
    class(CL_Problem_Scalar_1D), intent(in) :: problem
    real(RNP),                   intent(in) :: t         !< time
    real(RNP), contiguous,       intent(in) :: u(0:,:,:) !< u(x,t)

    real(RNP) :: rd(0:ubound(u,1), size(u,2), problem%nc)

    if (problem % nu_c > 0) then

      call problem % elliptic_op(1) % Apply( dx     = problem % dx   &
                                           , lambda = ZERO           &
                                           , nu     = problem % nu_c &
                                           , u      = u (:,:,1)      &
                                           , r      = rd(:,:,1)      )
      rd = -rd

    else if (allocated(problem % nu_v)) then

      rd = 0

    else

      rd = 0

    end if

  end function RHS_Diffusion

  !-----------------------------------------------------------------------------
  !> Get element-boundary traces

  subroutine GetElementBoundaryTraces(problem, t, u, ul, ur)
    class(CL_Problem_Scalar_1D), intent(in)  :: problem
    real(RNP),                   intent(in)  :: t         !< time
    real(RNP), contiguous,       intent(in)  :: u(0:,:,:) !< u(x,t)
    real(RNP), contiguous,       intent(out) :: ul(0:)    !< u⁻
    real(RNP), contiguous,       intent(out) :: ur(0:)    !< u⁺

    associate( po  =>  problem % eop % po  &
             , ne  =>  problem % ne        )

      ! element interface traces
      ul(1:ne)   = u(po, 1:ne, 1)
      ur(0:ne-1) = u( 0, 1:ne, 1)

      ! left boundary
      select case(problem%bc(1,1))
      case('D')
        ul(0) = 2 * BoundaryValue(problem, 1, t) - u(0,1,1)
      case('P')
        ul(0) = u(po, ne, 1)
      case default
        ul(0) = u( 0,  1, 1)
      end select

      ! right boundary
      select case(problem%bc(2,1))
      case('D')
        ur(ne) = 2 * BoundaryValue(problem, 2, t) - u(po,ne,1)
      case('P')
        ur(ne) = u( 0,  1, 1)
      case default
        ur(ne) = u(po, ne, 1)
      end select

    end associate

  end subroutine GetElementBoundaryTraces

  !-----------------------------------------------------------------------------
  !> Set element-boundary fluxes, reflective at Neumann boundaries

  subroutine SetBoundaryFluxes(problem, ql, qr)
    class(CL_Problem_Scalar_1D), intent(in)    :: problem
    real(RNP), contiguous,       intent(inout) :: ql(0:)   !< q⁻
    real(RNP), contiguous,       intent(inout) :: qr(0:)   !< q⁺

    associate( po  =>  problem % eop % po  &
             , ne  =>  problem % ne        )

      ! left boundary
      select case(problem%bc(1,1))
      case('D')
        ql(0) =  qr(0)
      case('P')
        ql(0) =  ql(ne)
      case default
        ql(0) = -qr(0)
      end select

      ! right boundary
      select case(problem%bc(2,1))
      case('D')
        qr(ne) =  ql(ne)
      case('P')
        qr(ne) =  qr(0)
      case default
        qr(ne) = -ql(ne)
      end select

    end associate

  end subroutine SetBoundaryFluxes

  !-----------------------------------------------------------------------------
  !> Dummy BoundaryValue function -- should not be used

  function BoundaryValue(problem, b, t) result(ub)
    class(CL_Problem_Scalar_1D), intent(in) :: problem
    integer,                     intent(in) :: b    !< boundary {1,2}
    real(RNP),                   intent(in) :: t    !< time
    real(RNP)                               :: ub

    ub = 0

    call Warning('BoundaryValue', 'ub not available', 'Scalar_Problem_1D')

    ! avoid compiler warnings on unused arguments
    if (problem%ne < 0 .or. b < 0 .or. t < 0) return

  end function BoundaryValue

  !-----------------------------------------------------------------------------
  !> Dummy BoundaryNormalFlux function -- should not be used

  function BoundaryNormalFlux(problem, b, t) result(qb)
    class(CL_Problem_Scalar_1D), intent(in) :: problem
    integer,                     intent(in) :: b    !< boundary {1,2}
    real(RNP),                   intent(in) :: t    !< time
    real(RNP)                               :: qb

    qb = 0

    call Warning('BoundaryNormalFlux', 'qb not available', 'Scalar_Problem_1D')

    ! avoid compiler warnings on unused arguments
    if (problem%ne < 0 .or. b < 0 .or. t < 0) return

  end function BoundaryNormalFlux

  !=============================================================================

end module CL__Problem__Scalar__1D

