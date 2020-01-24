!> summary:  Interface to scalar conservation problem
!> author:   Joerg Stiller
!> date:     2019/11/23
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Scalar_Problem_1D

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, HALF
  use Execution_Control
  use Conservation_Problem_1D

  implicit none
  private

  public :: ScalarProblem1D

  !-----------------------------------------------------------------------------
  !> Type for defining and handling 1D scalar problems

  type, abstract, extends(ConservationProblem1D) :: ScalarProblem1D

    real(RNP)              :: nu_0r = 0  !< constant regular  viscosity
    real(RNP)              :: nu_0s = 0  !< constant spectral viscosity
    real(RNP), allocatable :: nu_vr(:,:) !< variable regular  viscosity
    real(RNP), allocatable :: nu_vs(:,:) !< variable spectral viscosity

  contains

    procedure :: DiffusiveFlux
    procedure :: RHS_Convection
    procedure :: RHS_Diffusion
    procedure :: BoundaryValue
    procedure :: BoundaryNormalFlux

    procedure(NumericalConvectiveFlux), deferred :: NumericalConvectiveFlux

  end type ScalarProblem1D

  !=============================================================================

  abstract interface

    !---------------------------------------------------------------------------
    !> Numerical flux of scalar conservation problem

    elemental function NumericalConvectiveFlux(problem, ul, ur) result(hc)
      import
      class(ScalarProblem1D), intent(in) :: problem
      real(RNP),              intent(in) :: ul   !< left value
      real(RNP),              intent(in) :: ur   !< right value
      real(RNP)                          :: hc
    end function NumericalConvectiveFlux

  end interface

contains

  !-----------------------------------------------------------------------------
  !> Diffusive flux of scalar conservation problem

  function DiffusiveFlux(problem, u) result(fd)
    class(ScalarProblem1D), intent(in) :: problem
    real(RNP), contiguous,  intent(in) :: u(0:,:,:) !< u(x,t)
    real(RNP) :: fd(0:ubound(u,1), size(u,2), problem%nc)

    fd = 0

    ! constant viscosity .......................................................

    if (problem % nu_0r > 0) then
      associate( dx   => problem % dx      &
               , nu_0r => problem % nu_0r    &
               , Ds   => problem % eop % D )

        fd(:,:,1) = nu_0r * 2/dx * matmul(Ds, u(:,:,1))

      end associate
    end if

    ! variable viscosity .......................................................

    if (allocated(problem % nu_vr)) then

      ! TBD

    end if

    ! spectral viscosity .......................................................

    if (allocated(problem % nu_vs)) then

      ! TBD

    end if

  end function DiffusiveFlux

  !-----------------------------------------------------------------------------
  !> Convective contribution to RHS of DG-SEM formulation

  function RHS_Convection(problem, t, u) result(rc)
    class(ScalarProblem1D), intent(in) :: problem
    real(RNP),              intent(in) :: t         !< time
    real(RNP), contiguous,  intent(in) :: u(0:,:,:) !< u(x,t)
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
    class(ScalarProblem1D), intent(in) :: problem
    real(RNP),              intent(in) :: t         !< time
    real(RNP), contiguous,  intent(in) :: u(0:,:,:) !< u(x,t)

    real(RNP) :: rd(0:ubound(u,1), size(u,2), problem%nc)

    real(RNP), allocatable :: ul(:)     ! left  traces u⁻
    real(RNP), allocatable :: ur(:)     ! right traces u⁺
    real(RNP), allocatable :: ql(:)     ! left  flux   (ν∂u/∂x)⁻
    real(RNP), allocatable :: qr(:)     ! right flux   (ν∂u/∂x)⁺

    real(RNP), allocatable :: jmp_u(:)  ! [u]
    real(RNP), allocatable :: jmp_q(:)  ! [q]
    real(RNP), allocatable :: avg_q(:)  ! {q}

    real(RNP), allocatable :: Ds_i(:,:) ! std diff operator for interior fluxes
    real(RNP), allocatable :: Ds_b(:,:) ! std diff operator for boundary fluxes
    real(RNP), allocatable :: Ls(:,:)   ! std stiffness matrix
    integer   :: e, i

    rd = 0

    if (      problem % nu_0r <= 0             &
        .and. problem % nu_0s <= 0             &
        .and. .not. allocated(problem % nu_vr) &
        .and. .not. allocated(problem % nu_vs) ) return

    associate( po  =>  problem % eop % po  &
             , ne  =>  problem % ne        &
             , eop =>  problem % eop       )

      allocate(ul(0:ne), ur(0:ne), jmp_u(0:ne))
      allocate(ql(0:ne), qr(0:ne), jmp_q(0:ne), avg_q(0:ne))

      allocate(Ls(0:po, 0:po))
      allocate(Ds_i, Ds_b, mold = Ls)

      call GetElementBoundaryTraces(problem, t, u, ul, ur)
      jmp_u = ul - ur

      ! constant viscosity .....................................................

      if (problem % nu_0r > 0) then
        call ConstViscosityContrib(nu = problem % nu_0r, Ds = eop%D, Ls = eop%L)
      end if

      if (problem % nu_0s > 0) then
        call eop % Get_SVV_StandardDiffMatrix(Ds_b)
        call eop % Get_SVV_StandardStiffnessMatrix(Ls)
        call ConstViscosityContrib(nu = problem % nu_0s, Ds = Ds_b, Ls = Ls)
      end if

      if (problem % eop % hybrid) then
        if ((problem % nu_0r > 0) .or. (problem % nu_0s > 0)) then
          if (eop % Has_SVV()) then
            call eop % Get_SVV_StandardDiffMatrix(Ds_b)
            Ds_b = problem % nu_0s * Ds_b
          else
            Ds_b = 0
          end if

          Ds_b = Ds_b + problem % nu_0r * eop%D
          call AddHybridPenaltyTerm(Ds_b)
        end if
      end if

      ! variable viscosity .....................................................

      if (allocated(problem % nu_vr)) then
        ! TBD
      end if

      if (allocated(problem % nu_vs)) then
        ! TBD
      end if

      ! Neumann boundary flux ..................................................

      if (problem % bc(1,1) == 'N') then
        rd( 0, 1,1) = rd( 0, 1,1) - BoundaryNormalFlux(problem, 1, t)
      end if

      if (problem % bc(2,1) == 'N') then
        rd(po,ne,1) = rd(po,ne,1) + BoundaryNormalFlux(problem, 2, t)
      end if

    end associate

  contains

    subroutine ConstViscosityContrib(nu, Ds, Ls)
      real(RNP), intent(in) :: nu          !< constant diffusivity
      real(RNP), intent(in) :: Ds(0:,0:)   !< standard differentiation matrix
      real(RNP), intent(in) :: Ls(0:,0:)   !< standard stiffness matrix

      real(RNP) :: c, g, mu

      associate( po  =>  problem % eop % po  &
               , ne  =>  problem % ne        &
               , dx  =>  problem % dx        )

        mu = problem % eop % PenaltyFactor(dx)
        c  = nu * 2/dx
        g  = c / (4 * mu * nu)

        ! interior
        rd(:,:,1) = rd(:,:,1) - c * matmul(Ls, u(:,:,1))

        ! diffusive flux at element boundaries
        ql(1:ne  ) = c * matmul(Ds(po,:), u(:,:,1))
        qr(0:ne-1) = c * matmul(Ds( 0,:), u(:,:,1))
        call SetBoundaryFluxes(problem, ql, qr)
        avg_q = HALF * (ql + qr)

        ! [q] at interior interfaces, considered only with IP-H
        if (problem % eop % hybrid) then
          if (all(problem % bc(:,1) == 'P')) then
            jmp_q = ql - qr
          else
            jmp_q = [ ZERO, ql(1:ne-1) - qr(1:ne-1), ZERO ]
          end if
        else
          jmp_q = 0
        end if

        ! element boundary fluxes
        do e = 1, ne

          ! {ν∂v/∂x}[u]
          do i = 0, po
            rd(i,e,1) = rd(i,e,1) + HALF * c * ( Ds( 0, i) * jmp_u(e-1) &
                                               + Ds(po, i) * jmp_u(e  ) )
          end do

          ! [v]{ν∂u/∂x}
          rd( 0,e,1) = rd( 0,e,1) - avg_q(e-1)
          rd(po,e,1) = rd(po,e,1) + avg_q(e  )

          ! μ⟨ν⟩[v][u]
          rd( 0,e,1) = rd( 0,e,1) + mu * nu * jmp_u(e-1)
          rd(po,e,1) = rd(po,e,1) - mu * nu * jmp_u(e  )

        end do
      end associate

    end subroutine ConstViscosityContrib

    subroutine AddHybridPenaltyTerm(Bs)
      real(RNP), intent(in) :: Bs(0:,0:) !< standard "flux" operator (ν+νˢQ)D

      real(RNP) :: g, mu

      associate( po     =>  problem % eop % po  &
               , ne     =>  problem % ne        &
               , dx     =>  problem % dx        &
               , nu     =>  problem % nu_0r     &
               , nu_svv =>  problem % nu_0s     )

        mu = problem % eop % PenaltyFactor(dx)
        g  = (2/dx) / (4 * mu * (nu + nu_svv))

        ! diffusive flux at element boundaries
        ql(1:ne  ) = 2/dx * matmul(Bs(po,:), u(:,:,1))
        qr(0:ne-1) = 2/dx * matmul(Bs( 0,:), u(:,:,1))
        call SetBoundaryFluxes(problem, ql, qr)
        avg_q = HALF * (ql + qr)

        ! [q] at interior interfaces
        if (all(problem % bc(:,1) == 'P')) then
          jmp_q = ql - qr
        else
          jmp_q = [ ZERO, ql(1:ne-1) - qr(1:ne-1), ZERO ]
        end if

        ! element boundary fluxes
        do e = 1, ne

          ! 1/4μ⟨ν⟩ [ν∂v/∂x][ν∂u/∂x]
          do i = 0, po
            rd(i,e,1) = rd(i,e,1) - g * Bs( 0,i) * jmp_q(e-1) &
                                  + g * Bs(po,i) * jmp_q(e  )
          end do

        end do
      end associate

    end subroutine AddHybridPenaltyTerm

  end function RHS_Diffusion

  !-----------------------------------------------------------------------------
  !> Get element-boundary traces

  subroutine GetElementBoundaryTraces(problem, t, u, ul, ur)
    class(ScalarProblem1D), intent(in)  :: problem
    real(RNP),              intent(in)  :: t         !< time
    real(RNP), contiguous,  intent(in)  :: u(0:,:,:) !< u(x,t)
    real(RNP), contiguous,  intent(out) :: ul(0:)    !< u⁻
    real(RNP), contiguous,  intent(out) :: ur(0:)    !< u⁺

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
    class(ScalarProblem1D), intent(in)    :: problem
    real(RNP), contiguous,  intent(inout) :: ql(0:)   !< q⁻
    real(RNP), contiguous,  intent(inout) :: qr(0:)   !< q⁺

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
    class(ScalarProblem1D), intent(in) :: problem
    integer,                intent(in) :: b    !< boundary {1,2}
    real(RNP),              intent(in) :: t    !< time
    real(RNP)                          :: ub

    ub = 0

    call Warning('BoundaryValue', 'ub not available', 'Scalar_Problem_1D')

    ! avoid compiler warnings on unused arguments
    if (problem%ne < 0 .or. b < 0 .or. t < 0) return

  end function BoundaryValue

  !-----------------------------------------------------------------------------
  !> Dummy BoundaryNormalFlux function -- should not be used

  function BoundaryNormalFlux(problem, b, t) result(qb)
    class(ScalarProblem1D), intent(in) :: problem
    integer,                intent(in) :: b    !< boundary {1,2}
    real(RNP),              intent(in) :: t    !< time
    real(RNP)                          :: qb

    qb = 0

    call Warning('BoundaryNormalFlux', 'qb not available', 'Scalar_Problem_1D')

    ! avoid compiler warnings on unused arguments
    if (problem%ne < 0 .or. b < 0 .or. t < 0) return

  end function BoundaryNormalFlux

  !=============================================================================

end module Scalar_Problem_1D
