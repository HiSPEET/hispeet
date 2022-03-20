!> summary:  Abstract 3D elliptic operator
!> author:   Joerg Stiller
!> date:     2021/08/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @remark
!> TBP for adding BC to the RHS was removed but can be reactivated on demand
!===============================================================================

module DG__Diffusion_Operator__3D
  use Kind_Parameters, only: RNP, RDP, RSP
  use Constants      , only: ZERO, ONE, HALF
  use Array_Assignments                          ! required by
  use Array_Reductions                           ! CG_Method
  use XMPI, only: XMPI_Bcast                     !
  use DG__Element_Operators__1D
  use DG__Schwarz_Operator__3D
  use Spectral_Element_Mesh__3D
  use Spectral_Element_Boundary_Variable__3D

  !-----------------------------------------------------------------------------
  !> Base type for scalar diffusion operators for 3D DG-SEM

  type DG_DiffusionOperator_3D

    class(SpectralElementMesh_3D), pointer :: sem => null()
    real(RNP) :: lambda = 0                  !< Helmholtz parameter
    real(RNP) :: nu_pc  = 0                  !< constant physical diffusivity
    real(RNP) :: nu_sc  = 0                  !< constant spectral diffusivity
    real(RNP), allocatable :: nu_pv(:,:,:,:) !< variable physical diffusivity
    character, allocatable :: bc(:)          !< boundary conditions {P,D,N}
    type(DG_ElementOperators_1D) :: eop
    class(DG_SchwarzOperator_3D), allocatable :: schwarz

contains

    generic   :: Init_DG_DiffusionOperator_3D  =>  Init_C0, Init_CC, Init_V
    generic   :: SetDiffusivity  =>  SetDiffusivity_C, SetDiffusivity_V
    procedure :: Apply
!   procedure :: AddBC
    procedure :: CG_Method

    procedure, private :: Init_C0, Init_CC, Init_V
    procedure, private :: SetDiffusivity_C, SetDiffusivity_V

  end type DG_DiffusionOperator_3D

  ! constructors
  interface DG_DiffusionOperator_3D
    module procedure New_C0
    module procedure New_CC
    module procedure New_V
  end interface

  !=============================================================================
  ! Interfaces to separate module procedures

  interface

    !---------------------------------------------------------------------------
    !> Set constant physical and spectral diffusivities

    module subroutine SetDiffusivity_C(this, nu_p, nu_s)
      class(DG_DiffusionOperator_3D), intent(inout) :: this
      real(RNP),           intent(in) :: nu_p !< physical diffusivity
      real(RNP), optional, intent(in) :: nu_s !< spectral diffusivity [0]
    end subroutine SetDiffusivity_C

    !---------------------------------------------------------------------------
    !> Set variable physical diffusivity

    module subroutine SetDiffusivity_V(this, nu_p)
      class(DG_DiffusionOperator_3D), intent(inout) :: this
      real(RNP), contiguous, intent(in) :: nu_p(:,:,:,:) !< physical diffusivity
    end subroutine SetDiffusivity_V

    !---------------------------------------------------------------------------
    !> Application of the diffusion operator

    module subroutine Apply(this, u, r, f, bv)
      class(DG_DiffusionOperator_3D),  intent(in)  :: this
      real(RNP), contiguous,           intent(in)  :: u(:,:,:,:) !< operand
      real(RNP), contiguous,           intent(out) :: r(:,:,:,:) !< result
      real(RNP), contiguous, optional, intent(in)  :: f(:,:,:,:) !< RHS
      class(SpectralElementBoundaryVariable_3D), optional, intent(in) :: bv(:)
      !< boundary values
    end subroutine Apply

!?!    !---------------------------------------------------------------------------
!?!    !> Addition of BC to the RHS for regular mesh and constant ν
!?!    !>
!?!    !> `bv` is a boundary variable which contains the Dirichlet or Neumann
!?!    !> boundary values four each boundary. These values are applied to the
!?!    !> right hand side `f` according boundary type specified in `this % bc`.
!?!
!?!    module subroutine AddBC_RC(this, bv, f)
!?!      class(DG_DiffusionOperator_3D), intent(in) :: this
!?!      class(SpectralElementBoundaryVariable_3D), target, intent(in) :: bv(:)
!?!      real(RNP), contiguous, intent(inout) :: f(:,:,:,:)
!?!    end subroutine AddBC_RC
!?!
!?!    !-----------------------------------------------------------------------------
!?!    !> Addition of BC to the RHS for regular mesh and variable ν
!?!
!?!    module subroutine AddBC_RV(this, bv, f)
!?!      class(DG_DiffusionOperator_3D), intent(in) :: this
!?!      class(SpectralElementBoundaryVariable_3D), target, intent(in) :: bv(:)
!?!      real(RNP), contiguous, intent(inout) :: f(:,:,:,:)
!?!    end subroutine AddBC_RV
!?!
!?!    !---------------------------------------------------------------------------
!?!    !> Addition of BC to the RHS for deformed mesh and constant ν
!?!
!?!    module subroutine AddBC_DC(this, bv, f)
!?!      class(DG_DiffusionOperator_3D), intent(in) :: this
!?!      class(SpectralElementBoundaryVariable_3D), target, intent(in) :: bv(:)
!?!      real(RNP), contiguous, intent(inout) :: f(:,:,:,:)
!?!    end subroutine AddBC_DC

  end interface

contains

  !=============================================================================
  ! Constructors

  !-----------------------------------------------------------------------------
  !> New diffusion operator with constant physical diffusivity

  function New_C0(sem, dg_opt, lambda, nu_p, bc, schwarz_opt) result(this)
    class(SpectralElementMesh_3D), target, intent(in) :: sem
    class(DG_ElementOptions_1D), intent(in) :: dg_opt
    real(RNP), intent(in) :: lambda !< Helmholtz parameter
    real(RNP), intent(in) :: nu_p   !< physical diffusivity
    character, intent(in) :: bc(:)  !< BC {'D','N','P'}
    class(DG_SchwarzOptions_3D), optional, intent(in) :: schwarz_opt

    type(DG_DiffusionOperator_3D) :: this

    call Init_CC(this, sem, dg_opt, lambda, nu_p, ZERO, bc, schwarz_opt)

  end function New_C0

  !-----------------------------------------------------------------------------
  !> New diffusion operator with constant physical and spectral diffusivities

  function New_CC(sem, dg_opt, lambda, nu_p, nu_s, bc, schwarz_opt) result(this)
    class(SpectralElementMesh_3D), target, intent(in) :: sem
    class(DG_ElementOptions_1D), intent(in) :: dg_opt
    real(RNP), intent(in) :: lambda !< Helmholtz parameter
    real(RNP), intent(in) :: nu_p   !< physical diffusivity
    real(RNP), intent(in) :: nu_s   !< spectral diffusivity [0]
    character, intent(in) :: bc(:)  !< BC {'D','N','P'}
    class(DG_SchwarzOptions_3D), optional, intent(in) :: schwarz_opt

    type(DG_DiffusionOperator_3D) :: this

    call Init_CC(this, sem, dg_opt, lambda, nu_p, nu_s, bc, schwarz_opt)

  end function New_CC

  !-----------------------------------------------------------------------------
  !> New diffusion operator with variable physical diffusivity

  function New_V(sem, dg_opt, lambda, nu_p, bc, schwarz_opt) result(this)
    class(SpectralElementMesh_3D), target, intent(in) :: sem
    class(DG_ElementOptions_1D), intent(in) :: dg_opt
    real(RNP),             intent(in) :: lambda        !< Helmholtz parameter
    real(RNP), contiguous, intent(in) :: nu_p(:,:,:,:) !< variable physical ν
    character,             intent(in) :: bc(:)         !< BC {'D','N','P'}
    class(DG_SchwarzOptions_3D), optional, intent(in) :: schwarz_opt

    type(DG_DiffusionOperator_3D) :: this

    call Init_V(this, sem, dg_opt, lambda, nu_p, bc, schwarz_opt)

  end function New_V

  !=============================================================================
  ! Initialization procedures

  !-----------------------------------------------------------------------------
  !> Initialization with constant physical diffusivity

  subroutine Init_C0(this, sem, dg_opt, lambda, nu_p, bc, schwarz_opt)
    class(DG_DiffusionOperator_3D), intent(inout) :: this
    class(SpectralElementMesh_3D), target, intent(in) :: sem
    class(DG_ElementOptions_1D), intent(in) :: dg_opt
    real(RNP), intent(in) :: lambda !< Helmholtz parameter
    real(RNP), intent(in) :: nu_p   !< physical diffusivity
    character, intent(in) :: bc(:)  !< BC {'D','N','P'}
    class(DG_SchwarzOptions_3D), optional, intent(in) :: schwarz_opt

    call Init_CC(this, sem, dg_opt, lambda, nu_p, ZERO, bc, schwarz_opt)

  end subroutine Init_C0

  !-----------------------------------------------------------------------------
  !> Initialization with constant physical and spectral diffusivities

  subroutine Init_CC(this, sem, dg_opt, lambda, nu_p, nu_s, bc, schwarz_opt)
    class(DG_DiffusionOperator_3D), intent(inout) :: this
    class(SpectralElementMesh_3D), target, intent(in) :: sem
    class(DG_ElementOptions_1D), intent(in) :: dg_opt
    real(RNP), intent(in) :: lambda !< Helmholtz parameter
    real(RNP), intent(in) :: nu_p   !< physical diffusivity
    real(RNP), intent(in) :: nu_s   !< spectral diffusivity [0]
    character, intent(in) :: bc(:)  !< BC {'D','N','P'}
    class(DG_SchwarzOptions_3D), optional, intent(in) :: schwarz_opt

    this % sem     => sem
    this % eop     =  DG_ElementOperators_1D(dg_opt)
    this % lambda  =  lambda
    this % bc      =  bc

    call this % SetDiffusivity(nu_p, nu_s)

    if (present(schwarz_opt)) then
      this % schwarz = DG_SchwarzOperator_3D( schwarz_opt, this%eop, sem%mesh &
                                            , lambda, nu_p, nu_s, bc          )
    end if

  end subroutine Init_CC

  !-----------------------------------------------------------------------------
  !> Initialization with variable physical diffusivity

  subroutine Init_V(this, sem, dg_opt, lambda, nu_p, bc, schwarz_opt)
    class(DG_DiffusionOperator_3D), intent(inout) :: this
    class(SpectralElementMesh_3D), target, intent(in) :: sem
    class(DG_ElementOptions_1D), intent(in) :: dg_opt
    real(RNP), intent(in) :: lambda        !< Helmholtz parameter
    real(RNP), intent(in) :: nu_p(:,:,:,:) !< variable physical diffusivity
    character, intent(in) :: bc(:)         !< BC {'D','N','P'}
    class(DG_SchwarzOptions_3D), optional, intent(in) :: schwarz_opt

    this % sem     => sem
    this % eop     =  DG_ElementOperators_1D(dg_opt)
    this % lambda  =  lambda
    this % bc      =  bc

    call this % SetDiffusivity(nu_p)

    if (present(schwarz_opt)) then
      this % schwarz = DG_SchwarzOperator_3D( schwarz_opt, this%eop, sem%mesh &
                                            , lambda, nu_p, bc                )
    end if

  end subroutine Init_V


  !-----------------------------------------------------------------------------
  !> Conjugate gradient method, modified for r = Au - f

  subroutine CG_Method(this, u, f, bv, i_max, r_red, r_max, ni)

    class(DG_DiffusionOperator_3D), intent(in) :: this
    class(SpectralElementBoundaryVariable_3D), intent(in) :: bv(:) !< BC
    real(RNP),           intent(inout) :: u(:,:,:,:) !< approximate solution
    real(RNP),           intent(in)    :: f(:,:,:,:) !< right hand side
    integer,             intent(in)    :: i_max      !< max num iterations
    real(RNP), optional, intent(in)    :: r_red      !< min residual reduction
    real(RNP), optional, intent(in)    :: r_max      !< max admissible residual
    integer,   optional, intent(out)   :: ni         !< executed num iterations

    contiguous :: u, f

    ! local variables ..........................................................

    real(RNP), dimension(:,:,:,:), allocatable, save :: r, p, q
    real(RNP), save :: rr_term
    logical  , save :: converged

    real(RNP), parameter :: eps = epsilon(ONE) * 1e-3
    real(RNP) :: alpha, pq, rr, rr_old
    logical   :: singular
    integer   :: i

    ! initialization ...........................................................

    associate(mesh => this % sem % mesh)

      ! work space
      !$omp master
      allocate(r, mold = u)
      allocate(p, mold = u)
      allocate(q, mold = u)
      !$omp end master
      !$omp barrier

      singular = abs(this%lambda) < epsilon(ONE) .and. all(this%bc /= 'D')

      ! initial residual .......................................................

      call this % Apply(u, r, f, bv)
      if (singular) then
        call CalibrateArray(r, mesh%comm)
      end if
      call SetArray(p, r)

      rr = ScalarProduct(r, r, mesh%comm)

      !$omp single
      if (present(r_red)) then
        rr_term  = max(ZERO, sqrt(rr) * r_red)**2
        if (present(r_max)) then
          rr_term = max(rr_term, max(ZERO, r_max)**2)
        end if
      else
        rr_term = 0
      end if
      !$omp end single

      rr_old = 0

      ! iteration ..............................................................

      do i = 1, i_max

        ! termination check  . . . . . . . . . . . . . . . . . . . . . . . . . .

        ! MPI master decides about termination
        !$omp master
        if (mesh%part == 0) then
          converged = rr <= rr_term
        end if
        call XMPI_Bcast(converged, root=0, comm=mesh%comm)
        !$omp end master
        !$omp barrier

        if (converged) exit

        rr_old = rr

        ! next iteration . . . . . . . . . . . . . . . . . . . . . . . . . . . .

        ! operator application with no source and homogeneous BC
        call this % Apply(p, q)

        pq = ScalarProduct(p, q, mesh%comm)
!print '(A,G0)', 'rr_old = ', rr_old
!print '(A,G0)', 'pq1 = ', pq
        pq = sign(max(abs(pq),eps), pq)
!print '(A,G0)', 'pq2 = ', pq
        alpha = rr_old / pq
!print '(A,G0)', 'alpha = ', alpha
        call MergeArrays(ONE, u, -alpha, p)
!was!   call MergeArrays(ONE, u, alpha, p)

        if (mod(i,50) == 0) then
          ! compute true residual to get rid of round-off errors
          call this % Apply(u, r, f, bv)
          if (singular) then
            call CalibrateArray(r, mesh%comm)
          end if
        else
          call MergeArrays(ONE, r, -alpha, q)
        end if

        rr = ScalarProduct(r, r, mesh%comm)

        call  MergeArrays(rr/rr_old, p, ONE, r)

      end do

      if (present(ni)) ni = i - 1

      !$omp master
      deallocate(r, p, q)
      !$omp end master

    end associate

  end subroutine CG_Method

!?!  !-----------------------------------------------------------------------------
!?!  !> Addition of boundary conditions to the right hande side
!?!  !>
!?!  !> `bv` is a boundary variable which contains the Dirichlet or Neumann
!?!  !> boundary values four each boundary. These values are applied to the
!?!  !> right hand side `f` according boundary type specified in `this % bc`.
!?!
!?!  subroutine AddBC(this, bv, f)
!?!    class(DG_DiffusionOperator_3D), intent(in) :: this
!?!    class(SpectralElementBoundaryVariable_3D), target, intent(in) :: bv(:)
!?!    real(RNP), contiguous, intent(inout) :: f(:,:,:,:)
!?!
!?!    if (this % sem % mesh % regular) then
!?!      if (allocated(this % nu_pv)) then
!?!        ! regular variable
!?!        call AddBC_RV(this, bv, f)
!?!      else
!?!        ! regular constant
!?!        call AddBC_RC(this, bv, f)
!?!      end if
!?!    else
!?!      ! deformed constant
!?!      call AddBC_DC(this, bv, f)
!?!    end if
!?!
!?!  end subroutine AddBC

  !=============================================================================

end module DG__Diffusion_Operator__3D
