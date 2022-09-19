!> summary:  Abstract 3D elliptic operator
!> author:   Joerg Stiller
!> date:     2021/08/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @remark
!> TBP for adding BC to the RHS was removed but can be reactivated on demand
!===============================================================================

module DG__Elliptic_Operator__3D
  use Kind_Parameters, only: RNP, RDP, RSP
  use Constants      , only: ZERO, ONE, HALF
  use Execution_Control
  use Array_Assignments                          ! required by
  use Array_Reductions                           ! CG_Method
  use TPO__Schwarz__3D_CI
  use XMPI, only: XMPI_Bcast                     !
  use DG__Element_Operators__1D
  use DG__Schwarz_Operator__3D
  use Element_Transfer_Buffer__3D
  use Spectral_Element_Mesh__3D
  use Spectral_Element_Boundary_Variable__3D

  !-----------------------------------------------------------------------------
  !> Base type for scalar diffusion operators for 3D DG-SEM

  type DG_EllipticOperator_3D

    class(SpectralElementMesh_3D), pointer :: sem => null()
    real(RNP) :: lambda = 0                  !< Helmholtz parameter
    real(RNP) :: nu_pc  = 0                  !< constant physical diffusivity
    real(RNP) :: nu_sc  = 0                  !< constant spectral diffusivity
    real(RNP), allocatable :: nu_pv(:,:,:,:) !< variable physical diffusivity
    character, allocatable :: bc(:)          !< boundary conditions {P,D,N}
    type(DG_ElementOperators_1D) :: eop
    class(DG_SchwarzOperator_3D), allocatable :: schwarz

contains

    generic   :: Init_DG_EllipticOperator_3D  =>  Init_C0, Init_CC, Init_V
    generic   :: SetDiffusivity  =>  SetDiffusivity_C, SetDiffusivity_V
    procedure :: Apply
!   procedure :: AddBC
    procedure :: CG_Method
    procedure :: Schwarz_Method
    procedure :: SchwarzPCG_Method

    procedure, private :: Init_C0, Init_CC, Init_V
    procedure, private :: SetDiffusivity_C, SetDiffusivity_V

  end type DG_EllipticOperator_3D

  ! constructors
  interface DG_EllipticOperator_3D
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
      class(DG_EllipticOperator_3D), intent(inout) :: this
      real(RNP),           intent(in) :: nu_p !< physical diffusivity
      real(RNP), optional, intent(in) :: nu_s !< spectral diffusivity [0]
    end subroutine SetDiffusivity_C

    !---------------------------------------------------------------------------
    !> Set variable physical diffusivity

    module subroutine SetDiffusivity_V(this, nu_p)
      class(DG_EllipticOperator_3D), intent(inout) :: this
      real(RNP), contiguous, intent(in) :: nu_p(:,:,:,:) !< physical diffusivity
    end subroutine SetDiffusivity_V

    !---------------------------------------------------------------------------
    !> Application of the diffusion operator

    module subroutine Apply(this, u, r, f, bv)
      class(DG_EllipticOperator_3D),  intent(in)  :: this
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
!?!      class(DG_EllipticOperator_3D), intent(in) :: this
!?!      class(SpectralElementBoundaryVariable_3D), target, intent(in) :: bv(:)
!?!      real(RNP), contiguous, intent(inout) :: f(:,:,:,:)
!?!    end subroutine AddBC_RC
!?!
!?!    !-----------------------------------------------------------------------------
!?!    !> Addition of BC to the RHS for regular mesh and variable ν
!?!
!?!    module subroutine AddBC_RV(this, bv, f)
!?!      class(DG_EllipticOperator_3D), intent(in) :: this
!?!      class(SpectralElementBoundaryVariable_3D), target, intent(in) :: bv(:)
!?!      real(RNP), contiguous, intent(inout) :: f(:,:,:,:)
!?!    end subroutine AddBC_RV
!?!
!?!    !---------------------------------------------------------------------------
!?!    !> Addition of BC to the RHS for deformed mesh and constant ν
!?!
!?!    module subroutine AddBC_DC(this, bv, f)
!?!      class(DG_EllipticOperator_3D), intent(in) :: this
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

    type(DG_EllipticOperator_3D) :: this

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

    type(DG_EllipticOperator_3D) :: this

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

    type(DG_EllipticOperator_3D) :: this

    call Init_V(this, sem, dg_opt, lambda, nu_p, bc, schwarz_opt)

  end function New_V

  !=============================================================================
  ! Initialization procedures

  !-----------------------------------------------------------------------------
  !> Initialization with constant physical diffusivity

  subroutine Init_C0(this, sem, dg_opt, lambda, nu_p, bc, schwarz_opt)
    class(DG_EllipticOperator_3D), intent(inout) :: this
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
    class(DG_EllipticOperator_3D), intent(inout) :: this
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
    class(DG_EllipticOperator_3D), intent(inout) :: this
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

    class(DG_EllipticOperator_3D), intent(in) :: this
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

    ! skip empty partition
    if (this % sem % mesh % part < 0) return

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
        call CalibrateArray(r, mesh%comm_parts)
      end if
      call SetArray(p, r)

      rr = ScalarProduct(r, r, mesh%comm_parts)

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
        call XMPI_Bcast(converged, root=0, comm=mesh%comm_parts)
        !$omp end master
        !$omp barrier

        if (converged) exit

        rr_old = rr

        ! next iteration . . . . . . . . . . . . . . . . . . . . . . . . . . . .

        ! operator application with no source and homogeneous BC
        call this % Apply(p, q)

        pq = ScalarProduct(p, q, mesh%comm_parts)
        pq = sign(max(abs(pq),eps), pq)
        alpha = rr_old / pq
        call MergeArrays(ONE, u, -alpha, p)

        if (mod(i,50) == 0) then
          ! compute true residual to get rid of round-off errors
          call this % Apply(u, r, f, bv)
          if (singular) then
            call CalibrateArray(r, mesh%comm_parts)
          end if
        else
          call MergeArrays(ONE, r, -alpha, q)
        end if

        rr = ScalarProduct(r, r, mesh%comm_parts)

        call  MergeArrays(rr/rr_old, p, ONE, r)

      end do

      if (present(ni)) ni = i - 1

      !$omp master
      deallocate(r, p, q)
      !$omp end master

    end associate

  end subroutine CG_Method

  !-----------------------------------------------------------------------------
  !> Element-centered overlapping Schwarz method

  subroutine Schwarz_Method(this, u, f, bv, i_max, r_red, r_max, standby, ni)

    class(DG_EllipticOperator_3D), intent(in) :: this
    real(RNP),           intent(inout) :: u(:,:,:,:) !< approximate solution
    real(RNP),           intent(in)    :: f(:,:,:,:) !< right hand side

    class(SpectralElementBoundaryVariable_3D), optional, intent(in) :: bv(:)
    !< boundary values [homogeneous]

    integer,             intent(in)    :: i_max      !< max num iterations
    real(RNP), optional, intent(in)    :: r_red      !< min residual reduction
    real(RNP), optional, intent(in)    :: r_max      !< max admissible residual
    logical,   optional, intent(in)    :: standby    !< reuse workspace [F]
    integer,   optional, intent(out)   :: ni         !< executed num iterations

    contiguous :: u, f

    ! local variables ..........................................................

    real(RNP), allocatable, save :: r(:,:,:,:)     ! residual, including ghosts

    real(RDP), allocatable, save :: us_dp(:,:,:,:) ! solution of subsystems
    real(RDP), allocatable, save :: fs_dp(:,:,:,:) ! RHS of subsystems

    real(RSP), allocatable, save :: us_sp(:,:,:,:) ! solution of subsystems
    real(RSP), allocatable, save :: fs_sp(:,:,:,:) ! RHS of subsystems

    type(ElementTransferBuffer_3D), asynchronous, allocatable, save :: buf_r
    type(ElementTransferBuffer_3D), asynchronous, allocatable, save :: buf_us

    real(RNP), save :: rr_term
    logical  , save :: converged

    integer   :: wp = -1  ! working precision
    integer   :: no = -1  ! overlap
    integer   :: np = -1  ! element points per direction
    integer   :: ns = -1  ! subdomain points per direction

    integer   :: i, ne, ng, nl(3) = -1
    logical   :: reuse
    real(RNP) :: rr

    ! skip empty partition
    if (this % sem % mesh % part < 0) return

    associate(mesh => this % sem % mesh, schwarz => this % schwarz)

      ! initialization .........................................................

      ne = mesh % n_elem
      ng = mesh % n_ghost

      ! check for reusable workspace
      reuse = allocated(r) .and. wp == schwarz % wp .and. no == schwarz % no
      if (reuse) then
        reuse = size(r,1) == size(u,1) .and. size(r,4) == ne+ng
      end if

      !$omp master

      if (.not. reuse) then

        if (allocated( r      )) deallocate( r      )
        if (allocated( us_dp  )) deallocate( us_dp  )
        if (allocated( fs_dp  )) deallocate( fs_dp  )
        if (allocated( us_sp  )) deallocate( us_sp  )
        if (allocated( fs_sp  )) deallocate( fs_sp  )
        if (allocated( buf_r  )) deallocate( buf_r  )
        if (allocated( buf_us )) deallocate( buf_us )

        wp = schwarz % wp
        no = schwarz % no
        np = this % eop % po + 1

        nl = no
        ns = np + 2*no

        allocate(r(np, np, np, ne+ng), source = ZERO)

        if (wp == RSP) then
          allocate(us_sp(ns, ns, ns, ne+ng), source = 0E0)
          allocate(fs_sp(ns, ns, ns, ne)   , source = 0E0)
          buf_us = ElementTransferBuffer_3D(mesh, us_sp, nl)
        else
          allocate(us_dp(ns, ns, ns, ne+ng), source = 0D0)
          allocate(fs_dp(ns, ns, ns, ne)   , source = 0D0)
          buf_us = ElementTransferBuffer_3D(mesh, us_dp, nl)
        end if
        buf_r = ElementTransferBuffer_3D(mesh, r, nl)

      end if

      ! termination condition
      if (present(r_max)) then
        rr_term = max(ZERO, r_max)**2
      else
        rr_term = ZERO
      end if

      !$omp end master
      !$omp barrier

      ! Schwarz iterations .....................................................

      do i = 1, i_max

        call this % Apply(u, r(:,:,:,:ne), f, bv)

        ! termination check
        if (present(r_red)) then
          rr = ScalarProduct(r(:,:,:,:ne), r(:,:,:,:ne), mesh%comm_parts)
          !$omp master
          if (mesh%part == 0) then
            if (i == 1) then
              rr_term  = max(rr_term, max(ZERO, sqrt(rr) * r_red)**2)
            end if
            converged = rr <= rr_term
          end if
          call XMPI_Bcast(converged, root=0, comm=mesh%comm_parts)
          !$omp end master
          !$omp barrier
        end if
        if (converged) exit

        select case(wp)
        case(RSP)
          call schwarz % RestrictResidual(mesh, buf_r, r, fs_sp, sgn = -1)
          call TPO_Schwarz( schwarz % ops_sp % S      &
                          , schwarz % ops_sp % W      &
                          , schwarz % cfg             &
                          , schwarz % ops_sp % D_inv  &
                          , fs_sp                     &
                          , us_sp                     )
          call schwarz % MergeCorrections(mesh, buf_us, us_sp, u)
        case default
          call schwarz % RestrictResidual(mesh, buf_r, r, fs_dp, sgn = -1)
          call TPO_Schwarz( schwarz % ops_dp % S      &
                          , schwarz % ops_dp % W      &
                          , schwarz % cfg             &
                          , schwarz % ops_dp % D_inv  &
                          , fs_dp                     &
                          , us_dp                     )
          call schwarz % MergeCorrections(mesh, buf_us, us_dp, u)
        end select

      end do

      if (present(ni)) ni = min(i, i_max)

      ! clean-up ...............................................................

      ! keep workspace in case of standby
      if (present(standby)) then
        if (standby) return
      end if

      !$omp master
      if (allocated( r      )) deallocate( r      )
      if (allocated( us_dp  )) deallocate( us_dp  )
      if (allocated( fs_dp  )) deallocate( fs_dp  )
      if (allocated( us_sp  )) deallocate( us_sp  )
      if (allocated( fs_sp  )) deallocate( fs_sp  )
      if (allocated( buf_r  )) deallocate( buf_r  )
      if (allocated( buf_us )) deallocate( buf_us )
      !$omp end master

    end associate

  end subroutine Schwarz_Method

  !-----------------------------------------------------------------------------
  !> Schwarz-preconditioned conjugate gradient method

  subroutine SchwarzPCG_Method(this, u, f, bv, i_max, r_red, r_max, ni)

    class(DG_EllipticOperator_3D), intent(in) :: this
    class(SpectralElementBoundaryVariable_3D), intent(in) :: bv(:) !< BC
    real(RNP),           intent(inout) :: u(:,:,:,:) !< approximate solution
    real(RNP),           intent(in)    :: f(:,:,:,:) !< right hand side
    integer,             intent(in)    :: i_max      !< max num iterations
    real(RNP), optional, intent(in)    :: r_red      !< min residual reduction
    real(RNP), optional, intent(in)    :: r_max      !< max admissible residual
    integer,   optional, intent(out)   :: ni         !< executed num iterations

    contiguous :: u, f

    ! local variables ..........................................................

    real(RNP), dimension(:,:,:,:), allocatable, save :: r, p, q, s, z
    real(RNP), save :: rr_term
    logical  , save :: converged

    real(RNP), parameter :: eps = epsilon(ONE) * 1e-3
    real(RNP) :: alpha, beta, delta, rr !, rr_old
    logical   :: check_convergence, singular
    integer   :: i, i_max_

    ! skip empty partition
    if (this % sem % mesh % part < 0) return

    ! initialization ...........................................................

    associate(mesh => this % sem % mesh)

      ! work space
      !$omp master
      allocate(r, mold = u)
      allocate(p, mold = u)
      allocate(q, mold = u)
      allocate(s, mold = u)
      allocate(z, mold = u)
      !$omp end master
      !$omp barrier

      check_convergence = present(r_red) .or. present(r_max)
      singular = abs(this%lambda) < epsilon(ONE) .and. all(this%bc /= 'D')

      ! initial residual .......................................................

      ! r = Au - f
      call this % Apply(u, r, f, bv)
      if (singular) then
        call CalibrateArray(r, mesh%comm_parts)
      end if

      ! termination conditions
      if (check_convergence) then
        rr = ScalarProduct(r, r, mesh%comm_parts)
        !$omp master
        if (present(r_red)) then
          rr_term  = max(ZERO, sqrt(rr) * r_red)**2
        else
          rr_term = 0
        end if
        if (present(r_max)) then
          rr_term = max(rr_term, max(ZERO, r_max)**2)
        end if
        converged = rr <= rr_term
        call XMPI_Bcast(converged, root=0, comm=mesh%comm_parts)
        !$omp end master
        !$omp barrier
      else
        !$omp master
        converged = .false.
        !$omp end master
        !$omp barrier
      end if

      if (converged) then
        i_max_ = 0
        i      = 0
      else
        i_max_ = i_max
      end if

      ! iteration ...............................................................

      do i = 1, i_max_

        ! Schwarz preconditioner, z = (Aˢ)⁻¹ (Au - f) with homogeneous BC
        call SetArray(z, ZERO)
        call this % Schwarz_Method(z, r, i_max = 1, standby = i < i_max_)

        ! set/update search vector
        if (i == 1) then
          if (singular) then
            call CalibrateArray(z, mesh%comm_parts)
          end if
          call SetArray(p, z)                                 ! p = z
        else
          call SetArray(q, r)                                 ! q = r
          call MergeArrays(ONE, q, -ONE, s)                   ! q = r - s
          beta = ScalarProduct(q, z, mesh%comm_parts) / delta
          call MergeArrays(beta, p, ONE, z)                   ! p = beta p + z
        end if

        ! save old residual
        call SetArray(s, r)

        ! correction
        call this % Apply(p, q)
        delta = ScalarProduct(r, z, mesh%comm_parts)
        alpha = delta / ScalarProduct(p, q, mesh%comm_parts)
        call MergeArrays(ONE, u, -alpha, p)

        if (mod(i,50) == 0) then
          ! compute true residual to get rid of round-off errors
          call this % Apply(u, r, f, bv)
          if (singular) then
            call CalibrateArray(r, mesh%comm_parts)
          end if
        else
          call MergeArrays(ONE, r, -alpha, q)
        end if

        if (check_convergence) then
          rr = ScalarProduct(r, r, mesh%comm_parts)
          !$omp master
          converged = rr <= rr_term
          call XMPI_Bcast(converged, root=0, comm=mesh%comm_parts)
          !$omp end master
          !$omp barrier
        end if

        if (converged .or. i == i_max_) exit

      end do

      ! finalization ...........................................................

      if (present(ni)) ni = i

      !$omp master
      deallocate(p, q, r, s, z)
      !$omp end master

    end associate

  end subroutine SchwarzPCG_Method

!?!  !-----------------------------------------------------------------------------
!?!  !> Addition of boundary conditions to the right hande side
!?!  !>
!?!  !> `bv` is a boundary variable which contains the Dirichlet or Neumann
!?!  !> boundary values four each boundary. These values are applied to the
!?!  !> right hand side `f` according boundary type specified in `this % bc`.
!?!
!?!  subroutine AddBC(this, bv, f)
!?!    class(DG_EllipticOperator_3D), intent(in) :: this
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

end module DG__Elliptic_Operator__3D
