submodule(ML__DG__Elliptic_Solver__3D) MP_FAS_MGCG_Solver
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> FAS-MG solver for problems with constant diffusivity

  module subroutine FAS_MGCG_Solver_C( this, bc, lambda, nu, bv, f, u &
                                     , i_max, ni, r_2 )

    class(ML_DG_EllipticSolver_3D), intent(in) :: this

    character, intent(in) :: bc(:)  !< boundary conditions
    real(RNP), intent(in) :: lambda !< Helmholtz parameter
    real(RNP), intent(in) :: nu     !< diffusivity

    class(ML_BoundaryVariable_3D), intent(in)    :: bv !< boundary values
    class(ML_MeshVariable_3D),     intent(inout) :: f  !< RHS
    class(ML_MeshVariable_3D),     intent(inout) :: u  !< approx/final solution

    integer,   optional, intent(in)  :: i_max !< overrides max num cycles
    integer,   optional, intent(out) :: ni    !< num executed cycles
    real(RNP), optional, intent(out) :: r_2   !< Euclidean residual norm

    call FAS_MGCG_Solver_X( this, bc, lambda, nu, null(), bv, f, u &
                          , i_max, ni, r_2 )

  end subroutine FAS_MGCG_Solver_C

  !-----------------------------------------------------------------------------
  !> FAS-MG solver for problems with variable diffusivity

  module subroutine FAS_MGCG_Solver_V( this, bc, lambda, nu, bv, f, u &
                                     , i_max, ni, r_2 )

    class(ML_DG_EllipticSolver_3D), intent(in) :: this

    character, intent(in) :: bc(:)  !< boundary conditions
    real(RNP), intent(in) :: lambda !< Helmholtz parameter

    class(ML_MeshVariable_3D),     intent(in)    :: nu !< diffusivity
    class(ML_BoundaryVariable_3D), intent(in)    :: bv !< boundary values
    class(ML_MeshVariable_3D),     intent(inout) :: f  !< RHS
    class(ML_MeshVariable_3D),     intent(inout) :: u  !< approx/final solution

    integer,   optional, intent(in)  :: i_max !< overrides max num cycles
    integer,   optional, intent(out) :: ni    !< num executed cycles
    real(RNP), optional, intent(out) :: r_2   !< Euclidean residual norm

    call FAS_MGCG_Solver_X( this, bc, lambda, null(), nu, bv, f, u &
                          , i_max, ni, r_2)

  end subroutine FAS_MGCG_Solver_V

  !-----------------------------------------------------------------------------
  !> Generic FAS-MG solver for problems with constant or variable diffusivity
  !>
  !> Either `nu_0` or `nu_v` must be given.

  module subroutine FAS_MGCG_Solver_X( this, bc, lambda, nu_0, nu_v, bv, f, u &
                                     , i_max, ni, r_2 )

    class(ML_DG_EllipticSolver_3D), intent(in) :: this

    character, intent(in) :: bc(:)
      !< boundary conditions
    real(RNP), intent(in) :: lambda
      !< Helmholtz parameter
    real(RNP), optional, intent(in) :: nu_0
      !< constant diffusivity
    class(ML_MeshVariable_3D), optional, intent(in) :: nu_v
      !< variable diffusivity
    class(ML_BoundaryVariable_3D), intent(in) :: bv
      !< boundary values
    class(ML_MeshVariable_3D), intent(inout) :: f
      !< RHS
    class(ML_MeshVariable_3D), intent(inout) :: u
      !< approx/final solution
    integer, optional, intent(in):: i_max
      !< overrides max num cycles
    integer, optional, intent(out) :: ni
      !< num executed cycles
    real(RNP), optional, intent(out) :: r_2
      !< Euclidean residual norm

    ! internal variables :::::::::::::::::::::::::::::::::::::::::::::::::::::::

    type(ML_MeshVariable_3D), allocatable, save :: p, q, r, s, v, w
    logical, save :: converged

    real(RNP) :: rr, r_max, r_new, r_old
    real(RNP) :: alpha, beta, delta
    logical   :: check_convergence, singular
    integer   :: i, i_max_, l, l_top

    associate( sem        => this % ml_op % sem                        &
             , comm_parts => this % ml_op % sem(:) % mesh % comm_parts &
             , comm_world => this % ml_op % sem(1) % mesh % comm_world )

      ! prerequisites ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      l_top = size(sem)

      if (present(i_max)) then
        i_max_ = i_max
      else
        i_max_ = this % i_max
      end if

      check_convergence = i_max_ > 1 .and. max(this%r_red, this%r_max) > 0
      singular = lambda == ZERO .and. all(bc /= 'D'

      allocate(p, q, r, s, v, w)
      call p % Init(this%ml_op, nc = 1, l_top = l_top)
      call q % Init(this%ml_op, nc = 1, l_top = l_top)
      call r % Init(this%ml_op, nc = 1, l_top = l_top)
      call s % Init(this%ml_op, nc = 1, l_top = l_top)
      call v % Init(this%ml_op, nc = 1, l_top = l_top)
      call w % Init(this%ml_op, nc = 1, l_top = l_top)
      !$omp end master

      ! start ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      select case(this % start_method)
      case(START_CASC)
        call this % FAS_MG_Start_X(bc, lambda, nu_0, nu_v, bv, f, u, v, w, 0)
      case(START_FMG)
        call this % FAS_MG_Start_X(bc, lambda, nu_0, nu_v, bv, f, u, v, w, 1)
      end select

      ! flexible CG with MG preconditioner :::::::::::::::::::::::::::::::::::::

      ! r = f - Au, zero in frozen elements
      call this % FAS_MG_Residual_X(bc, lambda, nu_0, nu_v, bv, f, u, r)

      ! remove mean value in singular case
      call ML_CalibrateArray_3D(r)

      ! termination conditions
      if (check_convergence) then
        rr = ML_ScalarProduct_3D(r, r)
        r_old  = sqrt(rr)
        r_max  = max(r_old * this%r_red, this%r_max)
        !$omp master
        converged = r_old < r_max
        call XMPI_Bcast(converged, root = 0, comm = comm_world)
        if (log_level_outer_iteration > 0 .and. sem(l_top)%mesh%part == 0) then
          print '(A,T8,A,I5,A,ES12.5)', &
                '#MGCG','>>>  i  =',0,',  |r| =', sqrt(rr)
        end if
        !$omp end master
      else
        !$omp master
        converged = .false.
        !$omp end master
      end if
      !$omp barrier

      if (converged) then
        i_max_ = 0
        i      = 0
      end if

      ! iteration ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      do i = 1, i_max_

        ! FAS-MG preconditioner ................................................

        ! z = 0
        call ML_SetArray_3D(z, ZERO)

        ! z = MG(r, 0) with homogeneous boundary conditions
        call this % FAS_MG_Cycle_X( bc, lambda, nu_0, nu_v, f = r, u = z &
                                  , r = v, v = w, ni = ni, r_2 = r_2     )

        ! set/update search vector .............................................

        if (i == 1) then
          call ML_SetArray_3D(p, z)                         ! p = z
        else
          call ML_SetArray_3D(q, r)                         ! q = r
          call ML_MergeArrays_3D(ONE, q, -ONE, s)           ! q = r - s
          beta = ML_ScalarProduct_3D(q, z) / delta          ! β = (q,z) / δ
          call ML_MergeArrays_3D(beta, p, ONE, z)           ! p = β p + z
        end if

        ! update frozen elements (should not be necessary)
        !### p ###

        ! save old residual
        call ML_SetArray_3D(s, r)

        ! correction ...........................................................

        ! apply homogeneous operator: q = -Ap
        call this % FAS_MG_Residual_X(bc, lambda, nu_0, nu_v, u = p, r = q)

        delta = ML_ScalarProduct_3D(r, z)                   ! δ = (r,z)
        alpha = -delta / ScalarProduct(p, q)                ! α = δ / (p,Ap)
        call ML_MergeArrays_3D(ONE, u, alpha, p)            ! u = u + α p
        call ML_MergeArrays_3D(ONE, r, alpha, q)            ! r = r - α Ap

      end do

      ! finalization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      ! update frozen elements
      ### u ###

      !$omp master
      deallocate(p, q, r, s, v, w)
      !$omp end master

    end associate

  end subroutine FAS_MGCG_Solver_X

  !=============================================================================

end submodule MP_FAS_MGCG_Solver
