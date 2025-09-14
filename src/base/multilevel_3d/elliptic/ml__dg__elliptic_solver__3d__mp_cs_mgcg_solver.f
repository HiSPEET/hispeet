submodule(ML__DG__Elliptic_Solver__3D) MP_CS_MGCG_Solver
!### CHECK
use, intrinsic :: ieee_arithmetic
!### CHECK END
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> CS-MGCG solver for problems with global refinement and constant diffusivity
  !>
  !> Use `l_top` to specify a top level lower than `size(this%ml_op%sem)`

  module subroutine CS_MGCG_Solver_C( this, bc, lambda, nu, bv, f, u &
                                    , i_max, l_top, ni, r_2          )
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    character, intent(in) :: bc(:)                  !< boundary conditions
    real(RNP), intent(in) :: lambda                 !< Helmholtz parameter
    real(RNP), intent(in) :: nu                     !< diffusivity
    class(ML_BoundaryVariable_3D), intent(in) :: bv !< boundary values
    class(ML_MeshVariable_3D), intent(inout) :: f   !< RHS
    class(ML_MeshVariable_3D), intent(inout) :: u   !< approx/final solution
    integer,   optional, intent(in)  :: i_max       !< overrides max num cycles
    integer,   optional, intent(in)  :: l_top       !< top level
    integer,   optional, intent(out) :: ni          !< num executed cycles
    real(RNP), optional, intent(out) :: r_2         !< Euclidean residual norm

    call CS_MGCG_Solver_X &
             (this, bc, lambda, nu, null(), bv, f, u, i_max, l_top, ni, r_2)

  end subroutine CS_MGCG_Solver_C

  !-----------------------------------------------------------------------------
  !> CS-MG solver for problems with global refinement and variable diffusivity
  !>
  !> Use `l_top` to specify a top level lower than `size(this%ml_op%sem)`

  module subroutine CS_MGCG_Solver_V( this, bc, lambda, nu, bv, f, u &
                                    , i_max, l_top, ni, r_2          )
    class(ML_DG_EllipticSolver_3D), intent(in) :: this
    character, intent(in) :: bc(:)                  !< boundary conditions
    real(RNP), intent(in) :: lambda                 !< Helmholtz parameter
    class(ML_MeshVariable_3D), intent(in) :: nu     !< diffusivity
    class(ML_BoundaryVariable_3D), intent(in) :: bv !< boundary values
    class(ML_MeshVariable_3D), intent(inout) :: f   !< RHS
    class(ML_MeshVariable_3D), intent(inout) :: u   !< approx/final solution
    integer,   optional, intent(in)  :: i_max       !< overrides max num cycles
    integer,   optional, intent(in)  :: l_top       !< top level
    integer,   optional, intent(out) :: ni          !< num executed cycles
    real(RNP), optional, intent(out) :: r_2         !< Euclidean residual norm

    call CS_MGCG_Solver_X &
             (this, bc, lambda, null(), nu, bv, f, u, i_max, l_top, ni, r_2)

  end subroutine CS_MGCG_Solver_V

  !-----------------------------------------------------------------------------
  !> Generic CS-MGCG solver for problems with constant or variable diffusivity
  !>
  !> Either `nu_0` or `nu_v` must be given.

  subroutine CS_MGCG_Solver_X( this, bc, lambda, nu_0, nu_v, bv, f, u &
                             , i_max, l_top, ni, r_2 )

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
    integer, optional, intent(in) :: i_max
      !< overrides preset maximum number of cycles
    integer, optional, intent(in) :: l_top
      !< top level different from size(this%ml_op%sem)
    integer, optional, intent(out) :: ni
      !< number of executed cycles
    real(RNP), optional, intent(out) :: r_2
      !< Euclidean residual norm

    ! internal variables :::::::::::::::::::::::::::::::::::::::::::::::::::::::

    type(ML_MeshVariable_3D), allocatable, save :: r, z
    real(RNP), dimension(:,:,:,:), allocatable, save :: p, q, s
    logical, save :: converged, singular, singular_loc

    real(RNP) :: rr, r_max, r_new, r_old
    real(RNP) :: alpha, beta, delta
    logical   :: check_convergence
    integer   :: i, i_max_, l_top_

    ! prerequisites ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if (present(l_top)) then
      l_top_ = min(l_top, size(this%ml_op%sem))
    else
      l_top_ = size(this%ml_op%sem)
    end if

    if (present(i_max)) then
      i_max_ = i_max
    else
      i_max_ = this % i_max
    end if

    check_convergence = i_max_ > 1 .and. max(this%r_red, this%r_max) > 0

    !$omp master
    allocate(r, z)
    call r % Init(this%ml_op, nc = 1, l_top = l_top_)
    call z % Init(this%ml_op, nc = 1, l_top = l_top_)
    !$omp end master

    associate( mesh_top => this % ml_op % sem(l_top_) % mesh              &
             , comm_top => this % ml_op % sem(l_top_) % mesh % comm_parts &
             , ell_top  => this % elliptic_op(l_top_)                     &
             , f_top    => f    % level(l_top_) % val(:,:,:,:,1)          &
             , u_top    => u    % level(l_top_) % val(:,:,:,:,1)          &
             , r_top    => r    % level(l_top_) % val(:,:,:,:,1)          &
             , z_top    => z    % level(l_top_) % val(:,:,:,:,1)          &
             , bv_top   => bv   % level(l_top_) % var                     )

      ! initialization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      !$omp master

      allocate(p, mold = u_top)
      allocate(q, mold = u_top)
      allocate(s, mold = u_top)

      singular_loc = lambda == ZERO  &
                     .and. all(bc /= 'D') &
                     .and. mesh_top % n_elem_frozen == 0

      call XMPI_Allreduce(singular_loc, singular, MPI_LAND, comm_top)

      !$omp end master
      !$omp barrier

      ! flexible CG with MG preconditioner :::::::::::::::::::::::::::::::::::::

      ! initial residual .......................................................

      ! r = f - Au
      if (present(nu_0)) then
        call ell_top % Residual( bc, lambda, nu_0, f_top, bv_top, u_top, r_top )
      else
        call ell_top % Residual( bc, lambda                                 &
                               , nu = nu_v % level(l_top_) % val(:,:,:,:,1) &
                               , f  = f_top                                 &
                               , bv = bv_top                                &
                               , u  = u_top                                 &
                               , r  = r_top                                 )
      end if
      if (singular) then
        call CalibrateArray(r_top, comm_top)
      end if

      ! termination conditions
      if (check_convergence) then
        rr = ScalarProduct(r_top, r_top, comm_top)
        r_old  = sqrt(rr)
        r_max  = max(r_old * this%r_red, this%r_max)
        !$omp master
        converged = r_old < r_max
        call XMPI_Bcast(converged, root = 0, comm = comm_top)
        if (log_level_outer_iteration > 0 .and. mesh_top%part == 0) then
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

      ! iteration ...............................................................

      do i = 1, i_max_

        ! z = 0
        call ML_SetArray_3D(z, ZERO, l_top = l_top_)

        ! MG preconditioner: z = MG(r, 0)
        call this % CS_MG_Solver_X( bc, lambda, nu_0, nu_v, f = r, u = z &
                                  , i_max = 1, l_top = l_top_            )

        ! set/update search vector
        if (i == 1) then
          if (singular) then
            call CalibrateArray(z_top, comm_top)
          end if
          call SetArray(p, z_top)                           ! p = z
        else
          call SetArray(q, r_top)                           ! q = r
          call MergeArrays(ONE, q, -ONE, s)                 ! q = r - s
          beta = ScalarProduct(q, z_top, comm_top) / delta  ! β = (q,z) / δ
          call MergeArrays(beta, p, ONE, z_top)             ! p = beta p + z
        end if

        ! save old residual
        call SetArray(s, r_top)

        ! apply homogeneous operator: q = Ap
        if (present(nu_0)) then
          call ell_top % Apply(bc, lambda, nu_0, u=p, r=q)
        else
          call ell_top % Apply( bc, lambda                                 &
                              , nu = nu_v % level(l_top_) % val(:,:,:,:,1) &
                              , u  = p                                     &
                              , r  = q                                     )
        end if

        ! correction
        delta = ScalarProduct(r_top, z_top, comm_top)       ! δ = (r,z)
        alpha = delta / ScalarProduct(p, q, comm_top)       ! α = δ / (p,Ap)
        call MergeArrays(ONE, u_top,  alpha, p)             ! u = u + α p
        call MergeArrays(ONE, r_top, -alpha, q)             ! r = r - α Ap
!### CHECK
print '(99(G0,X))', 'delta =',delta
print '(99(G0,X))', 'alpha =',alpha
print '(99(G0,X))', 'max|p| =',maxval(abs(p))
print '(99(G0,X))', 'max|q| =',maxval(abs(q))
print '(99(G0,X))', 'any(ieee_is_nan(p)) =',any(ieee_is_nan(p))
print '(99(G0,X))', 'any(ieee_is_nan(q)) =',any(ieee_is_nan(q))
print '(99(G0,X))', 'any(ieee_is_nan(u_top)) =',any(ieee_is_nan(u_top))
print '(99(G0,X))', 'any(ieee_is_nan(r_top)) =',any(ieee_is_nan(r_top))
!### CHECK END

        if (check_convergence) then

          r_new = sqrt( ScalarProduct(r_top, r_top, comm_top) )

          !$omp master
          converged = r_new <= r_max
          call XMPI_Bcast(converged, root = 0, comm = comm_top)
          if (log_level_outer_iteration > 0 .and. mesh_top%part == 0) then
            print '(A,T8,A,I5,A,ES12.5)', &
                  '#MGCG','>>>  i  =',i,',  |r| =', r_new
          end if
          !$omp end master
          !$omp barrier

          r_old = r_new
        end if

        if (converged .or. i == i_max_) exit

      end do

      ! finalization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      if (present(ni)) then
        !$omp master
        ni = i
        !$omp end master
      end if

      if (present(r_2)) then
        if (.not. check_convergence) then
          r_new = sqrt( ScalarProduct(r_top, r_top, comm_top) )
        end if
        !$omp master
        r_2 = r_new
        !$omp end master
      end if

    end associate

    !$omp barrier
    !$omp master
    deallocate(p, q, r, s, z)
    !$omp end master

  end subroutine CS_MGCG_Solver_X

  !=============================================================================

end submodule MP_CS_MGCG_Solver
