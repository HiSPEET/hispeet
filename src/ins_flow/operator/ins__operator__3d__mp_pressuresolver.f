!> summary:  Incompressible Navier-Stokes pressure solver
!> author:   Joerg Stiller
!> date:     2022/09/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(INS__Operator__3D) MP_PressureSolver
  use TPO__AAA__3D
  use Mesh_Boundary__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Projection-based pressure solver

  module subroutine PressureSolver(this, tau, bv_u, v, f, p, precon, ni)

    ! arguments ................................................................

    class(INS_Operator_3D), intent(in) :: this !< INS operator

    real(RNP), intent(in) :: tau !< effective time step width

    class(BoundaryVariable_3D), intent(in) :: bv_u(:)
    !< boundary values
    !!   - Γᴰ :  [ v₁, v₂, v₃, - ]    (1:3 unchanged)
    !!   - Γᴼ :  [ - , - , - , p ]    (all unchanged)

    real(RNP), contiguous, intent(in)    :: v(:,:,:,:,:) !< preliminary velocity
    real(RNP), contiguous, intent(in)    :: f(:,:,:,:)   !< source at v-points
    real(RNP), contiguous, intent(inout) :: p(:,:,:,:)   !< pressure at v-points
    logical,     optional, intent(in)    :: precon !< switch preconditioner mode
    integer,     optional, intent(out)   :: ni     !< num iterations executed

    ! internal variables .......................................................

    real(RNP), allocatable, save :: w(:,:,:,:) ! workspace at v-points
    real(RNP), allocatable, save :: g(:,:,:,:) ! source at p-points
    real(RNP), allocatable, save :: q(:,:,:,:) ! pressure at p-points

    type(BoundaryVariable_3D), allocatable, save :: bv_p(:), bv_q(:)
    ! pressure boundary conditions at v- and p-points

    logical   :: mixed_order
    integer   :: i_max
    real(RNP) :: r_max, r_red
    integer   :: b, e, na, ne

    associate( po    => this % eop_u % po     &
             , pq    => this % eop_p % po     &
             , mesh  => this % mesh           &
             , bc_p  => this % problem % bc_p &
             , sem_p => this % sem_p          )

      ! initialization .........................................................

      mixed_order = pq /= po

      na = mesh % n_elem_active
      ne = mesh % n_elem

      i_max = this % i_max_p
      r_red = this % r_red
      r_max = this % r_max
      if (present(precon)) then
        if (precon) i_max = this % k_pre_p
      end if

      !$omp master
      allocate(w, mold = f)
      allocate(g(0:pq, 0:pq, 0:pq, 1:ne))
      allocate(bv_p(mesh % n_bound))
      do b = 1, mesh % n_bound
        call bv_p(b) % Init(mesh%boundary(b), po, nc = 1)
      end do
      if (pq /= po) then
        allocate(q, mold = g)
        allocate(bv_q(mesh % n_bound))
        do b = 1, mesh % n_bound
          call bv_q(b) % Init(mesh%boundary(b), pq, nc = 1)
        end do
      end if
      !$omp end master
      !$omp barrier

      ! build boundary values ..................................................

      call this % Get_PressureBoundaryValues(tau, v, bv_u, bv_p, bv_q)

      ! scale and project sources in velocity space ............................

      call this % sem_u % Get_DG_DiagonalMassMatrix(w)

      !$omp do
      do e = 1, ne
        w(:,:,:,e) = -1/tau * w(:,:,:,e) * f(:,:,:,e)
      end do

      ! solve ..................................................................

      if (mixed_order) then

        ! interpolation of initial pressure
        call TPO_AAA(this % iop_up % A, p, q)

        ! transfer sources to pressure space: w ≈ -1/τ ∫ 𝜑ᵖ f dx
        call TPO_AAA(transpose(this % iop_pu % A), w, g)

        select case(this % pressure_solver)
        case('AS')
          call this % elliptic_p % Schwarz_Method &
                          (bc_p, ZERO, ONE, q, g, bv_q, i_max, r_red, r_max, ni)
        case('CG')
          call this % elliptic_p % CG_Method &
                          (bc_p, ZERO, ONE, q, g, bv_q, i_max, r_red, r_max, ni)
        case('SPCG')
          call this % elliptic_p % SchwarzPCG_Method &
                          (bc_p, ZERO, ONE, q, g, bv_q, i_max, r_red, r_max, ni)
        case('MG','MGCG')
          call ML_PressureSolver(this, q, g, bv_q, ni)
        end select

        ! interpolate result to order po
        call TPO_AAA(this % iop_pu % A, q(:,:,:,1:na), p(:,:,:,1:na))

      else

        ! g = w
        select case(this % pressure_solver)
        case('AS')
          call this % elliptic_p % Schwarz_Method &
                          (bc_p, ZERO, ONE, p, w, bv_p, i_max, r_red, r_max, ni)
        case('CG')
          call this % elliptic_p % CG_Method &
                          (bc_p, ZERO, ONE, p, w, bv_p, i_max, r_red, r_max, ni)
        case('SPCG')
          call this % elliptic_p % SchwarzPCG_Method &
                          (bc_p, ZERO, ONE, p, w, bv_p, i_max, r_red, r_max, ni)
        case('MG','MGCG')
          call ML_PressureSolver(this, p, w, bv_p, i_max, ni)
        end select

      end if

      ! finalization ...........................................................

      !$omp master
      deallocate(w, g, bv_p)
      if (mixed_order) then
        deallocate(q, bv_q)
      end if
      !$omp end master

    end associate

  end subroutine PressureSolver

  !-----------------------------------------------------------------------------
  !> Multilevel pressure solver

  subroutine ML_PressureSolver(this, p, f, bv, i_max, ni)
    class(INS_Operator_3D),     intent(in)    :: this       !< INS operator
    real(RNP), contiguous,      intent(inout) :: p(:,:,:,:) !< pressure
    real(RNP), contiguous,      intent(in)    :: f(:,:,:,:) !< sources
    class(BoundaryVariable_3D), intent(in)    :: bv(:)      !< boundary values
    integer,                    intent(in)    :: i_max      !< max num cycles
    integer,          optional, intent(out)   :: ni         !< num iterations

    type(ML_MeshVariable_3D),     allocatable, save :: ml_p, ml_f
    type(ML_BoundaryVariable_3D), allocatable, save :: ml_bv

    integer :: b, l

    associate( l_top       => this % level       &
             , ml_solver_p => this % ml_solver_p )

      ! workspace ..............................................................

      !$omp master
      allocate(ml_p, ml_f, ml_bv)
      call ml_p  % Init(ml_solver_p % ml_op, nc = 1, l_top = l_top)
      call ml_f  % Init(ml_solver_p % ml_op, nc = 1, l_top = l_top)
      call ml_bv % Init(ml_solver_p % ml_op, nc = 1, l_top = l_top)
      !$omp end master
      !$omp barrier

      ! initialization .........................................................

      do l = 1, l_top - 1
        call SetArray(ml_p % level(l) % val, ZERO)
        call SetArray(ml_f % level(l) % val, ZERO)
        ! ml_bv already initialized with zero
      end do

      call SetArray(ml_p % level(l_top) % val(:,:,:,:,1), p)
      call SetArray(ml_f % level(l_top) % val(:,:,:,:,1), f)

      do b = 1, size(bv)
        call SetArray(ml_bv % level(l_top) % var(b) % val, bv(b) % val)
      end do

      ! solution ...............................................................

      select case(this % pressure_solver)
      case('MG')
        call ml_solver_p % CS_MG_Solver( bc     = this % problem % bc_p  &
                                       , lambda = ZERO                   &
                                       , nu     = ONE                    &
                                       , bv     = ml_bv                  &
                                       , f      = ml_f                   &
                                       , u      = ml_p                   &
                                       , i_max  = i_max                  &
                                       , l_top  = l_top                  &
                                       , ni     = ni                     )
      case('MGCG')
        call ml_solver_p % CS_MGCG_Solver( bc     = this % problem % bc_p  &
                                         , lambda = ZERO                   &
                                         , nu     = ONE                    &
                                         , bv     = ml_bv                  &
                                         , f      = ml_f                   &
                                         , u      = ml_p                   &
                                         , i_max  = i_max                  &
                                         , l_top  = l_top                  &
                                         , ni     = ni                     )
      end select

      ! copy result ............................................................

      call SetArray(p, ml_p % level(l_top) % val(:,:,:,:,1))

      ! finalization ...........................................................

      !$omp master
      deallocate(ml_p, ml_f, ml_bv)
      !$omp end master

    end associate

  end subroutine ML_PressureSolver

  !=============================================================================

end submodule MP_PressureSolver
