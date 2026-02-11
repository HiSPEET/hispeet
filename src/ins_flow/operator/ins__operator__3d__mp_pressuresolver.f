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

    real(RNP), allocatable, save :: mm(:,:,:,:) ! pressure mass matrix
    real(RNP), allocatable, save :: g (:,:,:,:) ! source at p-points
    real(RNP), allocatable, save :: q (:,:,:,:) ! pressure at p-points

    type(BoundaryVariable_3D), allocatable, save :: bv_p(:), bv_q(:)
    ! pressure boundary conditions at v- and p-points

    logical   :: mixed_order
    integer   :: i_max
    real(RNP) :: ct, r_max, r_red
    integer   :: b, e, na, ne

    associate( po    => this % eop_u % po     &
             , pq    => this % eop_p % po     &
             , mesh  => this % mesh           &
             , bc_p  => this % problem % bc_p &
             , sem_p => this % sem_p          )

      ! initialization .........................................................

      mixed_order = pq /= po
      ct = 1 / tau

      na = mesh % n_elem_active
      ne = mesh % n_elem

      i_max = this % i_max_p
      r_red = this % r_red
      r_max = this % r_max
      if (present(precon)) then
        if (precon) i_max = this % k_pre_p
      end if

      !$omp master
      allocate(mm(0:pq, 0:pq, 0:pq, 1:ne))
      allocate(g, mold = mm)
      allocate(bv_p(mesh % n_bound))
      do b = 1, mesh % n_bound
        call bv_p(b) % Init(mesh%boundary(b), po, nc = 1)
      end do
      if (pq /= po) then
        allocate(q, mold = mm)
        allocate(bv_q(mesh % n_bound))
        do b = 1, mesh % n_bound
          call bv_q(b) % Init(mesh%boundary(b), pq, nc = 1)
        end do
      end if
      !$omp end master
      !$omp barrier

      call sem_p % Get_DG_DiagonalMassMatrix(mm)

      ! build boundary values ..................................................

      call BuildPressureBV(this, ct, v, bv_u, bv_p, bv_q)

      ! solve ..................................................................

      if (mixed_order) then

        ! transfer current approximation and source to order pq
        call TPO_AAA(this % iop_up % A, p, q) ! interpolation of pressure
        call TPO_AAA(this % pop_up % A, f, g) ! L² projection of RHS
        !$omp do
        do e = 1, ne
          g(:,:,:,e) = -ct * mm(:,:,:,e) * g(:,:,:,e)
        end do

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

        !$omp do
        do e = 1, ne
          g(:,:,:,e) = -ct * mm(:,:,:,e) * f(:,:,:,e)
        end do

        select case(this % pressure_solver)
        case('AS')
          call this % elliptic_p % Schwarz_Method &
                          (bc_p, ZERO, ONE, p, g, bv_p, i_max, r_red, r_max, ni)
        case('CG')
          call this % elliptic_p % CG_Method &
                          (bc_p, ZERO, ONE, p, g, bv_p, i_max, r_red, r_max, ni)
        case('SPCG')
          call this % elliptic_p % SchwarzPCG_Method &
                          (bc_p, ZERO, ONE, p, g, bv_p, i_max, r_red, r_max, ni)
        case('MG','MGCG')
          call ML_PressureSolver(this, p, g, bv_p, i_max, ni)
        end select

      end if

      ! finalization ...........................................................

      !$omp master
      deallocate(mm, g, bv_p)
      if (mixed_order) then
        deallocate(q, bv_q)
      end if
      !$omp end master

    end associate

  end subroutine PressureSolver

  !-----------------------------------------------------------------------------
  !> Build pressure boundary values

  subroutine BuildPressureBV(ins_op, ct, v, bv_u, bv_p, bv_q)
    class(INS_Operator_3D), intent(in) :: ins_op
    !< time integration method
    real(RNP), intent(in) :: ct
    !< temporal scaling factor, usually ~ 1/dt
    real(RNP), contiguous, intent(in) :: v(:,:,:,:,:)
    !< preliminary velocity
    class(BoundaryVariable_3D), intent(in) :: bv_u(:)
    !< boundary values of flow variables
    class(BoundaryVariable_3D), intent(inout) :: bv_p(:)
    !< pressure boundary conditions on v-points
    class(BoundaryVariable_3D), optional, intent(inout) :: bv_q(:)
    !< pressure boundary conditions on v-points

    integer :: b

    do b = 1, ins_op % mesh % n_bound

      select case(ins_op % problem % bc_p(b))
      case('N')
        call BuildNeumannBV( boundary = ins_op % mesh % boundary(b)  &
                           , n        = ins_op % sem_u % metrics % n &
                           , ct       = ct                           &
                           , v        = v                            &
                           , vb       = bv_u(b) % val(:,:,:,1:3)     &
                           , dn_p     = bv_p(b) % val(:,:,:,1)       )
      case('D')
        call SetArray(bv_p(b) % val(:,:,:,1), bv_u(b) % val(:,:,:,4))
      end select

      if (present(bv_q)) then
        call InterpolateFaceData( A = ins_op % iop_up % A    &
                                , p = bv_p(b) % val(:,:,:,1) &
                                , q = bv_q(b) % val(:,:,:,1) )
      end if

    end do

  end subroutine BuildPressureBV

  !-----------------------------------------------------------------------------
  !> Build pressure boundary values from normal velocity conditions

  subroutine BuildNeumannBV(boundary, ct, n, v, vb, dn_p)
    class(MeshBoundary_3D), intent(in)  :: boundary
    real(RNP),              intent(in)  :: ct
    real(RNP), contiguous,  intent(in)  :: n (0:,0:,:,:,:)
    real(RNP), contiguous,  intent(in)  :: v (0:,0:,0:,:,:)
    real(RNP), contiguous,  intent(in)  :: vb(0:,0:,:,:)
    real(RNP), contiguous,  intent(out) :: dn_p(0:,0:,:)

    integer :: e, f, i, j, k, m, po

    po = ubound(v,1)

    !$omp do
    do f = 1, boundary % n_face

      e = boundary % face(f) % element_id
      m = boundary % face(f) % element_face

      select case(m)

      case(1,2)
        i = (m - 1) * po
        do k = 0, po
        do j = 0, po
          dn_p(j,k,f) = ct * ( n(j,k,m,e,1) * (v(i,j,k,e,1) - vb(j,k,f,1)) &
                             + n(j,k,m,e,2) * (v(i,j,k,e,2) - vb(j,k,f,2)) &
                             + n(j,k,m,e,3) * (v(i,j,k,e,3) - vb(j,k,f,3)) )
        end do
        end do

      case(3,4)
        j = (m - 3) * po
        do k = 0, po
        do i = 0, po
          dn_p(i,k,f) = ct * ( n(i,k,m,e,1) * (v(i,j,k,e,1) - vb(i,k,f,1)) &
                             + n(i,k,m,e,2) * (v(i,j,k,e,2) - vb(i,k,f,2)) &
                             + n(i,k,m,e,3) * (v(i,j,k,e,3) - vb(i,k,f,3)) )
        end do
        end do

      case(5,6)
        k = (m - 5) * po
        do j = 0, po
        do i = 0, po
          dn_p(i,j,f) = ct * ( n(i,j,m,e,1) * (v(i,j,k,e,1) - vb(i,j,f,1)) &
                             + n(i,j,m,e,2) * (v(i,j,k,e,2) - vb(i,j,f,2)) &
                             + n(i,j,m,e,3) * (v(i,j,k,e,3) - vb(i,j,f,3)) )
        end do
        end do

      end select

    end do

  end subroutine BuildNeumannBV

  !-----------------------------------------------------------------------------
  !> Interpolation of face data

  subroutine InterpolateFaceData(A, p, q)
    real(RNP), intent(in)  :: A(0:,0:)   !< 1D interpolation operator
    real(RNP), intent(in)  :: p(0:,0:,:) !< variable in velocity space
    real(RNP), intent(out) :: q(0:,0:,:) !< variable in pressure space

    real(RNP), allocatable :: w(:,:)
    real(RNP) :: tmp
    integer   :: po, pq, nf
    integer   :: f, i, j, k

    po = ubound(p,1)
    pq = ubound(q,1)
    nf = ubound(q,3)

    allocate(w(0:pq,0:po))

    !$omp do
    do f = 1, nf

      ! direction 1
      do j = 0, po
      do i = 0, pq
        tmp = 0
        do k = 0, po
          tmp = tmp + A(i,k) * p(k,j,f)
        end do
        w(i,j) = tmp
      end do
      end do

      ! direction 2
      do j = 0, pq
      do i = 0, pq
        tmp = 0
        do k = 0, po
          tmp = tmp + A(j,k) * w(i,k)
        end do
        q(i,j,f) = tmp
      end do
      end do

    end do

  end subroutine InterpolateFaceData

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
