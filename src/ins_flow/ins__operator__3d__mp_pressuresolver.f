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

  module subroutine PressureSolver( this, tau, bv_u, v, f, p &
                                  , i_max, r_red, r_max, ni  )

    ! arguments ................................................................

    class(INS_Operator_3D), intent(in) :: this !< INS operator

    real(RNP), intent(in) :: tau !<  effective time step width

    class(BoundaryVariable_3D), intent(inout) :: bv_u(:)
    !< boundary values
    !!   - components 1:3
    !!       * Γᴰ :  vᵇ  →  vᵇ         (unchanged)
    !!       * Γᴼ :  ×                 (unused)
    !!   - component 4
    !!       * Γᴰ :  ×   →  ∂p/∂n
    !!       * Γᴼ :  pᵇ  →  pᵇ         (unchanged)

    real(RNP), contiguous, intent(in)    :: v(:,:,:,:,:) !< preliminary velocity
    real(RNP), contiguous, intent(in)    :: f(:,:,:,:)   !< source at v-points
    real(RNP), contiguous, intent(inout) :: p(:,:,:,:)   !< pressure at v-points

    integer,               intent(in)    :: i_max  !< max num iterations
    real(RNP),   optional, intent(in)    :: r_red  !< min residual reduction
    real(RNP),   optional, intent(in)    :: r_max  !< max admissible residual
    integer,     optional, intent(out)   :: ni     !< executed num iterations

    ! internal variables .......................................................

    real(RNP), allocatable, save :: mm(:,:,:,:) ! pressure mass matrix
    real(RNP), allocatable, save :: g (:,:,:,:) ! source at p-points
    real(RNP), allocatable, save :: q (:,:,:,:) ! pressure at p-points

    type(BoundaryVariable_3D), allocatable, save :: bv_p(:), bv_q(:)
    ! pressure boundary conditions at v- and p-points

    logical   :: mixed_order
    real(RNP) :: ct
    integer   :: b, e

    associate( po          => this % eop_u % po  &
             , pq          => this % eop_p % po  &
             , mesh        => this % mesh        &
             , bc_p        => this % bc_p        &
             , sem_p       => this % sem_p       &
             , pressure_op => this % pressure_op )

      ! initialization .........................................................

      mixed_order = pq /= po
      ct = 1 / tau

      !$omp master
      allocate(mm(0:pq, 0:pq, 0:pq, 1:mesh%n_elem))
      allocate(g, mold = mm)
      allocate(bv_p(mesh % n_bound))
      do b = 1, mesh % n_bound
        call bv_u(b) % GetSlice(bv_p(b), first = 4, last = 4)
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

      ! build BC ...............................................................

      call BuildPressureBC(this, ct, v, bv_u, bv_p, bv_q)

      ! solve ..................................................................

      if (mixed_order) then
        ! transfer current approximation and source to order pq
        call TPO_AAA(this % iop_up % A, p, q) ! interpolation of pressure
        call TPO_AAA(this % pop_up % A, f, g) ! L² projection of RHS
        !$omp do
        do e = 1, mesh % n_elem
          g(:,:,:,e) = -ct * mm(:,:,:,e) * g(:,:,:,e)
        end do
        ! apply Schwarz-PCG with λ=0 and ν=1
        call pressure_op % SchwarzPCG_Method( ZERO, ONE, q, g, bv_q   &
                                            , i_max, r_red, r_max, ni )
        ! interpolate result to order po
        call TPO_AAA(this % iop_pu % A, q, p)
      else
        !$omp do
        do e = 1, mesh % n_elem
          g(:,:,:,e) = -ct * mm(:,:,:,e) * f(:,:,:,e)
        end do
        ! apply Schwarz-PCG with λ=0 and ν=1
        call pressure_op % SchwarzPCG_Method( ZERO, ONE, p, g, bv_p   &
                                            , i_max, r_red, r_max, ni )
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
  !> Build pressure boundary conditions

  subroutine BuildPressureBC(ins_op, ct, v, bv_u, bv_p, bv_q)
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

      select case(ins_op % bc_p(b))
      case('N')
        call BuildNeumannBC( boundary = ins_op % mesh % boundary(b)  &
                           , n        = ins_op % sem_u % metrics % n &
                           , ct       = ct                           &
                           , v        = v                            &
                           , vb       = bv_u(b) % val(:,:,:,1:3)     &
                           , dn_p     = bv_p(b) % val(:,:,:,1)       )
      case('D')
        ! Dirichlet BC already copied to bv_p(b)
      end select

      if (present(bv_q)) then
        call InterpolateFaceData( A = ins_op % iop_up % A    &
                                , p = bv_p(b) % val(:,:,:,1) &
                                , q = bv_q(b) % val(:,:,:,1) )
      end if

    end do

  end subroutine BuildPressureBC

  !-----------------------------------------------------------------------------
  !> Build pressure BC from normal velocity conditions

  subroutine BuildNeumannBC(boundary, ct, n, v, vb, dn_p)
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

  end subroutine BuildNeumannBC

  !-----------------------------------------------------------------------------
  !> Interpolation of face data

  pure subroutine InterpolateFaceData(A, p, q)
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

  !=============================================================================

end submodule MP_PressureSolver
