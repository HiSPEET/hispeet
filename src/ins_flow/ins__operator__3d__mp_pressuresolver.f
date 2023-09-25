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

  module subroutine PressureSolver( this, tau, bv_v, v, f, bv_p, p &
                                  , i_max, r_red, r_max, ni        )

    ! arguments ................................................................

    class(INS_Operator_3D), intent(in) :: this !< INS operator

    real(RNP), intent(in) :: tau !<  effective time step width

    class(BoundaryVariable_3D), intent(in) :: bv_v(:)
    !< velocity boundary values

    real(RNP), contiguous, intent(in)    :: v(:,:,:,:,:) !< preliminary velocity
    real(RNP), contiguous, intent(in)    :: f(:,:,:,:)   !< source at v-points

    class(BoundaryVariable_3D), intent(inout) :: bv_p(:)
    !< pressure boundary conditions at v-points

    real(RNP), contiguous, intent(inout) :: p(:,:,:,:) !< pressure at v-points

    integer,               intent(in)    :: i_max  !< max num iterations
    real(RNP),   optional, intent(in)    :: r_red  !< min residual reduction
    real(RNP),   optional, intent(in)    :: r_max  !< max admissible residual
    integer,     optional, intent(out)   :: ni     !< executed num iterations

    ! internal variables .......................................................

    real(RNP), allocatable, save :: mm(:,:,:,:) ! pressure mass matrix
    real(RNP), allocatable, save :: g (:,:,:,:) ! source at p-points
    real(RNP), allocatable, save :: q (:,:,:,:) ! pressure at p-points

    type(BoundaryVariable_3D), allocatable, save :: bv_q(:)
    ! pressure boundary conditions at p-points

    logical   :: mixed_order
    real(RNP) :: ct
    integer   :: b, e

    associate( po          => this % eop_v % po  &
             , pq          => this % eop_p % po  &
             , mesh        => this % mesh        &
             , bc_p        => this % bc_p        &
             , sem_p       => this % sem_p       &
             , laplacian_p => this % laplacian_p )

      ! initialization .........................................................

      mixed_order = pq /= po
      ct = 1 / tau

      !$omp master
      allocate(mm(0:pq, 0:pq, 0:pq, 1:mesh%n_elem))
      allocate(g, mold = mm)
      if (mixed_order) then
        allocate(q, mold = mm)
        allocate(bv_q(mesh % n_bound))
        do b = 1, mesh % n_bound
          call bv_q(b) % Create(mesh%boundary(b), pq, nc=1)
        end do
      end if
      !$omp end master
      !$omp barrier

      call sem_p % Get_DG_DiagonalMassMatrix(mm)

      ! build pressure BC ......................................................

      call BuildPressureBC(this, ct, v, bv_v, bv_p, bv_q)

      ! solve ..................................................................

      if (mixed_order) then
        ! interpolate current approximation and source to order pq
        call TPO_AAA(this % iop_vp % A, p, q)
        call TPO_AAA(this % iop_vp % A, f, g)
        ! apply mass matrix and time scale to source
        !$omp do
        do e = 1, mesh % n_elem
          g(:,:,:,e) = -ct * mm(:,:,:,e) * g(:,:,:,e)
        end do
        ! apply Schwarz-PCG with λ=0 and ν=1
!### CHECK
!!         call laplacian_p % CG_Method( ZERO, ONE, q, g, bv_q   &
!!                                     , i_max, r_red, r_max, ni )
!!         call laplacian_p % Schwarz_Method( ZERO, ONE, q, g, bv_q   &
!!                                          , i_max, r_red, r_max, ni )
        call laplacian_p % SchwarzPCG_Method( ZERO, ONE, q, g, bv_q   &
                                            , i_max, r_red, r_max, ni )
!### CHECK END
        ! interpolate result to order po
        call TPO_AAA(this % iop_pv % A, q, p)
      else
        ! apply mass matrix and time scale to source
        !$omp do
        do e = 1, mesh % n_elem
          g(:,:,:,e) = -ct * mm(:,:,:,e) * f(:,:,:,e)
        end do
        ! apply Schwarz-PCG with λ=0 and ν=1
        call laplacian_p % SchwarzPCG_Method( ZERO, ONE, p, g, bv_p   &
                                            , i_max, r_red, r_max, ni )
      end if

      ! finalization ...........................................................

      !$omp master
      deallocate(mm, g)
      if (mixed_order) then
        deallocate(q, bv_q)
      end if
      !$omp end master

    end associate

  end subroutine PressureSolver

  !-----------------------------------------------------------------------------
  !> Build pressure boundary conditions

  subroutine BuildPressureBC(ins_op, ct, v, bv_v, bv_p, bv_q)
    class(INS_Operator_3D), intent(in) :: ins_op
    !< time integration method
    real(RNP), intent(in) :: ct
    !< temporal scaling factor, usually ~ 1/dt
    real(RNP), contiguous, intent(in) :: v(:,:,:,:,:)
    !< preliminary velocity
    class(BoundaryVariable_3D), intent(in) :: bv_v(:)
    !< velocity boundary values
    class(BoundaryVariable_3D), intent(inout) :: bv_p(:)
    !< pressure boundary conditions on v-points
    class(BoundaryVariable_3D), optional, intent(inout) :: bv_q(:)
    !< pressure boundary conditions on p-points

    real(RNP), contiguous, pointer :: vb(:,:,:,:), pb(:,:,:), qb(:,:,:)
    integer :: b

    associate( boundary => ins_op % mesh % boundary     &
             , A        => ins_op % iop_vp % A          &
             , n        => ins_op % sem_v % metrics % n &
             , delta    => ins_op % delta_outflow       )

      nullify(qb)

      do b = 1, ins_op % mesh % n_bound

        select case(ins_op % bc_p(b))

        case('N')

          vb => bv_v(b) % val(:,:,:,1:3)
          pb => bv_p(b) % val(:,:,:,1)

          if (present(bv_q)) then
            qb => bv_q(b) % val(:,:,:,1)
          else
            qb => null()
          end if

          if (ins_op % mesh % regular) then
            call BuildNeumannBC_R(boundary(b), ct, A, v, vb, pb, qb)
          else
            call BuildNeumannBC_D(boundary(b), ct, A, n, v, vb, pb, qb)
          end if

        case('D')

          call BuildDirichletBC(boundary(b), delta, A, n, vb, pb, qb)

      end select

      end do
    end associate

  end subroutine BuildPressureBC

  !-----------------------------------------------------------------------------
  !> Build pressure BC from normal velocity conditions -- regular mesh

  subroutine BuildNeumannBC_R(boundary, ct, A, v, vb, dn_p, dn_q)
    class(MeshBoundary_3D),          intent(in)  :: boundary
    real(RNP),                       intent(in)  :: ct
    real(RNP),                       intent(in)  :: A (0:,0:)
    real(RNP), contiguous,           intent(in)  :: v (0:,0:,0:,:,:)
    real(RNP), contiguous,           intent(in)  :: vb(0:,0:,:,:)
    real(RNP), contiguous,           intent(out) :: dn_p(0:,0:,:)
    real(RNP), contiguous, optional, intent(out) :: dn_q(0:,0:,:)

    real(RNP), allocatable :: w(:,:)
    integer :: po, pq
    integer :: e, f, i, j, k, m, n

    po = ubound(dn_p, 1)

    if (present(dn_q)) then
      pq = ubound(dn_q, 1)
      allocate(w(0:pq,0:po))
    else
      pq = po
    end if

    !$omp do
    do f = 1, boundary % n_face

      e = boundary % face(f) % element_id
      m = boundary % face(f) % element_face

      select case(m)

      case(1,2)
        i = (m - 1) * po
        n = (m - 1) * 2 - 1
        do k = 0, po
        do j = 0, po
          dn_p(j,k,f) = ct * n * (v(i,j,k,e,1) - vb(j,k,f,1))
        end do
        end do

      case(3,4)
        j = (m - 3) * po
        n = (m - 3) * 2 - 1
        do k = 0, po
        do i = 0, po
          dn_p(i,k,f) = ct * n * (v(i,j,k,e,2) - vb(i,k,f,2))
        end do
        end do

      case(5,6)
        k = (m - 5) * po
        n = (m - 5) * 2 - 1
        do j = 0, po
        do i = 0, po
          dn_p(i,j,f) = ct * n * (v(i,j,k,e,3) - vb(i,j,f,3))
        end do
        end do

      end select

      if (pq /= po) then
        call InterpolateFaceData(po, pq, A, dn_p(:,:,f), dn_q(:,:,f), w)
      end if

    end do

  end subroutine BuildNeumannBC_R

  !-----------------------------------------------------------------------------
  !> Build pressure BC from normal velocity conditions -- deformed mesh

  subroutine BuildNeumannBC_D(boundary, ct, A, n, v, vb, dn_p, dn_q)
    class(MeshBoundary_3D),          intent(in)  :: boundary
    real(RNP),                       intent(in)  :: ct
    real(RNP),                       intent(in)  :: A (0:,0:)
    real(RNP), contiguous,           intent(in)  :: n (0:,0:,:,:,:)
    real(RNP), contiguous,           intent(in)  :: v (0:,0:,0:,:,:)
    real(RNP), contiguous,           intent(in)  :: vb(0:,0:,:,:)
    real(RNP), contiguous,           intent(out) :: dn_p(0:,0:,:)
    real(RNP), contiguous, optional, intent(out) :: dn_q(0:,0:,:)

    real(RNP), allocatable :: w(:,:)
    integer :: po, pq
    integer :: e, f, i, j, k, m

    po = ubound(vb, 1)

    if (present(dn_q)) then
      pq = ubound(dn_q, 1)
      allocate(w(0:pq,0:po))
    else
      pq = po
    end if

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

      if (pq /= po) then
        call InterpolateFaceData(po, pq, A, dn_p(:,:,f), dn_q(:,:,f), w)
      end if

    end do

  end subroutine BuildNeumannBC_D

  !-----------------------------------------------------------------------------
  !> Build pressure BC at outflow boundaries

  subroutine BuildDirichletBC(boundary, delta, A, n, vb, pb, qb)
    class(MeshBoundary_3D),          intent(in)  :: boundary
    real(RNP),                       intent(in)  :: delta
    real(RNP),                       intent(in)  :: A (0:,0:)
    real(RNP), contiguous,           intent(in)  :: n (0:,0:,:,:,:)
    real(RNP), contiguous,           intent(in)  :: vb(0:,0:,:,:)
    real(RNP), contiguous,           intent(out) :: pb(0:,0:,:)
    real(RNP), contiguous, optional, intent(out) :: qb(0:,0:,:)

    real(RNP), allocatable :: theta(:,:), vn(:,:), vv(:,:), w(:,:)
    real(RNP) :: v_max = 0 ! shared with OpenMP !
    real(RNP) :: cv
    integer   :: po, pq
    integer   :: e, f, m

    po = ubound(vb, 1)

    if (present(qb)) then
      pq = ubound(qb, 1)
      allocate(w(0:pq,0:po))
    else
      pq = po
    end if

    allocate(theta(0:po,0:po))
    allocate(vn, mold = theta)
    allocate(vv, mold = theta)

    ! determine maximum velocity ...............................................

    !$omp master
    v_max = 0
    !$omp end master

    !$omp do reduction(max:v_max)
    do f = 1, boundary % n_face
      e = boundary % face(f) % element_id
      m = boundary % face(f) % element_face

      v_max = max(v_max, maxval(abs( n(:,:,m,e,1) * vb(:,:,f,1)  &
                                   + n(:,:,m,e,2) * vb(:,:,f,2)  &
                                   + n(:,:,m,e,3) * vb(:,:,f,3) ))
    end do

    ! complete pressure BC .....................................................

    !$omp do
    do f = 1, boundary % n_face

      e = boundary % face(f) % element_id
      m = boundary % face(f) % element_face

      vn = n(:,:,m,e,1) * vb(:,:,f,1)
         + n(:,:,m,e,2) * vb(:,:,f,2)
         + n(:,:,m,e,3) * vb(:,:,f,3)

      vv = vb(:,:,m,e,1) * vb(:,:,f,1)
         + vb(:,:,m,e,2) * vb(:,:,f,2)
         + vb(:,:,m,e,3) * vb(:,:,f,3)

      theta = HALF * tanh(cv * vn)

      pb = pb - HALF * (vv + vn) * theta

      if (pq /= po) then
        call InterpolateFaceData(po, pq, A, pb(:,:,f), qb(:,:,f), w)
      end if

    end do

  end subroutine BuildDirichletBC

  !-----------------------------------------------------------------------------
  !> Interpolation of face data

  pure subroutine InterpolateFaceData(po, pq, A, p, q, w)
    integer,   intent(in)    :: po           !< polynomial order of p
    integer,   intent(in)    :: pq           !< polynomial order of q
    real(RNP), intent(in)    :: A(0:pq,0:po) !< 1D interpolation operator
    real(RNP), intent(in)    :: p(0:po,0:po) !< variable in velocity space
    real(RNP), intent(out)   :: q(0:pq,0:pq) !< variable in pressure space
    real(RNP), intent(inout) :: w(0:pq,0:po) !< workspace

    real(RNP) :: tmp
    integer   :: i, j, k

    ! direction 1
    do j = 0, po
    do i = 0, pq
      tmp = 0
      do k = 0, po
        tmp = tmp + A(i,k) * p(k,j)
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
      q(i,j) = tmp
    end do
    end do

  end subroutine InterpolateFaceData

  !=============================================================================

end submodule MP_PressureSolver
