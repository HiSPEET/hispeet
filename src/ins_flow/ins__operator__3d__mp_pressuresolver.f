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
    !< pressure boundary values at p-points

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
    ! pressure boundary values at p-points

    logical   :: mixed_order
    real(RNP) :: cs
    integer   :: b, e

    associate( po          => this % eop_v % po  &
             , pq          => this % eop_p % po  &
             , mesh        => this % mesh        &
             , bc_p        => this % bc_p        &
             , sem_p       => this % sem_p       &
             , laplacian_p => this % laplacian_p )

      ! initialization .........................................................

      mixed_order = pq /= po
      cs = 1 / tau

      !$omp master
      allocate(mm(0:pq, 0:pq, 0:pq, 1:mesh%n_elem))
      allocate(g, mold = mm)
      if (mixed_order) then
        allocate(q, mold = mm)
        allocate(bv_q(mesh % n_bound))
        do b = 1, mesh % n_bound
          bv_q(b) = BoundaryVariable_3D(mesh%boundary(b), pq, nc=1)
        end do
      end if
      !$omp end master
      !$omp barrier

      call sem_p % Get_DG_DiagonalMassMatrix(mm)

      ! build pressure BC ......................................................

      call BuildPressureBC(this, cs, v, bv_v, bv_p, bv_q)

      ! solve ..................................................................

      if (mixed_order) then
        ! interpolate current approximation and source to order pq
        call TPO_AAA(this % iop_vp % A, p, q)
        call TPO_AAA(this % iop_vp % A, f, g)
        ! apply mass matrix and time scale to source
        !$omp do
        do e = 1, mesh % n_elem
          g(:,:,:,e) = -cs * mm(:,:,:,e) * g(:,:,:,e)
        end do
        ! apply Schwarz-PCG with λ=0 and ν=1
        call laplacian_p % SchwarzPCG_Method( ZERO, ONE, q, g, bv_q   &
                                            , i_max, r_red, r_max, ni )
        ! interpolate result to order po
        call TPO_AAA(this % iop_pv % A, q, p)
      else
        ! apply mass matrix and time scale to source
        !$omp do
        do e = 1, mesh % n_elem
          g(:,:,:,e) = -cs * mm(:,:,:,e) * f(:,:,:,e)
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

  subroutine BuildPressureBC(ins_op, cs, v, bv_v, bv_p, bv_q)
    class(INS_Operator_3D), intent(in) :: ins_op
    !< time integration method
    real(RNP), intent(in) :: cs
    !< scaling factor, usually ~ 1/dt
    real(RNP), contiguous, intent(in) :: v(:,:,:,:,:)
    !< preliminary velocity
    class(BoundaryVariable_3D), intent(in) :: bv_v(:)
    !< velocity boundary values
    class(BoundaryVariable_3D), intent(inout) :: bv_p(:)
    !< pressure boundary values on v-points
    class(BoundaryVariable_3D), optional, intent(inout) :: bv_q(:)
    !< pressure boundary values on p-points

    real(RNP), contiguous, pointer :: vb(:,:,:,:), hp(:,:,:), hq(:,:,:)
    integer :: b

    associate( boundary => ins_op % mesh % boundary     &
             , A        => ins_op % iop_vp % A          &
             , n        => ins_op % sem_v % metrics % n )

      nullify(hq)

      do b = 1, ins_op % mesh % n_bound

        select case(ins_op % bc_p(b))

        case('N')

          vb => bv_v(b) % val(:,:,:,1:3)
          hp => bv_p(b) % val(:,:,:,1)

          if (present(bv_q)) then
            hq => bv_q(b) % val(:,:,:,1)
          end if

          if (ins_op % mesh % regular) then
            call PressureBC_NormalVelocity_R(boundary(b), cs, A, v, vb, hp, hq)
          else
            call PressureBC_NormalVelocity_D(boundary(b), cs, A, n, v, vb, hp, hq)
          end if

      end select

      end do
    end associate

  end subroutine BuildPressureBC

  !-----------------------------------------------------------------------------
  !> Build pressure BC from normal velocity conditions -- regular mesh

  subroutine PressureBC_NormalVelocity_R(boundary, cs, A, v, vb, hp, hq)
    class(MeshBoundary_3D),          intent(in)  :: boundary
    real(RNP),                       intent(in)  :: cs
    real(RNP),                       intent(in)  :: A (0:,0:)
    real(RNP), contiguous,           intent(in)  :: v (0:,0:,0:,:,:)
    real(RNP), contiguous,           intent(in)  :: vb(0:,0:,:,:)
    real(RNP), contiguous,           intent(out) :: hp(0:,0:,:)
    real(RNP), contiguous, optional, intent(out) :: hq(0:,0:,:)

    real(RNP), allocatable :: w(:,:)
    integer :: po, pq
    integer :: e, f, i, j, k, m, n

    po = ubound(hp, 1)

    if (present(hq)) then
      pq = ubound(hq, 1)
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
          hp(j,k,f) = cs * n * (v(i,j,k,e,1) - vb(j,k,f,1))
        end do
        end do

      case(3,4)
        j = (m - 3) * po
        n = (m - 3) * 2 - 1
        do k = 0, po
        do i = 0, po
          hp(i,k,f) = cs * n * (v(i,j,k,e,2) - vb(i,k,f,2))
        end do
        end do

      case(5,6)
        k = (m - 5) * po
        n = (m - 5) * 2 - 1
        do j = 0, po
        do i = 0, po
          hp(i,j,f) = cs * n * (v(i,j,k,e,3) - vb(i,j,f,3))
        end do
        end do

      end select

      if (pq /= po) then
        call InterpolateFaceData(po, pq, A, hp(:,:,f), hq(:,:,f), w)
      end if

    end do

  end subroutine PressureBC_NormalVelocity_R

  !-----------------------------------------------------------------------------
  !> Build pressure BC from normal velocity conditions -- deformed mesh

  subroutine PressureBC_NormalVelocity_D(boundary, cs, A, n, v, vb, hp, hq)
    class(MeshBoundary_3D),          intent(in)  :: boundary
    real(RNP),                       intent(in)  :: cs
    real(RNP),                       intent(in)  :: A (0:,0:)
    real(RNP), contiguous,           intent(in)  :: n (0:,0:,:,:,:)
    real(RNP), contiguous,           intent(in)  :: v (0:,0:,0:,:,:)
    real(RNP), contiguous,           intent(in)  :: vb(0:,0:,:,:)
    real(RNP), contiguous,           intent(out) :: hp(0:,0:,:)
    real(RNP), contiguous, optional, intent(out) :: hq(0:,0:,:)

    real(RNP), allocatable :: w(:,:)
    integer :: po, pq
    integer :: e, f, i, j, k, m

    po = ubound(vb, 1)

    if (present(hq)) then
      pq = ubound(hq, 1)
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
          hp(j,k,f) = cs * ( n(j,k,m,e,1) * (v(i,j,k,e,1) - vb(j,k,f,1)) &
                           + n(j,k,m,e,2) * (v(i,j,k,e,2) - vb(j,k,f,2)) &
                           + n(j,k,m,e,3) * (v(i,j,k,e,3) - vb(j,k,f,3)) )
        end do
        end do

      case(3,4)
        j = (m - 3) * po
        do k = 0, po
        do i = 0, po
          hp(i,k,f) = cs * ( n(i,k,m,e,1) * (v(i,j,k,e,1) - vb(i,k,f,1)) &
                           + n(i,k,m,e,2) * (v(i,j,k,e,2) - vb(i,k,f,2)) &
                           + n(i,k,m,e,3) * (v(i,j,k,e,3) - vb(i,k,f,3)) )
        end do
        end do

      case(5,6)
        k = (m - 5) * po
        do j = 0, po
        do i = 0, po
          hp(i,j,f) = cs * ( n(i,j,m,e,1) * (v(i,j,k,e,1) - vb(i,j,f,1)) &
                           + n(i,j,m,e,2) * (v(i,j,k,e,2) - vb(i,j,f,2)) &
                           + n(i,j,m,e,3) * (v(i,j,k,e,3) - vb(i,j,f,3)) )
        end do
        end do

      end select

      if (pq /= po) then
        call InterpolateFaceData(po, pq, A, hp(:,:,f), hq(:,:,f), w)
      end if

    end do

  end subroutine PressureBC_NormalVelocity_D

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
