submodule(INS__Time_Integrator__3D) MP_PressureSolver
  use TPO__AAA__3D
  implicit none

contains

  subroutine PressureSolver(this, dt, bv_v, v, f_v, p_v)
    class(INS_TimeIntegrator_3D), intent(in) :: this
    !< time integration method
    real(RNP), intent(in) :: dt
    !< time step, possibly scaled by a factor
    class(SpectralElementBoundaryVariable_3D),intent(in) :: bv_v(:)
    !< velocity boundary conditions
    real(RNP), contiguous, intent(in) :: v(:,:,:,:,:)
    !< preliminary velocity
    real(RNP), contiguous, intent(in) :: f_v(:,:,:,:)
    !< source at velocity points
    real(RNP), contiguous, intent(inout) :: p_v(:,:,:,:)
    !< pressure at velocity points

    ! internal variables .......................................................

    real(RNP), allocatable, save :: p(:,:,:,:) ! pressure at p-points
    real(RNP), allocatable, save :: f(:,:,:,:) ! source at p-points
    class(SpectralElementBoundaryVariable_3D), allocatable, save :: bv_p(:)
    ! pressure boundary values at p-points

    character, allocatable :: bc_p(:) ! pressure boundary conditions
    logical :: mixed_order
    integer :: b

    associate( bc_v    => this % problem % bc_v       &
             , po_v    => this % ins_op % eop_v % po  &
             , po_p    => this % ins_op % eop_p % po  &
             , sem_p   => this % ins_op % sem_p       &
             , mesh    => this % ins_op % mesh        )

      ! initialization .........................................................

      mixed_order = po_p /= po_v

      !$omp master
      allocate(p(0:po_p, 0:po_p, 0:po_p, 1:mesh%n_elem))
      if (mixed_order) then
        allocate(f, mold = p)
      end if
      bv_p = SpectralElementBoundaryVariable_3D(sem_p, mesh%boundary, nc=1)
      !$omp end master
      !$omp barrier

      allocate(bc_p(mesh%n_bound))

      ! build pressure BC ......................................................

      call BuildPressureBC(this, dt, bv_v, v, bc_p, bv_p)

      ! solve ..................................................................

      if (mixed_order) then

        call TPO_AAA(this % ins_op % iop_vp % A, f_v, f)
        ! solve ...
        call TPO_AAA(this % ins_op % iop_pv % A, p, p_v)

      else
        ! solve ...
      end if

      ! finalization ...........................................................

      !$omp master
      deallocate(bv_p, p)
      if (mixed_order) then
        deallocate(f)
      end if
      !$omp end master

    end associate

  end subroutine PressureSolver_MO

  !-----------------------------------------------------------------------------
  !>

  subroutine BuildPressureBC(ins_ti, dt, bv_v, v, bc_p, bv_p)
    class(INS_TimeIntegrator_3D), intent(in) :: ins_ti
    !< time integration method
    real(RNP), intent(in) :: dt
    !< time step, possibly scaled by a factor
    class(SpectralElementBoundaryVariable_3D), intent(in) :: bv_v(:)
    !< velocity boundary values
    real(RNP), contiguous, intent(in) :: v(:,:,:,:,:)
    !< preliminary velocity
    character, intent(out) :: bc_p
    !< pressure boundary conditions
    class(SpectralElementBoundaryVariable_3D), intent(inout) :: bv_p(:)
    !< pressure boundary values

    integer :: b

    associate( problem => ins_ti % problem &
             , ins_op  => ins_ti % ins_op  )

      do b = 1, ins_op % mesh % n_bound

        select case(problem % bc_v(b))

        case('D') ! Dirichlet conditions for normal velocity ...................

          bc_p(b) = 'N'

          associate( boundary => ins_op % mesh % boundary(b)  &
                   , A        => ins_op % iop_vp % A          &
                   , n        => ins_op % sem_v % metrics % n &
                   , vb       => bv_v(b) % val(:,:,:,1:3)     &
                   , hp       => bv_p(b) % val(:,:,:,1)       )

            if (mesh % regular) then
              call PressureBC_NormalVelocity_R(boundary, dt, A, v, vb, hp)
            else
              call PressureBC_NormalVelocity_D(boundary, dt, A, n, v, vb, hp)
            end if
          end associate

        case('P') ! periodic conditions ........................................

          bc_p(b) = 'P'

        end select

      end do
    end associate

  end subroutine BuildPressureBC

  !-----------------------------------------------------------------------------
  !> Build pressure BC from normal velocity conditions -- regular mesh

  subroutine PressureBC_NormalVelocity_R(boundary, dt, A, v, vb, hp)
    class(MeshBoundary_3D), intent(in)  :: boundary
    real(RNP),              intent(in)  :: dt
    real(RNP),              intent(in)  :: A (0:,0:)
    real(RNP), contiguous,  intent(in)  :: v (0:,0:,0:,:,:)
    real(RNP), contiguous,  intent(in)  :: vb(0:,0:,:,:)
    real(RNP), contiguous,  intent(out) :: hp(0:,0:,:)

    real(RNP), allocatable :: hv(:,:), w(:,:)
    real(RNP) :: c
    integer :: po_p, po_v
    integer :: f, i, j, k, m

    po_v = ubound(vb, 1)
    po_p = ubound(hp, 1)

    allocate(hv(0:po_v,0:po_v), w(0:po_p,0:po_v))

    c = 1 / dt

    !$omp do
    do f = 1, boundary % n_face

      e = boundary % face(f) % element_id
      m = boundary % face(f) % element_face

      select case(m)

      case(1,2)
        i = (m - 1) * po_v
        n = (m - 1) * 2 - 1
        do k = 0, po_v
        do j = 0, po_v
          hv(j,k) = c * n * (v(i,j,k,e,1) - vb(j,k,f,1))
        end do
        end do

      case(3,4)
        j = (m - 3) * po_v
        n = (m - 3) * 2 - 1
        do k = 0, po_v
        do i = 0, po_v
          hv(i,k) = c * n * (v(i,j,k,e,2) - vb(i,k,f,2))
        end do
        end do

      case(5,6)
        k = (m - 5) * po_v
        n = (m - 5) * 2 - 1
        do j = 0, po_v
        do i = 0, po_v
          hv(i,j) = c * n * (v(i,j,k,e,3) - vb(i,j,f,3))
        end do
        end do

      end select

      if (po_v == po_p) then
        hp(:,:,f) = hv
      else
        call InterpolateToPressureSpace(po_v, po_p, A, hv, hp(:,:,f), w)
      end if

    end do

  end subroutine PressureBC_NormalVelocity_R

  !-----------------------------------------------------------------------------
  !> Build pressure BC from normal velocity conditions -- deformed mesh

  subroutine PressureBC_NormalVelocity_D(boundary, dt, A, n, v, vb, hp)
    class(MeshBoundary_3D), intent(in)  :: boundary
    real(RNP),              intent(in)  :: dt
    real(RNP),              intent(in)  :: A (0:,0:)
    real(RNP), contiguous,  intent(in)  :: n (0:,0:,:,:,:)
    real(RNP), contiguous,  intent(in)  :: v (0:,0:,0:,:,:)
    real(RNP), contiguous,  intent(in)  :: vb(0:,0:,:,:)
    real(RNP), contiguous,  intent(out) :: hp(0:,0:,:)

    real(RNP), allocatable :: hv(:,:), w(:,:)
    real(RNP) :: c
    integer :: po_p, po_v
    integer :: f, i, j, k, m

    po_v = ubound(vb, 1)
    po_p = ubound(hp, 1)

    allocate(hv(0:po_v,0:po_v), w(0:po_p,0:po_v))

    c = 1 / dt

    !$omp do
    do f = 1, boundary % n_face

      e = boundary % face(f) % element_id
      m = boundary % face(f) % element_face

      select case(m)

      case(1,2)
        i = (m - 1) * po_v
        do k = 0, po_v
        do j = 0, po_v
          hv(j,k) = c * ( n(j,k,m,e,1) * (v(i,j,k,e,1) - vb(j,k,f,1)) &
                        + n(j,k,m,e,2) * (v(i,j,k,e,2) - vb(j,k,f,2)) &
                        + n(j,k,m,e,3) * (v(i,j,k,e,3) - vb(j,k,f,3)) )
        end do
        end do

      case(3,4)
        j = (m - 3) * po_v
        do k = 0, po_v
        do i = 0, po_v
          hv(i,k) = c * ( n(i,k,m,e,1) * (v(i,j,k,e,1) - vb(i,k,f,1)) &
                        + n(i,k,m,e,2) * (v(i,j,k,e,2) - vb(i,k,f,2)) &
                        + n(i,k,m,e,3) * (v(i,j,k,e,3) - vb(i,k,f,3)) )
        end do
        end do

      case(5,6)
        k = (m - 5) * po_v
        do j = 0, po_v
        do i = 0, po_v
          hv(i,j) = c * ( n(i,j,m,e,1) * (v(i,j,k,e,1) - vb(i,j,f,1)) &
                        + n(i,j,m,e,2) * (v(i,j,k,e,2) - vb(i,j,f,2)) &
                        + n(i,j,m,e,3) * (v(i,j,k,e,3) - vb(i,j,f,3)) )
        end do
        end do

      end select

      if (po_v == po_p) then
        hp(:,:,f) = hv
      else
        call InterpolateToPressureSpace(po_v, po_p, A, hv, hp(:,:,f), w)
      end if

    end do

  end subroutine PressureBC_NormalVelocity_D

  !-----------------------------------------------------------------------------
  !> Interpolate face data from velocity into pressure space

  pure subroutine InterpolateToPressureSpace(po_v, po_p, A, uv, up, w)
    integer,   intent(in)    :: po_v              !< polynomial order of velocity
    integer,   intent(in)    :: po_p              !< polynomial order of pressure
    real(RNP), intent(in)    :: A (0:po_p,0:po_v) !< 1D interpolation operator
    real(RNP), intent(in)    :: uv(0:po_v,0:po_v) !< variable in velocity space
    real(RNP), intent(out)   :: up(0:po_p,0:po_p) !< variable in pressure space
    real(RNP), intent(inout) :: w (0:po_p,0:po_v) !< workspace

    real(RNP) :: tmp
    integer   :: i, j, k

    ! direction 1
    do j = 0, po_v
    do i = 0, po_p
      tmp = 0
      do k = 0, po_v
        tmp = tmp + A(i,k) * uv(k,j)
      end do
      w(i,j) = tmp
    end do
    end do

    ! direction 2
    do j = 0, po_p
    do i = 0, po_p
      tmp = 0
      do k = 0, po_v
        tmp = tmp + A(j,k) * w(i,k)
      end do
      up(i,j) = tmp
    end do
    end do

  end subroutine InterpolateToPressureSpace

  !=============================================================================

end submodule MP_PressureSolver
