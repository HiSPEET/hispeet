submodule(INS__Time_Integrator__3D) MP_PressureSolver

contains

  subroutine PressureSolver(this, dt, div_F_v, bv_v, v_0, p, w)
    class(INS_TimeIntegrator_3D), intent(in) :: this
    real(RNP), intent(in) :: dt
    real(RNP), contiguous, intent(in) :: v_0(:,:,:,:,:)
    real(RNP), contiguous, intent(in) :: div_F_v(:,:,:,:)
    class(SpectralElementBoundaryVariable_3D),intent(in) :: bv_v(:)
    real(RNP), contiguous, intent(inout) :: p(:,:,:,:)
    !< pressure (nv,nv,nv,n_elem)
    real(RNP), contiguous, intent(inout) :: w(:,:,:,:,:)
    !< workspace (nv,nv,nv,n_elem,3)

    class(SpectralElementBoundaryVariable_3D), allocatable, save :: bv_p(:)

    character, allocatable :: bc_p(:)
    integer :: b

    associate( bc_v    => this % problem % bc_v          &
             , sem_p   => this % ins_op % sem_p          &
             , n_bound => this % ins_op % mesh % n_bound &
             )

      ! initialization .........................................................

      !$omp master
      bv_p = SpectralElementBoundaryVariable_3D
      !$omp end master
      !$omp barrier

      ! build pressure BC ......................................................

      ! 1) precompute (v_0 - F_v)/dt on Dirichlet faces
      ! 2) extract normal component
      ! 3) apply 2D interpolation
      !
      ! x) interpolate RHS, remove constant
      ! y) solve
      ! z) interpolate result


      allocate(bc_p(n_bound))

      do b = 1, n_bound
        select case(bc_v)
        case('D')
          bc_p(b) = 'N'
        case('P')
          bc_p(b) = 'P'
        end select
      end do

    ! compute boundary conditions
    !   - difference to BC
    !   – extract normal component
    !   – interpolate
    !   -



      ! finalization ...........................................................

      !$omp master
      deallocate(bv_p)
      !$omp end master

    end associate

  end subroutine PressureSolver_MO

  !-----------------------------------------------------------------------------
  !>

  subroutine BuildPressureBC(ins_ti, dt, v, bv_v, bv_p, bc_p)
    class(INS_TimeIntegrator_3D), intent(in) :: ins_ti
    real(RNP), intent(in) :: dt
    real(RNP), contiguous, intent(in) :: v(:,:,:,:,:)
    class(SpectralElementBoundaryVariable_3D), intent(in) :: bv_v(:)
    class(SpectralElementBoundaryVariable_3D), intent(inout) :: bv_p(:)
    character, intent(out) :: bc_p

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
    real(RNP),              intent(in)  :: A(0:,0:)
    real(RNP), contiguous,  intent(in)  :: v (0:,0:,0:,:,:)
    real(RNP), contiguous,  intent(in)  :: vb(0:,0:,:,:)
    real(RNP), contiguous,  intent(out) :: hp(0:,0:,:)

    real(RNP), allocatable :: hv(:,:), w(:,:)
    real(RNP) :: c
    integer :: pp, pv
    integer :: f, i, j, k, m

    pv = ubound(vb, 1)
    pp = ubound(hp, 1)

    allocate(hv(0:pv,0:pv), w(0:pp,0:pv))

    c = 1 / dt

    !$omp do
    do f = 1, boundary % n_face

      e = boundary % face(f) % element_id
      m = boundary % face(f) % element_face

      select case(m)

      case(1,2)
        i = (m - 1) * pv
        n = (m - 1) * 2 - 1
        do k = 0, pv
        do j = 0, pv
          hv(j,k) = c * n * (v(i,j,k,e,1) - vb(j,k,f,1))
        end do
        end do

      case(3,4)
        j = (m - 3) * pv
        n = (m - 3) * 2 - 1
        do k = 0, pv
        do i = 0, pv
          hv(i,k) = c * n * (v(i,j,k,e,2) - vb(i,k,f,2))
        end do
        end do

      case(5,6)
        k = (m - 5) * pv
        n = (m - 5) * 2 - 1
        do j = 0, pv
        do i = 0, pv
          hv(i,j) = c * n * (v(i,j,k,e,3) - vb(i,j,f,3))
        end do
        end do

      end select

      if (pv == pp) then
        hp(:,:,f) = hv
      else
        call InterpolateToPressureSpace(pv, pp, A, hv, hp(:,:,f), w)
      end if

    end do

  end subroutine PressureBC_NormalVelocity_R

  !-----------------------------------------------------------------------------
  !> Build pressure BC from normal velocity conditions -- deformed mesh

  subroutine ressureBC_NormalVelocity_D(boundary, dt, A, n, v, vb, hp)
    class(MeshBoundary_3D), intent(in) :: boundary
    real(RNP),             intent(in)  :: dt
    real(RNP),             intent(in)  :: A(0:,0:)
    real(RNP), contiguous, intent(in)  :: n (0:,0:,:,:,:)
    real(RNP), contiguous, intent(in)  :: v (0:,0:,0:,:,:)
    real(RNP), contiguous, intent(in)  :: vb(0:,0:,:,:)
    real(RNP), contiguous, intent(out) :: hp(0:,0:,:)

    real(RNP), allocatable :: hv(:,:), w(:,:)
    real(RNP) :: c
    integer :: pp, pv
    integer :: f, i, j, k, m

    pv = ubound(vb, 1)
    pp = ubound(hp, 1)

    allocate(hv(0:pv,0:pv), w(0:pp,0:pv))

    c = 1 / dt

    !$omp do
    do f = 1, boundary % n_face

      e = boundary % face(f) % element_id
      m = boundary % face(f) % element_face

      select case(m)

      case(1,2)
        i = (m - 1) * pv
        do k = 0, pv
        do j = 0, pv
          hv(j,k) = c * ( n(j,k,m,e,1) * (v(i,j,k,e,1) - vb(j,k,f,1)) &
                        + n(j,k,m,e,2) * (v(i,j,k,e,2) - vb(j,k,f,2)) &
                        + n(j,k,m,e,3) * (v(i,j,k,e,3) - vb(j,k,f,3)) )
        end do
        end do

      case(3,4)
        j = (m - 3) * pv
        do k = 0, pv
        do i = 0, pv
          hv(i,k) = c * ( n(i,k,m,e,1) * (v(i,j,k,e,1) - vb(i,k,f,1)) &
                        + n(i,k,m,e,2) * (v(i,j,k,e,2) - vb(i,k,f,2)) &
                        + n(i,k,m,e,3) * (v(i,j,k,e,3) - vb(i,k,f,3)) )
        end do
        end do

      case(5,6)
        k = (m - 5) * pv
        do j = 0, pv
        do i = 0, pv
          hv(i,j) = c * ( n(i,j,m,e,1) * (v(i,j,k,e,1) - vb(i,j,f,1)) &
                        + n(i,j,m,e,2) * (v(i,j,k,e,2) - vb(i,j,f,2)) &
                        + n(i,j,m,e,3) * (v(i,j,k,e,3) - vb(i,j,f,3)) )
        end do
        end do

      end select

      if (pv == pp) then
        hp(:,:,f) = hv
      else
        call InterpolateToPressureSpace(pv, pp, A, hv, hp(:,:,f), w)
      end if

    end do

  end subroutine PressureBC_NormalVelocity_D

  !-----------------------------------------------------------------------------
  !> Interpolate face data from velocity into pressure space

  pure subroutine InterpolateToPressureSpace(pv, pp, A, uv, up, w)
    integer,   intent(in)    :: pv            !< polynomial order of velocity
    integer,   intent(in)    :: pp            !< polynomial order of pressure
    real(RNP), intent(in)    :: A (0:pp,0:pv) !< 1D interpolation operator
    real(RNP), intent(in)    :: uv(0:pv,0:pv) !< variable in velocity space
    real(RNP), intent(out)   :: up(0:pp,0:pp) !< variable in pressure space
    real(RNP), intent(inout) :: w (0:pp,0:pv) !< workspace

    real(RNP) :: tmp
    integer   :: i, j, k

    ! direction 1
    do j = 0, pv
    do i = 0, pp
      tmp = 0
      do k = 0, pv
        tmp = tmp + A(i,k) * uv(k,j)
      end do
      w(i,j) = tmp
    end do
    end do

    ! direction 2
    do j = 0, pp
    do i = 0, pp
      tmp = 0
      do k = 0, pv
        tmp = tmp + A(j,k) * w(i,k)
      end do
      up(i,j) = tmp
    end do
    end do

  end subroutine InterpolateToPressureSpace

  !=============================================================================

end submodule MP_PressureSolver
