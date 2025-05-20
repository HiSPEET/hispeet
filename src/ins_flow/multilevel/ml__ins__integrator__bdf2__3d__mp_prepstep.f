!> summary:  Preparation for BDF2 time step
!> author:   Joerg Stiller
!> date:     2025/05/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(ML__INS__Integrator__BDF2__3D): MP_PrepStep
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Preparation of a BDF2 time step for one level

  subroutine PrepLevel( ins_op, dt, u, u1, f_c1, f_d1    &
                      , tau, f_d, f, bv_u, mu, nu, first )

    ! arguments ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    class(INS_Operator_3D), intent(in) :: ins_op
    real(RNP), intent(in) :: dt
    real(RNP), contiguous, intent(in) :: u(:,:,:,:,:)
    real(RNP), contiguous, intent(inout) :: u1(:,:,:,:,:)
    real(RNP), contiguous, intent(inout) :: f_c1(:,:,:,:,:)
    real(RNP), contiguous, intent(inout) :: f_d1(:,:,:,:,:)
    real(RNP), intent(out) :: tau
    real(RNP), contiguous, intent(out) :: f_d(:,:,:,:,:)
    real(RNP), contiguous, intent(out) :: f(:,:,:,:,:)
    class(BoundaryVariable_3D), intent(inout) :: bv_u(:)
    real(RNP), contiguous, optional, intent(out) :: mu(:,:,:,:)
    real(RNP), contiguous, optional, intent(out) :: nu(:,:,:,:)
    logical, optional, intent(in) :: first

    ! internal variables :::::::::::::::::::::::::::::::::::::::::::::::::::::::

    real(RNP), allocatable, save :: inv_mm(:,:,:,:) ! inv diagonal mass matrix
    real(RNP), allocatable, save :: v0 (:,:,:,:,:)  ! initial velocity
    real(RNP), allocatable, save :: f_c(:,:,:,:,:)  ! convection term
    real(RNP), allocatable, save :: f_d(:,:,:,:,:)  ! diffusion term
    real(RNP), allocatable, save :: vp (:,:,:,:,:)  ! velocity traces v⁺
    real(RNP), allocatable, save :: sp (:,:,:,:,:)  ! viscous flux traces s⁺

    ! boundary points and values ...............................................

    type(BoundaryVariable_3D), allocatable, save :: bv_x(:)
    type(BoundaryVariable_3D), allocatable, save :: bv_v(:)
    type(BoundaryVariable_3D), allocatable, save :: bv_p(:)
    type(BoundaryVariable_3D), allocatable, save :: bv_dp(:)

    ! IMEX BDF2 coefficients ...................................................

    real(RNP), parameter :: gamma_0 =  3 * HALF
    real(RNP), parameter :: alpha_0 =  2
    real(RNP), parameter :: alpha_1 = -HALF
    real(RNP), parameter :: beta_0  =  2
    real(RNP), parameter :: beta_1  = -1

    ! auxiliary ................................................................

    real(RNP), allocatable :: w(:,:,:)
    real(RNP) :: a0, a1
    integer   :: b, d, e, np, po

    associate( problem => ins_op % problem &
             , mesh    => ins_op % mesh    &
             , sem_u   => ins_op % sem_u   &
             , v       => u(:,:,:,:,1:3)   )

      ! initialization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      ! dimensions ............................................................

      po = ins_op % eop_u % po
      np = po + 1

      !$omp master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! workspace ..............................................................

      allocate( inv_mm  (np, np, np, mesh % n_elem   ), source = ZERO )
      allocate( v0     (np, np, np, mesh % n_elem, 3), source = ZERO )
      allocate( f_c     (np, np, np, mesh % n_elem, 3), source = ZERO )
      allocate( f_d     (np, np, np, mesh % n_elem, 3), source = ZERO )
      allocate( vp      (np, np,  6, mesh % n_elem, 3), source = ZERO )
      allocate( sp      (np, np,  6, mesh % n_elem, 3), source = ZERO )

      ! boundary points and values .............................................

      allocate( bv_x  (mesh % n_bound) )
      allocate( bv_v  (mesh % n_bound) )
      allocate( bv_p  (mesh % n_bound) )
      allocate( bv_dp (mesh % n_bound) )

      do b = 1, mesh % n_bound
        call bv_x(b) % Init(mesh % boundary(b), po, nc = 3)
        call bv_x(b) % Extract(sem_u % metrics % x)
        call bv_u(b) % GetSlice(first=1, last=3, slice = bv_v (b))
        call bv_u(b) % GetSlice(first=4, last=4, slice = bv_p (b))
        call bv_u(b) % GetSlice(first=5, last=5, slice = bv_dp(b))
      end do

      !$omp end master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !$omp barrier

      allocate(w(np,np,np))

      ! inverse diagonal mass matrix ...........................................

      call sem_u % Get_DG_DiagonalMassMatrix(inv_mm)

      !$omp do
      do e = 1, mesh % n_elem
        inv_mm(:,:,:,e) = 1 / inv_mm(:,:,:,e)
      end do
      !$omp end do nowait

      ! initial velocity and effective step size ...............................

      if (first) then
        !$omp do collapse(2)
        do d = 1, 3
        do e = 1, mesh % n_elem
          v0  (:,:,:,e,d) = v(:,:,:,e,d)
          v1(:,:,:,e,d) = v(:,:,:,e,d)
        end do
        end do
        !$omp end do nowait
        tau = dt
      else
        a0 = alpha_0 / gamma_0
        a1 = alpha_1 / gamma_0
        !$omp do collapse(2)
        do d = 1, 3
        do e = 1, mesh % n_elem
          v0  (:,:,:,e,d) = a0 * v(:,:,:,e,d) + a1 * v1(:,:,:,e,d)
          v1(:,:,:,e,d) = v(:,:,:,e,d)
        end do
        end do
        !$omp end do nowait
        tau = dt / gamma_0
      end if

      ! BC at time t₀ ...........................................................

      ! fetch required boundary values
      do b = 1, mesh % n_bound
        select case(problem % bc_v(b))
        case('D')
          call problem % GetBoundaryValues(b, bv_x(b) % val, t_0, bv_u(b) % val)
        end select
      end do

      ! viscosity ..............................................................

      call problem % GetExternalSources(sem_u % metrics % x, t, f)

      if (present(mu) .and. present(nu)) then
        call problem % GetViscosity(sem_u % metrics % x, t, u, nu)
        !$omp do
        do e = 1, mesh % n_elem
          mu(:,:,:,e) = ins_op%mu_0
        end do
        !$omp end do nowait
      end if

      ! RHS contributions ......................................................

      ! diffusion term using rotational form with extrapolation: s⁺ = s⁻ at ∂Ωᴼ
      call ins_op % GetDiffusionTerm( mu, nu, v, vp, sp, f_d, bv_u &
                                    , xout = .true., form = 2      )

      ! convection term
      if (problem % stokes) then
        call SetArray(f_c, ZERO, multi = .true.)
      else
        call ins_op % GetConvectionTerm(v, vp, f_c)
      end if

      if (first) then
        a0 = 1
        a1 = 0
      else
        a0 = beta_0
        a1 = beta_1
      end if

      !$omp do
      do e = 1, mesh % n_elem
        do d = 1, 3

          w = inv_mm(:,:,:,e) * f_c(:,:,:,e,d)
          f_c     (:,:,:,e,d) = a0 * w + a1 * f_c1(:,:,:,e,d)
          f_c1 (:,:,:,e,d) = w

          w = inv_mm(:,:,:,e) * f_d(:,:,:,e,d)
          f_d     (:,:,:,e,d) = a0 * w + a1 * f_d1(:,:,:,e,d)
          f_d1 (:,:,:,e,d) = w

          f(:,:,:,e,d) = 1/tau * v0(:,:,:,e,d) + f_c(:,:,:,e,d) + f(:,:,:,e,d)

        end do
      end do

      ! update boundary conditions .............................................

      do b = 1, mesh % n_bound
        select case(problem % bc_v(b))
        case('D')
          call problem % GetBoundaryValues(b, bv_x(b) % val, t, bv_u(b) % val)
        case('O')
          associate(dp => bv_dp(b) % val(:,:,:,1) )
            call bv_p(b) % MergeNormalTrace(sem_u, cb=ZERO, ct=-ONE, vt=sp)
            call ins_op % GetBackflowPenalty(problem, b, v, dp)
            call MergeArrays(ONE, bv_p(b) % val(:,:,:,1), ONE, dp)
          end associate
        end select
      end do

      ! cleanup ................................................................

      deallocate(inv_mm, v0, f_c, f_d, vp, sp)
      deallocate(bv_x, bv_v, bv_p, bv_dp)

    end associate

  end subroutine PrepLevel

  !=============================================================================

end submodule MP_PrepStep
