!> summary:  Preparation for BDF2 time step
!> author:   Joerg Stiller
!> date:     2025/05/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module INS__Integrator__BDF2__PrepStep__3D
  use Kind_Parameters
  use Boundary_Variable__3D
  use INS__Operator__3D
  implicit none
  private

contains

  !-----------------------------------------------------------------------------
  !> Preparation of a BDF2 time step for one level

  subroutine BDF2_PrepStep_3D( ins_op, t0, dt, u0, u1, f_c1, f_d1 &
                             , tau, f_d, f, bv_u, mu, nu, first   )

    ! arguments ................................................................

    class(INS_Operator_3D), intent(in)    :: ins_op
    real(RNP),              intent(in)    :: t0              !< current time tⁿ
    real(RNP),              intent(in)    :: dt              !< time step ∆t
    real(RNP), contiguous,  intent(in)    :: u0(:,:,:,:,:)   !< uⁿ
    real(RNP), contiguous,  intent(inout) :: u1(:,:,:,:,:)   !< uⁿ⁻¹ → uⁿ
    real(RNP), contiguous,  intent(inout) :: f_c1(:,:,:,:,:) !< f_cⁿ⁻¹ → f_cⁿ
    real(RNP), contiguous,  intent(inout) :: f_d1(:,:,:,:,:) !< f_dⁿ⁻¹ → f_dⁿ
    real(RNP),              intent(out)   :: tau             !< τ
    real(RNP), contiguous,  intent(out)   :: f_d(:,:,:,:,:)  !< approx f_dⁿ⁺¹
    real(RNP), contiguous,  intent(out)   :: f(:,:,:,:,:)    !< approx fⁿ⁺¹
    class(BoundaryVariable_3D), intent(inout) :: bv_u(:)     !< BV at tⁿ⁺¹
    real(RNP), contiguous, optional, intent(out) :: mu(:,:,:,:) !< approx μⁿ⁺¹
    real(RNP), contiguous, optional, intent(out) :: nu(:,:,:,:) !< approx νⁿ⁺¹
    logical, optional, intent(in) :: first !< T/F for Euler/BDF2 [F]

    ! internal variables .......................................................

    real(RNP), allocatable, save :: inv_mm(:,:,:,:) ! inv diagonal mass matrix
    real(RNP), allocatable, save :: f_c(:,:,:,:,:)  ! convection term
    real(RNP), allocatable, save :: u  (:,:,:,:,:)  ! extrapolated variables
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
    real(RNP) :: a0, a1, b0, b1, c0, c1, t
    integer   :: b, d, e, np, po

    associate( problem => ins_op % problem &
             , mesh    => ins_op % mesh    &
             , sem_u   => ins_op % sem_u   &
             , v       => u0(:,:,:,:,1:3)  )

      ! initialization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      t = t0 + dt

      po = ins_op % eop_u % po
      np = po + 1
      nc = problem % nc

      !$omp master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! workspace ..............................................................

      allocate( inv_mm (np, np, np, mesh % n_elem   ), source = ZERO )
      allocate( f_c    (np, np, np, mesh % n_elem, 3), source = ZERO )
      allocate( vp     (np, np,  6, mesh % n_elem, 3), source = ZERO )
      allocate( sp     (np, np,  6, mesh % n_elem, 3), source = ZERO )

      allocate( u, source = u0 )

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
        call bv_u(b) % GetSlice(first=3, last=3, slice = bv_dp(b))
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

      ! effective step size and coefficients ...................................

      if (first) then
        tau = dt
        a0  = 1
        a1  = 0
        b0  = 1
        b1  = 0
      else
        tau = dt / gamma_0
        a0  = alpha_0 / (gamma_0
        a1  = alpha_1 / (gamma_0
        b0  = beta_0
        b1  = beta_1
      end if
      c0 = a0 / tau
      c1 = a1 / tau

      ! viscosity ..............................................................

      if (ins_op % HasVariableViscosity()) then
        ! extrapolated flo variables: u = β₀u₀ + β₁u₁
        call MergeArrays(b0, u, b1, u1, multi=.true.)
        call ins_op % GetVariableViscosity(t, u, mu, nu)
      end if

      ! external source contribution ...........................................

      call problem % GetExternalSources(sem_u % metrics % x, t, f)

      ! convective and diffusive contributions .................................

      ! BC at time t₀
      do b = 1, mesh % n_bound
        select case(problem % bc_v(b))
        case('D')
          ! bv%val(*,1:3) = [ v₁, v₂, v₃ ]
          call problem % GetBoundaryValues(b, bv_x(b) % val, t_0, bv_u(b) % val)
        end select
      end do

      ! diffusion term using rotational form with extrapolation: s⁺ = s⁻ at ∂Ωᴼ
      call ins_op % GetDiffusionTerm( mu, nu, v, vp, sp, f_d, bv_u &
                                    , xout = .true., form = 2      )

      ! convection term
      if (problem % stokes) then
        call SetArray(f_c, ZERO, multi = .true.)
      else
        call ins_op % GetConvectionTerm(v, vp, f_c)
      end if

      !$omp do
      do e = 1, mesh % n_elem
        do d = 1, 3

          ! extrapolated unweighted convection term
          w = inv_mm(:,:,:,e) * f_c(:,:,:,e,d)
          f_c  (:,:,:,e,d) = a0 * w + a1 * f_c1(:,:,:,e,d)
          f_c1 (:,:,:,e,d) = w

          ! extrapolated unweighted diffusion term
          w = inv_mm(:,:,:,e) * f_d(:,:,:,e,d)
          f_d  (:,:,:,e,d) = a0 * w + a1 * f_d1(:,:,:,e,d)
          f_d1 (:,:,:,e,d) = w

          ! unweighted RHS with no diffusion and pressure terms
          f(:,:,:,e,d) = f(:,:,:,e,d)
                       + c0 * u (:,:,:,e,d) &
                       + c1 * u1(:,:,:,e,d) &
                       + f_c(:,:,:,e,d)

        end do
      end do

      ! boundary values at time t ..............................................

      do b = 1, mesh % n_bound
        select case(problem % bc_v(b))
        case('D')
          ! bv%val(*,1:3) = [ v₁, v₂, v₃ ]
          call problem % GetBoundaryValues(b, bv_x(b) % val, t, bv_u(b) % val)
        case('O')
          ! bv%val(*,3:4) = [ ∆p, p ]
          associate(dp => bv_dp(b) % val(:,:,:,1) )
            call bv_p(b) % MergeNormalTrace(sem_u, cb=ZERO, ct=-ONE, vt=sp)
            call ins_op % GetBackflowPenalty(problem, b, v, dp)
            call MergeArrays(ONE, bv_p(b) % val(:,:,:,1), ONE, dp)
          end associate
        end select
      end do

      ! save u₀ to u₁ ..........................................................

      call SetArray(u1, u0, multi = .true.)

      ! cleanup ................................................................

      deallocate(inv_mm, f_c, u, vp, sp)
      deallocate(bv_x, bv_v, bv_p, bv_dp)

    end associate

  end subroutine BDF2_PrepStep_3D

  !=============================================================================

end module INS__Integrator__BDF2__PrepStep__3D
