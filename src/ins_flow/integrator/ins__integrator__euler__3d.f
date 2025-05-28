!> summary:  Euler method for incompressible flows
!> author:   Joerg Stiller
!> date:     2022/09/09
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @note
!>   - The present version implements the IMEX Euler method
!===============================================================================

module INS__Integrator__Euler__3D
  use Kind_Parameters
  use Constants
  use Array_Assignments
  use XMPI

  use Trace_Operators__3D
  use Element_Face_Transfer_Buffer__3D
  use Boundary_Variable__3D

  use INS__Integrator__3D
  use INS__Problem__3D
  use INS__Operator__3D

  implicit none
  private

  public :: INS_Integrator_Euler_3D
  public :: INS_Integrator_Euler_Options_3D

  !-----------------------------------------------------------------------------
  !> Euler method for incompressible flows

  type, extends(INS_Integrator_3D) :: INS_Integrator_Euler_3D
  contains
    procedure, non_overridable :: Init_INS_Integrator_Euler_3D
    procedure :: TimeStep
  end type INS_Integrator_Euler_3D

  ! constructor
  interface INS_Integrator_Euler_3D
    module procedure New_INS_Integrator_Euler_3D
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing Euler time-integrator options (none, so far)

  type, extends(INS_IntegratorOptions_3D) :: &
    INS_Integrator_Euler_Options_3D
  contains
    procedure :: Bcast => Bcast_INS_Integrator_Euler_Options
  end type INS_Integrator_Euler_Options_3D

contains

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type INS_Integrator_Euler_3D

  function New_INS_Integrator_Euler_3D(problem, ins_op, opt) result(this)
    class(INS_Problem_3D),                  intent(in) :: problem
    class(INS_Operator_3D),                 intent(in) :: ins_op
    class(INS_Integrator_Euler_Options_3D), intent(in) :: opt
    type(INS_Integrator_Euler_3D) :: this

    call Init_INS_Integrator_Euler_3D(this, problem, ins_op, opt)

  end function New_INS_Integrator_Euler_3D

  !-----------------------------------------------------------------------------
  !> Initialization of a INS_Integrator_Euler_3D object

  subroutine Init_INS_Integrator_Euler_3D(this, problem, ins_op, opt)
    class(INS_Integrator_Euler_3D),         intent(inout) :: this
    class(INS_Problem_3D),                  intent(in)    :: problem
    class(INS_Operator_3D),                 intent(in)    :: ins_op
    class(INS_Integrator_Euler_Options_3D), intent(in)    :: opt

    ! intialize parent type
    call this % Init_INS_Integrator_3D(problem, ins_op, opt)

    this % name = 'Euler method'

  end subroutine Init_INS_Integrator_Euler_3D

  !-----------------------------------------------------------------------------
  !> Execution of an Euler time step

  subroutine TimeStep(this, t, dt, u, standby)
    class(INS_Integrator_Euler_3D), intent(inout) :: this
    real(RNP),             intent(inout) :: t            !< time t₀ → t
    real(RNP),             intent(in)    :: dt           !< step size ∆t = t-t₀
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:,:) !< u(x,t₀) → u(x,t)
    logical,     optional, intent(in)    :: standby      !< reuse workspace [F]

    ! internal variables .......................................................

    real(RNP), allocatable, save :: inv_mm(:,:,:,:) ! inverse diagonal mass matrix
    real(RNP), allocatable, save :: f_c(:,:,:,:,:)  ! convection term
    real(RNP), allocatable, save :: f_d(:,:,:,:,:)  ! viscous diffusion term
    real(RNP), allocatable, save :: f  (:,:,:,:,:)  ! source term / RHS
    real(RNP), allocatable, save :: vp (:,:,:,:,:)  ! outer velocity traces v⁺
    real(RNP), allocatable, save :: sp (:,:,:,:,:)  ! outer viscous flux traces s⁺
    real(RNP), allocatable, save :: mu (:,:,:,:)    ! variable bulk diffusivity μ
    real(RNP), allocatable, save :: nu (:,:,:,:)    ! variable shear diffusivity ν

    ! boundary points and values
    type(BoundaryVariable_3D), allocatable, save :: bv_x(:), bv_u(:)
    type(BoundaryVariable_3D), allocatable, save :: bv_v(:), bv_p(:), bv_dp(:)
    ! components of bv_u
    !   'D' - Dirichlet:  [ v₁, v₂, v₃, - ]
    !   'O' - Outflow:    [ - , - , ∆p, p ]

    ! control
    real(RNP), save :: t_0 = -huge(ONE)
    logical,   save :: first

    ! auxiliary
    integer :: b, e, d, np, po

    associate( problem => this % problem        &
             , ins_op  => this % ins_op         &
             , mesh    => this % ins_op % mesh  &
             , sem_u   => this % ins_op % sem_u &
             , v       => u(:,:,:,:,1:3)        )

      ! initialization .........................................................

      po = ins_op % eop_u % po
      np = po + 1

      !$omp master

      if (allocated(f)) then
        if (any(shape(f(:,:,:,:,1:3)) /= shape(v))) then
          deallocate(inv_mm, f_c, f_d, f, vp, sp)
          deallocate(bv_x, bv_u, bv_v, bv_p, bv_dp)
          if (allocated(mu)) deallocate(mu)
          if (allocated(nu)) deallocate(nu)
        end if
      end if

      if (allocated(f)) then

        first = abs(t_0 - t) > epsilon(ONE)

      else

        first = .true.

        allocate( inv_mm (np, np, np, mesh % n_elem   ), source = ZERO )
        allocate( f_c    (np, np, np, mesh % n_elem, 3), source = ZERO )
        allocate( f_d    (np, np, np, mesh % n_elem, 3), source = ZERO )
        allocate( f      (np, np, np, mesh % n_elem, 3), source = ZERO )
        allocate( vp     (np, np,  6, mesh % n_elem, 3), source = ZERO )
        allocate( sp     (np, np,  6, mesh % n_elem, 3), source = ZERO )

        if (problem % HasVariableProperties()) then
          allocate( mu(np, np, np, mesh % n_elem), source = this%ins_op%mu_0 )
          allocate( nu(np, np, np, mesh % n_elem), source = ZERO )
        end if

        allocate(bv_x (mesh % n_bound) )
        allocate(bv_u (mesh % n_bound) )
        allocate(bv_v (mesh % n_bound) )
        allocate(bv_p (mesh % n_bound) )
        allocate(bv_dp(mesh % n_bound) )
        do b = 1, mesh % n_bound
          call bv_x(b) % Init(mesh % boundary(b), po, nc = 3)
          call bv_u(b) % Init(mesh % boundary(b), po, nc = 4)
          call bv_x(b) % Extract(sem_u % metrics % x)
          call bv_u(b) % GetSlice(first=1, last=3, slice = bv_v (b))
          call bv_u(b) % GetSlice(first=4, last=4, slice = bv_p (b))
          call bv_u(b) % GetSlice(first=3, last=3, slice = bv_dp(b))
        end do

      end if

      t_0 = t
      t   = t + dt

      !$omp end master
      !$omp barrier

      ! inverse diagonal mass matrix
      call sem_u % Get_DG_DiagonalMassMatrix(inv_mm)
      !$omp workshare
      inv_mm = 1 / inv_mm
      !$omp end workshare nowait

      ! BC at time t₀ ..........................................................

      do b = 1, mesh % n_bound
        select case(problem % bc_v(b))
        case('D')
          call problem % GetBoundaryValues(b, bv_x(b)%val, t_0, bv_u(b)%val)
        end select
      end do

      ! RHS ....................................................................

      ! external sources, f_s
      call problem % GetExternalSources(sem_u % metrics % x, t, f)

      if (ins_op % HasVariableViscosity()) then
        call ins_op % GetVariableViscosity(t, u, mu, nu)
      end if

      ! rotational diffusion term F_d0 with extrapolation: s⁺ = s⁻ at ∂Ωᴼ
      call ins_op % GetDiffusionTerm( mu, nu, v, vp, sp, f_d, bv_u &
                                    , xout = .true., form = 2      )

      ! convection term F_c
      if (problem % stokes) then
        call SetArray(f_c, ZERO, multi = .true.)
      else
        call ins_op % GetConvectionTerm(v, vp, f_c)
      end if

      !$omp do
      do e = 1, mesh % n_elem
        do d = 1, 3
          f_c(:,:,:,e,d) = inv_mm(:,:,:,e) * f_c(:,:,:,e,d) ! f_c  = M⁻¹ F_c
          f_d(:,:,:,e,d) = inv_mm(:,:,:,e) * f_d(:,:,:,e,d) ! f_d0 = M⁻¹ F_d0
          f(:,:,:,e,d) = 1/dt * v(:,:,:,e,d) + f_c(:,:,:,e,d) + f(:,:,:,e,d)
        end do
      end do

      ! boundary conditions at time t = t₀ + ∆t.................................

      do b = 1, mesh % n_bound
        select case(problem % bc_v(b))
        case('D')
          ! v = vᵇ  on Ωᴰ
          call problem % GetBoundaryValues(b, bv_x(b) % val, t, bv_u(b) % val)
        case('O')
          associate(pb => bv_p(b) % val(:,:,:,1), dp => bv_dp(b) % val(:,:,:,1))
            ! pᵇ = -n⋅s⁺ + ∆pᵇ, ∆pᵇ = -E(v,n)
            call bv_p(b) % MergeNormalTrace(sem_u, cb = ZERO, ct = -ONE, vt = sp)
            call ins_op % GetBackflowPenalty(problem, b, v, dp)
            call MergeArrays(ONE, pb, ONE, dp)
          end associate
        end select
      end do

      ! extrapolation-projection-diffusion step ................................

      call ins_op % StokesSolver(dt, f, bv_u, mu, nu, u, f_d0 = f_d)

      ! cleanup ................................................................

      if (present(standby)) then
        if (standby) return
      end if

      !$omp master
      deallocate(inv_mm, f_c, f_d, f, vp, sp)
      deallocate(bv_x, bv_u, bv_v, bv_p, bv_dp)
      if (allocated(mu)) deallocate(mu)
      if (allocated(nu)) deallocate(nu)
      !$omp end master

    end associate

  end subroutine TimeStep

  !=============================================================================
  ! TBP of INS_Integrator_Euler_Options_3D

  !-----------------------------------------------------------------------------
  !> MPI broadcasting of Euler time-integrator options

  subroutine Bcast_INS_Integrator_Euler_Options(this, root, comm)
    class(INS_Integrator_Euler_Options_3D), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    call this % INS_IntegratorOptions_3D % Bcast(root, comm)

  end subroutine Bcast_INS_Integrator_Euler_Options

  !=============================================================================

end module INS__Integrator__Euler__3D
