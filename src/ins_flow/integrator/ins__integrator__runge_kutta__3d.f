!> summary:  Runge-Kutta method for incompressible flows
!> author:   Joerg Stiller
!> date:     2023/02/08
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module INS__Integrator__Runge_Kutta__3D
  use Kind_Parameters
  use Constants
  use Array_Assignments
  use IMEX_Runge_Kutta_Method
  use XMPI

  use TPO__Div__3D
  use TPO__Grad__3D
  use Trace_Operators__3D
  use Element_Face_Transfer_Buffer__3D
  use Boundary_Variable__3D

  use INS__Integrator__3D
  use INS__Problem__3D
  use INS__Operator__3D

  implicit none
  private

  public :: INS_Integrator_RungeKutta_3D
  public :: INS_Integrator_RungeKutta_Options_3D

  !-----------------------------------------------------------------------------
  !> Runge-Kutta method for incompressible flows

  type, extends(INS_Integrator_3D) :: INS_Integrator_RungeKutta_3D
    type(IMEX_RK_Method) :: imex_rk !< IMEX Runge-Kutta method
  contains
    procedure, non_overridable :: Init_INS_Integrator_RungeKutta_3D
    procedure :: TimeStep
  end type INS_Integrator_RungeKutta_3D

  ! constructor
  interface INS_Integrator_RungeKutta_3D
    module procedure New_INS_Integrator_RungeKutta_3D
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing Runge-Kutta time-integrator options

  type, extends(INS_IntegratorOptions_3D) :: &
    INS_Integrator_RungeKutta_Options_3D
    integer   :: n_stage = 5 !< number of stages
    integer   :: method  = 1 !< RK method selector, if more than one exist
  contains
    procedure :: Bcast => Bcast_Integrator_RungeKutta_Options
  end type INS_Integrator_RungeKutta_Options_3D

contains

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type INS_Integrator_RungeKutta_3D

  function New_INS_Integrator_RungeKutta_3D(problem, ins_op, opt) result(this)
    class(INS_Problem_3D),                       intent(in) :: problem
    class(INS_Operator_3D),                      intent(in) :: ins_op
    class(INS_Integrator_RungeKutta_Options_3D), intent(in) :: opt
    type(INS_Integrator_RungeKutta_3D) :: this

    call Init_INS_Integrator_RungeKutta_3D(this, problem, ins_op, opt)

  end function New_INS_Integrator_RungeKutta_3D

  !-----------------------------------------------------------------------------
  !> Initialization of a INS_Integrator_RungeKutta_3D object

  subroutine Init_INS_Integrator_RungeKutta_3D(this, problem, ins_op, opt)
    class(INS_Integrator_RungeKutta_3D),         intent(inout) :: this
    class(INS_Problem_3D),                       intent(in)    :: problem
    class(INS_Operator_3D),                      intent(in)    :: ins_op
    class(INS_Integrator_RungeKutta_Options_3D), intent(in)    :: opt

    ! intialize parent type
    call this % Init_INS_Integrator_3D(problem, ins_op, opt)

    ! initialize RK method
    call this % imex_rk % Init_IMEX_RK_Method(opt % n_stage, opt % method)

    this % name = 'Runge-Kutta method: '// trim(this % imex_rk % name)

  end subroutine Init_INS_Integrator_RungeKutta_3D

  !-----------------------------------------------------------------------------
  !> Execution of an Runge-Kutta time step

  subroutine TimeStep(this, t, dt, u, standby)
    class(INS_Integrator_RungeKutta_3D), intent(inout) :: this
    real(RNP),             intent(inout) :: t            !< time t₀ → t
    real(RNP),             intent(in)    :: dt           !< step size ∆t = t-t₀
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:,:) !< u(x,t₀) → u(x,t)
    logical,     optional, intent(in)    :: standby      !< reuse workspace [F]

    ! internal variables .......................................................

    real(RNP), allocatable, save :: inv_mm(:,:,:,:) ! inv diagonal mass matrix
    real(RNP), allocatable, save :: u_i(:,:,:,:,:)  ! stage solution uᵢ
    real(RNP), allocatable, save :: vp (:,:,:,:,:)  ! velocity traces v⁺
    real(RNP), allocatable, save :: sp (:,:,:,:,:)  ! viscous flux traces s⁺
    real(RNP), allocatable, save :: mu (:,:,:,:)    ! variable bulk diffusivity μ
    real(RNP), allocatable, save :: nu (:,:,:,:)    ! variable shear diffusivity ν

    ! stage contributions to RHS and BC
    real(RNP), allocatable, save :: f_c     (:,:,:,:,:,:) ! convection
    real(RNP), allocatable, save :: f_d     (:,:,:,:,:,:) ! diffusion, standard
    real(RNP), allocatable, save :: f_d_rot (:,:,:,:,:,:) ! diffusion, rotational
    real(RNP), allocatable, save :: f_s     (:,:,:,:,:,:) ! sources

    ! boundary points and values
    type(BoundaryVariable_3D), allocatable, save :: bv_x(:), bv_u(:)  &
                                                  , bv_v(:), bv_p(:)  &
                                                  , bv_dp(:), bv_po(:,:)
    ! components of bv_u
    !   'D' - Dirichlet:  [ v₁, v₂, v₃, - ]
    !   'O' - Outflow:    [ - , - , ∆p, p ]

    ! control
    real(RNP), save :: t_0 = -huge(ONE)
    logical,   save :: globally_stiffly_accurate, reuse

    ! auxiliary
    real(RNP), allocatable :: w(:,:,:)
    real(RNP) :: ca, c_ex, c_im, t_i, tau

    integer   :: b, e, d, i, j, np, po

    associate( problem => this % problem           &
             , ins_op  => this % ins_op            &
             , mesh    => this % ins_op % mesh     &
             , sem_u   => this % ins_op % sem_u    &
             , v_0     => u(:,:,:,:,1:3)           &
             , p_0     => u(:,:,:,:, 4 )           &
             , a_im    => this % imex_rk % a_im    &
             , a_ex    => this % imex_rk % a_ex    &
             , b_im    => this % imex_rk % b_im    &
             , b_ex    => this % imex_rk % b_ex    &
             , c       => this % imex_rk % c       &
             , n_stage => this % imex_rk % n_stage &
             )

      !-------------------------------------------------------------------------
      ! initialization

      po = ins_op % eop_u % po
      np = po + 1

      !$omp master

      if (allocated(u_i)) then
        if (any(shape(u_i) /= shape(u))) then
          deallocate(u_i, vp, sp, inv_mm, f_c, f_d, f_d_rot, f_s)
          deallocate(bv_x, bv_u, bv_v, bv_p, bv_dp, bv_po)
          if (allocated(mu)) deallocate(mu)
          if (allocated(nu)) deallocate(nu)
        end if
      end if

      globally_stiffly_accurate = c(n_stage) == ONE                   &
                                  .and. all(b_ex == a_ex(n_stage,:))  &
                                  .and. all(b_im == a_im(n_stage,:))

      reuse = globally_stiffly_accurate .and. t == t_0 .and. allocated(u_i)

      if (.not. allocated(u_i)) then

        allocate( u_i    (np, np, np, mesh % n_elem, 4), source = ZERO )
        allocate( vp     (np, np,  6, mesh % n_elem, 3), source = ZERO )
        allocate( sp     (np, np,  6, mesh % n_elem, 3), source = ZERO )
        allocate( inv_mm (np, np, np, mesh % n_elem   ), source = ZERO )

        if (problem % HasVariableProperties()) then
          allocate( mu(np, np, np, mesh % n_elem), source = this%ins_op%mu_0 )
          allocate( nu(np, np, np, mesh % n_elem), source = ZERO )
        end if

        allocate( f_c     (np, np, np, mesh % n_elem, 3, n_stage), source = ZERO )
        allocate( f_d     (np, np, np, mesh % n_elem, 3, n_stage), source = ZERO )
        allocate( f_d_rot (np, np, np, mesh % n_elem, 3, n_stage), source = ZERO )
        allocate( f_s     (np, np, np, mesh % n_elem, 3, n_stage), source = ZERO )

        allocate( bv_x  (mesh % n_bound) )
        allocate( bv_u  (mesh % n_bound) )
        allocate( bv_v  (mesh % n_bound) )
        allocate( bv_p  (mesh % n_bound) )
        allocate( bv_dp (mesh % n_bound) )
        allocate( bv_po (mesh % n_bound, n_stage) )

        do b = 1, mesh % n_bound
          call bv_x(b) % Init(mesh % boundary(b), po, nc = 3)
          call bv_u(b) % Init(mesh % boundary(b), po, nc = 4)
          call bv_x(b) % Extract(sem_u % metrics % x)
          call bv_u(b) % GetSlice(first=1, last=3, slice = bv_v (b))
          call bv_u(b) % GetSlice(first=4, last=4, slice = bv_p (b))
          call bv_u(b) % GetSlice(first=3, last=3, slice = bv_dp(b))
          do i = 1, n_stage
            call bv_po(b,i) % Init(mesh % boundary(b), po, nc = 1)
          end do
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

      allocate(w(np,np,np))

      !-------------------------------------------------------------------------
      ! stage 1

      t_i = t_0

      if (reuse) then

        call SetArray( f_s(:,:,:,:,:,1), f_s(:,:,:,:,:,n_stage), multi = .true. )
        call SetArray( f_c(:,:,:,:,:,1), f_c(:,:,:,:,:,n_stage), multi = .true. )
        call SetArray( f_d(:,:,:,:,:,1), f_d(:,:,:,:,:,n_stage), multi = .true. )
        call SetArray( f_d_rot (:,:,:,:,:,1)       &
                     , f_d_rot (:,:,:,:,:,n_stage) &
                     , multi = .true.              )

        ! normal stress vector and backflow penalty at ∂Ωᴼ
        do b = 1, mesh % n_bound
          if (problem % bc_v(b) /= 'O') cycle
          call SetArray(bv_po(b,1)%val(:,:,:,1), bv_po(b,n_stage)%val(:,:,:,1))
          call ins_op%GetBackflowPenalty(problem, b, v_0, bv_dp(b)%val(:,:,:,1))
        end do

      else
        associate(v => u(:,:,:,:,1:3))

          ! fetch required boundary values
          do b = 1, mesh % n_bound
            select case(problem % bc_v(b))
            case('D')
              call problem % GetBoundaryValues(b, bv_x(b)%val, t_i, bv_u(b)%val)
            end select
          end do

          ! source term
          call problem % GetExternalSources( sem_u % metrics % x, t_i &
                                           , f_s(:,:,:,:,:,1)         )

          ! variable viscosity
          if (ins_op % HasVariableViscosity()) then
            call ins_op % GetVariableViscosity(t_i, u, mu, nu)
          end if

          ! diffusion term using standard form with extrapolation at ∂Ωᴼ
          call ins_op % GetDiffusionTerm( mu, nu, v, vp, sp &
                                        , f_d(:,:,:,:,:,1)  &
                                        , bv_u              &
                                        , xout = .true.     )

          ! diffusion term using rotational form with extrapolation at ∂Ωᴼ
          call ins_op % GetDiffusionTerm( mu, nu, v, vp, sp    &
                                        , f_d_rot(:,:,:,:,:,1) &
                                        , bv_u                 &
                                        , xout = .true.        &
                                        , form = 2             )

          ! convective RHS
          if (problem % stokes) then
            call SetArray(f_c(:,:,:,:,:,1), ZERO, multi = .true.)
          else
            call ins_op % GetConvectionTerm(v, vp, f_c(:,:,:,:,:,1))
          end if

          !$omp do
          do e = 1, mesh % n_elem
            do d = 1, 3
              f_c    (:,:,:,e,d,1) = inv_mm(:,:,:,e) * f_c    (:,:,:,e,d,1)
              f_d    (:,:,:,e,d,1) = inv_mm(:,:,:,e) * f_d    (:,:,:,e,d,1)
              f_d_rot(:,:,:,e,d,1) = inv_mm(:,:,:,e) * f_d_rot(:,:,:,e,d,1)
            end do
          end do

          ! normal stress vector and backflow penalty at ∂Ωᴼ
          do b = 1, mesh % n_bound
            if (problem % bc_v(b) /= 'O') cycle
            call bv_po(b,1)%MergeNormalTrace(sem_u, cb=ZERO, ct=-ONE, vt=sp)
            call ins_op%GetBackflowPenalty(problem, b, v, bv_dp(b)%val(:,:,:,1))
          end do

        end associate
      end if

      !-------------------------------------------------------------------------
      ! stages 2 to n_stage

      Stages: do i = 2, n_stage
                                                  ! workspace for
        associate( f_d0 => f_d_rot(:,:,:,:,:,i) & ! approximate diffusion term
                 , f    => f_d    (:,:,:,:,:,i) ) ! Stokes RHS

          t_i = t_0 + c(i) * dt
          tau = dt * a_im(i,i)

          ! uᵢ = u₀ ............................................................
          ! not really needed, provides initial approximation for p_i

          call SetArray(u_i, u, multi=.true.)

          ! sources and boundary conditions ....................................

          call problem % GetExternalSources( sem_u % metrics % x, t_i &
                                           , f_s(:,:,:,:,:,i)         )

          do b = 1, mesh % n_bound
            select case(problem % bc_v(b))
            case('D')
              call problem % GetBoundaryValues(b, bv_x(b)%val, t_i, bv_u(b)%val)
            case('O')
              associate(pb => bv_p(b) % val(:,:,:,1))
                ! initialize outlet pressure BC with penalty
                call SetArray(pb, bv_dp(b) % val(:,:,:,1))
                ! merge normal stress vectors from previous stages using
                ! weights identical to the those for viscous terms
                do j = 1, i-1
                  ca = a_ex(i,j) / a_im(i,i)
                  call MergeArrays(ONE, pb, ca, bv_po(b,j) % val(:,:,:,1))
                end do
              end associate
            end select
          end do

          ! RHS and f for projection step ......................................

          !$omp do collapse(2)
          do e = 1, mesh % n_elem
          do d = 1, 3

            c_ex = a_ex(i,1) / a_im(i,i)
            c_im = a_im(i,1) / a_im(i,i)

            f_d0(:,:,:,e,d) = c_ex * f_d_rot(:,:,:,e,d,1) &
                            - c_im * f_d    (:,:,:,e,d,1)

            f(:,:,:,e,d) = c_ex *   f_c (:,:,:,e,d,1) &
                         + c_im * ( f_s (:,:,:,e,d,1) &
                                  + f_d (:,:,:,e,d,1) )

            do j = 2, i-1

              c_ex = a_ex(i,j) / a_im(i,i)
              c_im = a_im(i,j) / a_im(i,i)

              f_d0(:,:,:,e,d) = f_d0(:,:,:,e,d)              &
                              + c_ex * f_d_rot (:,:,:,e,d,j) &
                              - c_im * f_d     (:,:,:,e,d,j)

              f(:,:,:,e,d) = f(:,:,:,e,d)               &
                           + c_ex *   f_c(:,:,:,e,d,j)  &
                           + c_im * ( f_s(:,:,:,e,d,j)  &
                                    + f_d(:,:,:,e,d,j)  )

            end do

            f(:,:,:,e,d) = 1/tau * v_0(:,:,:,e,d) &
                         + f(:,:,:,e,d)           &
                         + f_s(:,:,:,e,d,i)

          end do
          end do

          ! projection-diffusion step ..........................................

          call ins_op % StokesSolver(tau, f, bv_u, mu, nu, u_i, f_d0)

        end associate

        ! RHS contributions ....................................................

        associate(v => u_i(:,:,:,:,1:3))

          ! variable viscosity
          if (ins_op % HasVariableViscosity()) then
            call ins_op % GetVariableViscosity(t_i, u_i, mu, nu)
          end if

          ! diffusion term using standard form with extrapolation at ∂Ωᴼ
          call ins_op % GetDiffusionTerm( mu, nu, v, vp, sp &
                                        , f_d(:,:,:,:,:,i)  &
                                        , bv_u              &
                                        , xout = .true.     )

          ! diffusion term using rotational form with extrapolation at ∂Ωᴼ
          call ins_op % GetDiffusionTerm( mu, nu, v, vp, sp    &
                                        , f_d_rot(:,:,:,:,:,i) &
                                        , bv_u                 &
                                        , xout = .true.        &
                                        , form = 2             )

          ! convection term
          if (problem % stokes) then
            call SetArray(f_c(:,:,:,:,:,i), ZERO, multi = .true.)
          else
            call ins_op % GetConvectionTerm(v, vp, f_c(:,:,:,:,:,i))
          end if

          !$omp do
          do e = 1, mesh % n_elem
            do d = 1, 3
              f_c     (:,:,:,e,d,i) = inv_mm(:,:,:,e) * f_c     (:,:,:,e,d,i)
              f_d     (:,:,:,e,d,i) = inv_mm(:,:,:,e) * f_d     (:,:,:,e,d,i)
              f_d_rot (:,:,:,e,d,i) = inv_mm(:,:,:,e) * f_d_rot (:,:,:,e,d,i)
            end do
          end do

          ! normal stress vector at ∂Ωᴼ
          do b = 1, mesh % n_bound
            if (problem % bc_v(b) /= 'O') cycle
            call bv_po(b,i) % MergeNormalTrace(sem_u, cb=ZERO, ct=-ONE, vt=sp)
          end do

          ! source term, already done ;)

        end associate

      end do Stages

      !-------------------------------------------------------------------------
      ! assembly

      call SetArray(u, u_i, multi=.true.)

      if (.not. globally_stiffly_accurate) then
        associate( v       =>  u  (:,:,:,:,1:3)   &
                 , div_v   =>  f_c(:,:,:,:,1  ,1) &
                 , p       =>  f_c(:,:,:,:,2  ,1) &
                 , pp      =>  vp (:,:,:,:,1)     &
                 , grad_p  =>  f_s(:,:,:,:,1:3,1) )

          !$omp do collapse(2)
          do e = 1, mesh % n_elem
          do d = 1, 3
            do i = 1, n_stage
              v(:,:,:,e,d) = v(:,:,:,e,d)                                 &
                  + dt * (b_ex(i) - a_ex(n_stage,i)) *   f_c(:,:,:,e,d,i) &
                  + dt * (b_im(i) - a_im(n_stage,i)) * ( f_d(:,:,:,e,d,i) &
                                                       + f_s(:,:,:,e,d,i) )
            end do
          end do
          end do

          ! projection .........................................................

          ! update boundary conditions, if required
          if (t_i /= t) then
            do b = 1, mesh % n_bound
              call problem % GetBoundaryValues(b, bv_x(b)%val, t, bv_u(b)%val)
            end do
          end if

          ! velocity divergence
          call GetOuterVectorTraces_3D(mesh, v, vp)    ! vp = v⁺ on Γᴵ and v⁻ on ∂Ω
          call ins_op % ApplyEssentialBC(bv_u, vp, vp) ! vp = v⁺ on ∂Ω
          call TPO_Div(ins_op % eop_u, sem_u, v, vp, div_v)

          ! pressure potential
          call SetArray(p, ZERO)
          call ins_op % PressureSolver(ONE, bv_u, v, div_v, p)

          ! pressure correction
          call GetOuterTraces_3D(mesh, p, pp)
          call TPO_Grad(ins_op % eop_u, sem_u, p, pp, grad_p)
          call MergeArrays(ONE, v, -ONE, grad_p, multi=.true.)

        end associate
      end if

      !-------------------------------------------------------------------------
      ! finalization

      !$omp master
      t_0 = t
      !$omp end master

      if (present(standby)) then
        if (standby) return
      end if

      !$omp master
      deallocate(u_i, vp, sp, inv_mm, f_c, f_d, f_d_rot, f_s)
      deallocate(bv_x, bv_u, bv_v, bv_p, bv_dp, bv_po)
      if (allocated(mu)) deallocate(mu)
      if (allocated(nu)) deallocate(nu)
      !$omp end master

    end associate

  end subroutine TimeStep

  !=============================================================================
  ! TBP of INS_Integrator_RungeKutta_Options_3D

  !-----------------------------------------------------------------------------
  !> MPI broadcasting of Runge-Kutta time-integrator options

  subroutine Bcast_Integrator_RungeKutta_Options(this, root, comm)
    class(INS_Integrator_RungeKutta_Options_3D), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    type(MPI_Request) :: request(6)
    integer :: n

    call this % INS_IntegratorOptions_3D % Bcast(root, comm)

    n = 1
    call XMPI_Ibcast( this % n_stage, root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % method , root, comm, request(n) )

    call MPI_Waitall( n, request, MPI_STATUSES_IGNORE )

  end subroutine Bcast_Integrator_RungeKutta_Options

  !=============================================================================

end module INS__Integrator__Runge_Kutta__3D
