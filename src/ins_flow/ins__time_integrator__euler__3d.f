!> summary:  Euler method for incompressible flows
!> author:   Joerg Stiller
!> date:     2022/09/09
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module INS__Time_Integrator__Euler__3D
  use Kind_Parameters
  use Constants
  use XMPI

  use Convert_Traces__3D
  use Element_Face_Transfer_Buffer__3D
  use Spectral_Element_Boundary_Variable__3D

  use TPO__Div__3D_D             ! replace with generic module once available
  use TPO__INS_Convection__3D

  use INS__Time_Integrator__3D
  use INS__Problem__3D
  use INS__Operator__3D

  implicit none
  private

  public :: INS_TimeIntegrator_Euler_3D
  public :: INS_TimeIntegrator_Euler_Options_3D

  !-----------------------------------------------------------------------------
  !> Euler method for incompressible flows

  type, extends(INS_TimeIntegrator_3D) :: INS_TimeIntegrator_Euler_3D
  contains
    procedure, non_overridable :: Init_INS_TimeIntegrator_Euler_3D
    procedure :: TimeStep
  end type INS_TimeIntegrator_Euler_3D

  ! constructor
  interface INS_TimeIntegrator_Euler_3D
    module procedure New_INS_TimeIntegrator_Euler_3D
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing Euler time-integrator options (none, so far)

  type, extends(INS_TimeIntegratorOptions_3D) :: &
    INS_TimeIntegrator_Euler_Options_3D
  end type INS_TimeIntegrator_Euler_Options_3D

contains

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type INS_TimeIntegrator_Euler_3D

  function New_INS_TimeIntegrator_Euler_3D(problem, ins_op, opt) result(this)
    class(INS_Problem_3D),  intent(in) :: problem
    class(INS_Operator_3D), intent(in) :: ins_op
    class(INS_TimeIntegrator_Euler_Options_3D), optional, intent(in) :: opt
    type(INS_TimeIntegrator_Euler_3D) :: this

    call Init_INS_TimeIntegrator_Euler_3D(this, problem, ins_op, opt)

  end function New_TimeIntegrator_Euler

  !-----------------------------------------------------------------------------
  !> Initialization of a INS_TimeIntegrator_Euler_3D object

  subroutine Init_INS_TimeIntegrator_Euler_3D(this, problem, ins_op, opt)
    class(INS_TimeIntegrator_Euler_3D), intent(inout) :: this
    class(INS_Problem_3D),              intent(in)    :: problem
    class(INS_Operator_3D),             intent(in)    :: ins_op
    class(INS_TimeIntegrator_Euler_Options_3D), optional, intent(in) :: opt

    ! intialize parent type
    call this % Init_INS_TimeIntegrator_3D(problem, ins_op, opt)
    this % name = 'Euler method'

  end subroutine Init_INS_TimeIntegrator_Euler_3D

  !---------------------------------------------------------------------------
  !> Execution of an Euler time step

  subroutine TimeStep(this, t, dt, u, standby)
    class(INS_TimeIntegrator_Euler_3D), intent(inout) :: this
    real(RNP), intent(inout) :: t  !< time t₀ → t
    real(RNP), intent(in)    :: dt !< step size ∆t = t-t₀
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:,:) !< u(x,t₀) → u(x,t)
    logical, optional, intent(in) :: standby !< reuse workspace T/F [F]
  end subroutine TimeStep

  !---------------------------------------------------------------------------
  !> Generic convection step
  !>
  !> @note
  !> Interface and structure are experimental and not yet optimized ;)

  subroutine ProjectionStep(this, dt, u_0, bv_v, F_c, F_s, u)
    class(INS_TimeIntegrator_Euler_3D), intent(inout) :: this

    real(RNP), intent(in) :: dt
    !< step size ∆t = t-t₀

    real(RNP), contiguous, intent(in) :: u_0(:,:,:,:,:)
    !< u(x,t₀)

    class(SpectralElementBoundaryVariable_3D),intent(in) :: bv_v(:)
    !< bv_v(1:nb): boundary values at time t
    !<  - v         at Dirichlet faces
    !<  - s = n⋅τ   at traction faces   (no yet supported)
    !<  - 0         at free slip faces  (no yet supported)

    real(RNP), contiguous, intent(inout) :: F_c(:,:,:,:,:)
    !< convective contribution to RHS

    real(RNP), contiguous, intent(in) :: F_s(:,:,:,:,:)
    !< source contribution to RHS

    real(RNP), contiguous, intent(inout) :: u(:,:,:,:,:)
    !< u(x,t): (v,p) current → new

    ! internal variables .......................................................

    ! space
    real(RNP), allocatable, save :: inv_M(:,:,:,:) ! inverse diagonal mass matrix
    real(RNP), allocatable, save :: F_d(:,:,:,:,:) ! RHS/diffusion
    real(RNP), allocatable, save :: v_i(:,:,:,:,:) ! intermediate velocity
    real(RNP), allocatable, save :: vm (:,:,:,:,:) ! inner velocity traces v⁻
    real(RNP), allocatable, save :: vp (:,:,:,:,:) ! outer velocity traces v⁺
    real(RNP), allocatable, save :: sp (:,:,:,:,:) ! outer viscous flux traces s⁺
    real(RNP), allocatable, save :: w  (:,:,:,:,:) ! work

    type(ElementFaceTransferBuffer_3D), asynchronous, allocatable, save :: vm_buf

    ! handle for velocity component boundary values
    ! type(SpectralElementBoundaryVariable_3D), allocatable, save :: bv_v(:,:)

    ! options
    logical :: eval_conv   = .true.  ! compute F_c
    logical :: div_with_bc = .true.  ! compute divergence using essential BCs

    ! auxiliary
    real(RNP) :: inv_mm
    integer   :: n_bound, n_elem, n_ghost
    integer   :: po, np
    integer   :: e, d, i, j, k

    associate( problem => this % problem                  &
             , ins_op  => this % ins_op                   &
             , mesh    => this % ins_op % mesh            &
             , n_elem  => this % ins_op % mesh % n_elem   &
             , n_ghost => this % ins_op % mesh % n_ghost  &
             , n_bound => this % ins_op % mesh % n_bound  &
             , v => u(:,:,:,:,1:3)                        &
             )

      ! initialization .........................................................

      po = ins_op % eop_v % po
      np = po + 1

      !$omp master

      allocate( inv_M (np, np, np, n_elem   ), source = ZERO )
      allocate( F_d   (np, np, np, n_elem, 3), source = ZERO )
      allocate( v_i   (np, np, np, n_elem, 3), source = ZERO )
      allocate( w     (np, np, np, n_elem, 4), source = ZERO )

      allocate( vm    (np, np, 6, n_elem + n_ghost, 3), source = ZERO )
      allocate( vp    (np, np, 6, n_elem          , 3), source = ZERO )
      allocate( sp    (np, np, 6, n_elem          , 3), source = ZERO )

      vm_buf = ElementFaceTransferBuffer_3D(mesh, vm)

      !$omp end master
      !$omp barrier

      ! inverse diagonal mass matrix
      call ins_op % sem_v % Get_DG_DiagonalMassMatrix(inv_M)
      inv_M = 1 / inv_M

      ! viscous term: F_d = ∇·τ ................................................
      ! so far ν is constant and boundaries are periodic or have Dirichlet BC

      call ins_op % SetVelocityBC(problem % bc_v, bv_v, vp, sp)

      call ins_op % GetDiffusionTerm( problem % bc_v   &
                                    , problem % nu_ref &
                                    , v, vp, sp, F_d   )
      !$omp do
      do e = 1, ne
        do k = 1, np
        do j = 1, np
        do i = 1, np
          F_d(i,j,k,e,1) = inv_mm(i,j,k,e) * F_d(i,j,k,e,1)
          F_d(i,j,k,e,2) = inv_mm(i,j,k,e) * F_d(i,j,k,e,2)
          F_d(i,j,k,e,3) = inv_mm(i,j,k,e) * F_d(i,j,k,e,3)
        end do
        end do
        end do
      end do

      ! convective term: F_c = -∇·vv ...........................................

      if (eval_conv) then
        associate(metrics => ins_op % sem_q % metrics)

          call TPO_INS_Convection( nv   = ins_op % eop_v % po + 1  &
                                 , nq   = ins_op % sop_q % po + 1  &
                                 , ne   = n_elem                   &
                                 , D_v  = ins_op % eop_v  % D      &
                                 , I_vq = ins_op % iop_vq % A      &
                                 , w_q  = ins_op % sop_q  % w      &
                                 , Jd_q = metrics % Jd             &
                                 , Ji_q = metrics % Ji             &
                                 , a_q  = metrics % a              &
                                 , n_q  = metrics % n              &
                                 , v    = v  (:,:,:,:,1:3)         &
                                 , vp   = vp (:,:,:,:,1:3)         &
                                 , F_c  = F_c(:,:,:,:,1:3)         )
          !$omp do
          do e = 1, ne
            do k = 1, np
            do j = 1, np
            do i = 1, np
              F_c(i,j,k,e,1) = inv_mm(i,j,k) * F_c(i,j,k,e,1)
              F_c(i,j,k,e,2) = inv_mm(i,j,k) * F_c(i,j,k,e,2)
              F_c(i,j,k,e,3) = inv_mm(i,j,k) * F_c(i,j,k,e,3)
            end do
            end do
            end do
          end do
        end associate
      end if

      ! extrapolation step, neglecting pressure ................................

      !$omp do collapse(2)
      do d = 1, 3
      do e = 1, ne

        ! intermediate velocity
        do k = 1, np
        do j = 1, np
        do i = 1, np
            v_i(i,j,k,e,d) = u_0(i,j,k,e,d) + dt * ( F_c(i,j,k,e,d) &
                                                   + F_d(i,j,k,e,d) &
                                                   + F_s(i,j,k,e,d) )
        end do
        end do
        end do

        ! inner traces: vm = v_i⁻
        do k = 1, np
        do j = 1, np
          vm(j,k,1,e,d) = v_i( 1,j,k,e,d)
          vm(j,k,2,e,d) = v_i(np,j,k,e,d)
        end do
        end do
        do k = 1, np
        do i = 1, np
          vm(i,k,3,e,d) = v_i(i, 1,k,e,d)
          vm(i,k,4,e,d) = v_i(i,np,k,e,d)
        end do
        end do
        do j = 1, np
        do i = 1, np
          vm(i,j,5,e,d) = v_i(i,j, 1,e,d)
          vm(i,j,6,e,d) = v_i(i,j,np,e,d)
        end do
        end do

      end do
      end do

      ! generate outer traces ..................................................

      ! start transfer of inner traces
      call vm_buf % Transfer(mesh, vm, tag=100)

      ! optionally inject the boundary values
      if (div_with_bc) then
        call ins_op % SetVelocityBC(problem % bc_v, bv_v, vm)
      end if

      ! finalize transfer and merge remote inner traces
      call vm_buf % Merge(vm)

      ! convert inner to outer traces such that vp = v_i⁺
      call ConvertTraces(mesh, vm, vp)

      ! divergence and pressure .................................................

      ! divergence of intermediate velocity
      call TPO_Div( Ms = ins_op % eop_v % w            &
                  , Ds = ins_op % eop_v % D            &
                  , Jd = ins_op % sem_v % metrics % Jd &
                  , Ji = ins_op % sem_v % metrics % Ji &
                  , a  = ins_op % sem_v % metrics % a  &
                  , n  = ins_op % sem_v % metrics % n  &
                  , u  = v (:,:,:,:,1:3)               &
                  , up = vp(:,:,:,:,1:3)               &
                  , v  = w (:,:,:,:,4)                 )

      ! Review INS_Operator_3D
      !   - availability of viscous diffusion term
      !   - availability of pressure term / operator / solver

      ! roadmap/wishlist/questions
      !   - derive the right pressure boundary conditions
      !   - compute pressure BC and set pressure boundary variable
      !   - solve pressure equation
      !   - verify result

      ! finalization ...........................................................

      !$omp master
      deallocate(inv_M, vp, sp, F_d, v_i, w, bv_v)
      deallocate(vm_buf)
      !$omp end master

    end associate

  end subroutine ProjectionStep

  !=============================================================================

end module INS__Time_Integrator__Euler__3D
