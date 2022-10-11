!> summary:  Generic projection-diffusion step for incompressible flows
!> author:   Joerg Stiller
!> date:     2022/09/24
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @todo
!>   - validate
!>   - revise interface
!>   - add standby mode
!===============================================================================

submodule (INS__Time_Integrator__3D) MP_ProjectionStep
  use Array_Assignments
  use TPO__Div__3D
  use TPO__Grad__3D
  use Convert_Traces__3D
  use Element_Face_Transfer_Buffer__3D
  use Spectral_Element_Boundary_Variable__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !>

  module subroutine ProjectionStep( this, tau, v_0, F_c, F_d, Q, bv_u, u &
                                  , i_max_p, i_max_v, r_red, r_max       )


    ! arguments ................................................................

    class(INS_TimeIntegrator_3D), intent(in) :: this
    real(RNP), intent(in) :: tau
    real(RNP), contiguous, intent(in) :: v_0(:,:,:,:,:)
    real(RNP), contiguous, intent(in) :: F_c(:,:,:,:,:)
    real(RNP), contiguous, intent(in) :: F_d(:,:,:,:,:)
    real(RNP), contiguous, intent(in) :: Q  (:,:,:,:,:)
    class(SpectralElementBoundaryVariable_3D), intent(inout) :: bv_u(:)
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:,:)

    integer :: i_max_p !< max num iterations of pressure solver
    integer :: i_max_v !< max num iterations of viscous diffusion solver
    real(RNP), optional :: r_red !< minimum L² residual reduction to reach
    real(RNP), optional :: r_max !< maximum L² residual allowed

    ! internal variables .......................................................

    real(RNP), allocatable, save :: pm (:,:,:,:)   ! inner pressure traces p⁻
    real(RNP), allocatable, save :: pp (:,:,:,:)   ! outer pressure traces p⁺
    real(RNP), allocatable, save :: vm (:,:,:,:,:) ! inner velocity traces v⁻
    real(RNP), allocatable, save :: vp (:,:,:,:,:) ! outer velocity traces v⁺
    real(RNP), allocatable, save :: w  (:,:,:,:,:) ! work

    type(ElementFaceTransferBuffer_3D), asynchronous, allocatable, save :: pm_buf
    type(ElementFaceTransferBuffer_3D), asynchronous, allocatable, save :: vm_buf

    type(SpectralElementBoundaryVariable_3D), allocatable, save :: bv_v(:)
    type(SpectralElementBoundaryVariable_3D), allocatable, save :: bv_p(:)

    integer :: np
    integer :: d, e, i, j, k

    associate( ins_op  => this % ins_op                   &
             , mesh    => this % ins_op % mesh            &
             , n_elem  => this % ins_op % mesh % n_elem   &
             , n_ghost => this % ins_op % mesh % n_ghost  &
             , v       => u(:,:,:,:,1:3)                  &
             , p       => u(:,:,:,:,4)                    )

      ! Initialization .........................................................

      np = size(v,1)

      !$omp master

      allocate( pm (np, np,  6, n_elem + n_ghost   ), source = ZERO )
      allocate( pp (np, np,  6, n_elem             ), source = ZERO )
      allocate( vm (np, np,  6, n_elem + n_ghost, 3), source = ZERO )
      allocate( vp (np, np,  6, n_elem          , 3), source = ZERO )
      allocate( w  (np, np, np, n_elem          , 4), source = ZERO )

      vm_buf = ElementFaceTransferBuffer_3D(mesh, vm)
      pm_buf = ElementFaceTransferBuffer_3D(mesh, pm)

      ! handles for velocity and pressure boundary values, based on pointers
      allocate(bv_v(mesh % n_bound), bv_p(mesh % n_bound))
      call bv_u % GetSlice(bv_v, first = 1, last = 3)
      call bv_u % GetSlice(bv_p, first = 4, last = 4)

      !$omp end master
      !$omp barrier

      ! Extrapolation step .....................................................

      ! Computes the intermediate velocity v' = v₀ + τ(Fc + Fd + Q) and extracts
      ! the inner traces v'⁻ needed in the next step.

      !$omp do collapse(2)
      do d = 1, 3
      do e = 1, mesh % n_elem

        ! intermediate velocity v = v'
        do k = 1, np
        do j = 1, np
        do i = 1, np
            v(i,j,k,e,d) = v_0(i,j,k,e,d) + tau * ( F_c(i,j,k,e,d) &
                                                  + F_d(i,j,k,e,d) &
                                                  + Q  (i,j,k,e,d) )
        end do
        end do
        end do

        ! inner traces: vm = v'⁻
        do k = 1, np
        do j = 1, np
          vm(j,k,1,e,d) = v( 1,j,k,e,d)
          vm(j,k,2,e,d) = v(np,j,k,e,d)
        end do
        end do
        do k = 1, np
        do i = 1, np
          vm(i,k,3,e,d) = v(i, 1,k,e,d)
          vm(i,k,4,e,d) = v(i,np,k,e,d)
        end do
        end do
        do j = 1, np
        do i = 1, np
          vm(i,j,5,e,d) = v(i,j, 1,e,d)
          vm(i,j,6,e,d) = v(i,j,np,e,d)
        end do
        end do

      end do
      end do

      ! pressure computation ...................................................

      associate(div_v => w(:,:,:,:,4))

        ! generate outer traces of intermediate velocity
        call vm_buf % Transfer(mesh, vm, tag=100)    ! transfer v⁻ from masters
        call ins_op % SetVelocityBC(bv_v, vm)        ! set boundary values
        call vm_buf % Merge(vm)                      ! merge received traces
        call ConvertInnerToOuterTraces(mesh, vm, vp) ! vm → vp = v⁺

        ! divergence of intermediate velocity
        call TPO_Div(ins_op % eop_v, ins_op % sem_v, v, vp, div_v)

        ! solve pressure equation
        call ins_op % PressureSolver( tau, bv_v, v, div_v, bv_p, p &
                                    , i_max_p, r_red, r_max        )

      end associate

      ! pressure correction ....................................................

      associate(grad_p => w(:,:,:,:,1:3))

        ! generate outer traces of pressure -- preliminary assuming Neumann BC
        call pm_buf % Transfer(mesh, vp, tag=100)    ! transfer p⁻ from masters
        call pm_buf % Merge(pm)                      ! merge received traces
        call ConvertInnerToOuterTraces(mesh, pm, pp) ! pm → pp = p⁺

        ! pressure gradient
        call TPO_Grad(ins_op % eop_v, ins_op % sem_v, p, pp, grad_p)

        ! correction: v = v - τ∇p
        call MergeArrays(ONE, v, -tau, grad_p, multi=.true.)

      end associate

      ! diffusive correction ...................................................

      associate(f => w(:,:,:,:,1:3))

        !$omp do collapse(2)
        do d = 1, 3
        do e = 1, mesh % n_elem
          do k = 1, np
          do j = 1, np
          do i = 1, np
            f(i,j,k,e,d) = 1/tau * v(i,j,k,e,d) - F_d(i,j,k,e,d)
          end do
          end do
          end do
        end do
        end do

        call ins_op % DiffusionSolver(tau, f, bv_v, v, i_max_v, r_red, r_max)

      end associate

      ! cleanup ................................................................

      !$omp master
      deallocate(w)
      deallocate(pm, pp, vm, vp)
      deallocate(pm_buf, vm_buf)
      deallocate(bv_v, bv_p)
      !$omp end master

    end associate

  end subroutine ProjectionStep

  !=============================================================================

end submodule MP_ProjectionStep
