!> summary:  Generic projection-diffusion step for incompressible flows
!> author:   Joerg Stiller
!> date:     2022/09/24
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @todo
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

  module subroutine ProjectionStep(this, tau, v_0, F_c, F_d, Q, bv_u, u)
    class(INS_TimeIntegrator_3D), intent(in) :: this
    real(RNP), intent(in) :: tau
    real(RNP), contiguous, intent(in) :: v_0(:,:,:,:,:)
    real(RNP), contiguous, intent(in) :: F_c(:,:,:,:,:)
    real(RNP), contiguous, intent(in) :: F_d(:,:,:,:,:)
    real(RNP), contiguous, intent(in) :: Q  (:,:,:,:,:)
    class(SpectralElementBoundaryVariable_3D), intent(inout) :: bv_u(:)
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:,:)

    ! internal variables .......................................................

!### input arguments to be incorporated
integer   :: i_max_p = 4
real(RNP) :: r_red_p = 1E-6
real(RNP) :: r_max_p = 1E-12

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

      np = size(v,1)

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

      associate( eop   => ins_op % eop_v           &
               , met   => ins_op % sem_v % metrics &
               , div_v => w(:,:,:,:,4)             )

        ! generate outer traces of intermediate velocity
        call vm_buf % Transfer(mesh, vm, tag=100)    ! transfer v⁻ from masters
        call ins_op % SetVelocityBC(bv_v, vm)        ! set boundary values
        call vm_buf % Merge(vm)                      ! merge received traces
        call ConvertInnerToOuterTraces(mesh, vm, vp) ! vm → vp = v⁺

        ! divergence of intermediate velocity
        if (mesh % regular) then
          call TPO_Div(eop%w, eop%D, mesh%dx, v, vp, div_v)
        else
          call TPO_Div(eop%w, eop%D, met%Jd, met%Ji, met%a, met%n, v, vp, div_v)
        end if

        ! solve pressure equation
        call ins_op % PressureSolver( tau, bv_v, v, div_v, bv_p, p &
                                    , i_max_p, r_red_p, r_max_p    )

      end associate

      ! pressure correction ....................................................

      associate( eop    => ins_op % eop_v           &
               , met    => ins_op % sem_v % metrics &
               , grad_p => w(:,:,:,:,1:3)           )

        ! generate outer traces of pressure -- preliminary assuming Neumann BC
        call pm_buf % Transfer(mesh, vp, tag=100)    ! transfer p⁻ from masters
        call pm_buf % Merge(pm)                      ! merge received traces
        call ConvertInnerToOuterTraces(mesh, pm, pp) ! pm → pp = p⁺

        ! pressure gradient
        if (mesh % regular) then
          call TPO_Grad(eop%w, eop%D, mesh%dx, p, pp, grad_p)
        else
          call TPO_Grad(eop%w, eop%D, met%Jd, met%Ji, met%a, met%n, p, pp, grad_p)
        end if

        ! correction: v = v - τ∇p
        call MergeArrays(ONE, v, -tau, grad_p, multi=.true.)

      end associate
?
      ! diffusive correction ...................................................

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
