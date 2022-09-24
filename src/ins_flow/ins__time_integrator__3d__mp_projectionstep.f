!> summary:  Generic projection-diffusion step for incompressible flows
!> author:   Joerg Stiller
!> date:     2022/09/24
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @todo
!>   -  add standby mode
!===============================================================================

submodule (INS__Time_Integrator__3D) MP_ProjectionStep
  implicit none

contains

  !-----------------------------------------------------------------------------
  !>

  module subroutine ProjectionStep(this, tau, v_0, F_c, F_d, Q, bv_u, u)
    class(INS_TimeIntegratorOptions_3D), intent(in) :: this
    real(RNP), intent(in) :: tau
    real(RNP), contiguous, intent(in) :: v_0(:,:,:,:,:)
    real(RNP), contiguous, intent(in) :: F_c(:,:,:,:,:)
    real(RNP), contiguous, intent(in) :: F_d(:,:,:,:,:)
    real(RNP), contiguous, intent(in) :: Q  (:,:,:,:,:)
    class(SpectralElementBoundaryVariable_3D), intent(inout) :: bv_u(:)
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:,:)

    type(SpectralElementBoundaryVariable_3D), allocatable, save :: bv_v(:)
    type(SpectralElementBoundaryVariable_3D), allocatable, save :: bv_p(:)

    integer :: np
    integer :: d, e, i, j, k

    associate( ins_op => this % ins_op        &
             , mesh   => this % ins_op % mesh &
             , v      => u(:,:,:,:,1:3)       &
             , p      => u(:,:,:,:,4)         )

      ! Initialization .........................................................

      !$omp master

      ! handles for velocity and pressure boundary values, based on pointers
      allocate(bv_v(mesh % n_bound), bv_p(mesh % n_bound))
      call bv_u % GetSlice(bv_v, first = 1, last = 3)
      call bv_u % GetSlice(bv_p, first = 4, last = 4)

      !$omp end master
      !$omp barrier

      np = ???

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

      ! start transfer of traces from masters to ghosts
      call vm_buf % Transfer(mesh, vm, tag=100)

      ! inject velocity boundary conditions into boundary traces
      call ins_op % SetVelocityBC(ins_op % bc_v, bv_v, vm)

      ! merge received traces in to ghost entries
      call vm_buf % Merge(vm)

      ! transform inner to outer traces: v'⁻ → v'⁺
      call ConvertInnerToOuterTraces(mesh, vm, vp)

      ! divergence of intermediate velocity
      associate( Ms => ins_op % eop_v % w            &
               , Ds => ins_op % eop_v % D            &
               , Jd => ins_op % sem_v % metrics % Jd &
               , Ji => ins_op % sem_v % metrics % Ji &
               , a  => ins_op % sem_v % metrics % a  &
               , n  => ins_op % sem_v % metrics % n  )

        if (mesh % regular) then
          call TPO_Div(Ms, Ds, mesh%dx, v, vp, w(:,:,:,:,4) ??? )
        else
          call TPO_Div(Ms, Ds, Jd, Ji, a, v, vp, w(:,:,:,:,4) ??? )
        end if
      end associate

      ! solve pressure equation
      call ins_op % PressureSolver( tau, bv_v, v, f ???, bv_p, p &
                                  , i_max ???, r_red ???, r_max ???, ni ??? )

      ! pressure correction ....................................................

      ! diffusive correction ...................................................

      ! cleanup ................................................................

      !$omp master
      deallocate(bv_v, bv_p)
      !$omp end master

    end associate

  end subroutine ProjectionStep

  !=============================================================================

end submodule MP_ProjectionStep
