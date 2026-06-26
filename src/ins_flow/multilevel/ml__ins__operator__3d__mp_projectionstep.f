!> summary:  Multilevel projection step
!> author:   Joerg Stiller
!> date:     2026/06/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule (ML__INS__Operator__3D) MP_ProjectionStep
  use TPO__AAA__3D
  use TPO__Div__3D
  use Trace_Operators__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Performs

  module subroutine MP_ProjectionStep(this, tau, bv_u, u)
    class(ML_INS_Operator_3D), intent(in) :: this
    real(RNP), intent(in) :: tau
      !< effective time step
    class(ML_BoundaryVariable_3D), intent(in) :: bv_u
      !< boundary values
    class(ML_MeshVariable_3D), intent(inout) :: u
      !< solution

    ! internal variables :::::::::::::::::::::::::::::::::::::::::::::::::::::::

    type(ML_MeshVariable_3D),     allocatable, save :: w_u, w_q
    type(ML_MeshVariable_3D),     allocatable, save :: f, q
    type(ML_BoundaryVariable_3D), allocatable, save :: bv_q

    integer :: l_top
    integer :: l
    logical :: mixed_order

    ! initialization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    !$omp master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(w_u, f, q, bv_q)
    call w_u  % Init(this%ml_op_u, nc = 3)
    call f    % Init(this%ml_op_p, nc = 1)
    call q    % Init(this%ml_op_p, nc = 1)
    call bv_q % Init(this%ml_op_p, nc = 1)
    !$omp end master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    l_top = size(this % ins_op)

    mixed_order = this % ins_op(l_top) % eop_u % po /=
                  this % ins_op(l_top) % eop_p % po

    ! mass matrix and its inverse
    do l = 1, l_top
      associate( sem_u  => this % ml_op_u % sem_u(l)       &
                 mm_    => w_u % level(l) % val(:,:,:,:,1) &
                 mm_inv => w_u % level(l) % val(:,:,:,:,2) )

        call sem_u % Get_DG_DiagonalMassMatrix(mm)

        !$omp do
        do e = 1, sem_u % mesh % n_elem
          mm_inv(:,:,:,e) = 1 / mm(:,:,:,e)
        end do
        !$omp end do nowait

      end associate
    end do

    ! divergence :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    do l = 1, l_top
      associate( sem_u => this % ml_op_u % sem_u_(l)        &
                 v_    => u   % level(l) % val(:,:,:,:,1:3) &
                 div_v => w_u % level(l) % val(:,:,:,:,3)   )
        block
          real(RNP), allocatable, save :: vp(:,:,:,:,:) ! v⁺
          integer :: na, ne, np

          np = size(v_, 1)
          ne = sem_u % mesh % n_elem

          !$omp master
          allocate(vp(np,np,np,ne,3))
          !$omp end master
          ! no barrier required ;)

          call GetOuterVectorTraces_3D(sem_u%mesh, v_, vp)
          call TPO_Div(sem_u_%eop, sem_u, v_, vp, div_v)

          !$omp master
          deallocate(vp)
          !$omp end master
        end block
      end associate
    end do

    ! pressure :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! start values .............................................................

    do l = 1, l_top
      associate( iop_up => this % ins_op(l) % iop_up     &
                 p_     => u % level(l) % val(:,:,:,:,4) &
                 q_     => q % level(l) % val(:,:,:,:,1) )

        if (mixed_order) then
          call TPO_AAA(iop_up % A, p_, q_)
        else
          call SetArray(q_, p_)
        end if

      end associate
    end do

    ! boundary values ..........................................................

    do l = 1, l_top
      associate( ins_op => this % ins_op(l)                &
                 mesh   => this % ins_op(l)) % mesh        &
                 v      => u % level(l) % val(:,:,:,:,1:3) )
        block
          type(BoundaryVariable_3D), allocatable, save :: bv_u_(:)
          type(BoundaryVariable_3D), allocatable, save :: bv_p_(:), bv_q_(:)
          integer :: b

          !$omp master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          allocate(bv_u_(mesh % n_bound))
          allocate(bv_q_(mesh % n_bound))
          do b = 1, mesh % n_bound
            call bv_u % level(l) % var(b) % GetSlice(bv_u_(b), first=1, last=4)
            call bv_q % level(l) % var(b) % GetSlice(bv_q_(n), first=1, last=1)
          end do
          if (mixed_order) then
            allocate(bv_p_(mesh % n_bound))
            do b = 1, mesh % n_bound
              call bv_p_(b) % Init(mesh%boundary(b), ins_op%eop_u%po, nc = 1)
            end do
          end if
          !$omp end master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !$omp barrier

          if (mixed_order) then
            call ins_op % GetPressureBoundaryValues(tau, v, bv_u_, bv_p_, bv_q_)
          else
            call ins_op % GetPressureBoundaryValues(tau, v, bv_u_, bv_q_)
          end if

          !$omp master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          deallocate(bv_u_, bv_q_)
          if (allocated(bv_p_)) deallocate(bv_p_)
          !$omp end master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        end block
      end associate
    end do

    ! RHS ......................................................................

    do l = 1, l_top
      associate( ins_op  => this % ins_op(l)                 &
                 mm      => w_u  % level(l) % val(:,:,:,:,1) &
                 div_v   => w_u  % level(l) % val(:,:,:,:,3) &
                 f_p     => w_u  % level(l) % val(:,:,:,:,3) &
                 f_q     => f    % level(l) % val(:,:,:,:,1) )

        f_p = -1/tau * mm * div_v

        if (mixed_order) then
          ! transfer sources to pressure space: w ≈ -1/τ ∫ 𝜑ᵖ f dx
          call TPO_AAA(transpose(ins_op%iop_pu%A), f_p, f_q)
        else
          f_q = f_p
        end if

      end associate
    end do

    ! solution ...............................................................


    ! correction :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! finalization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    !$omp master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    deallocate(w_u, f, q, bv_q)
    !$omp end master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine MP_ProjectionStep

  !=============================================================================

end submodule MP_ProjectionStep
