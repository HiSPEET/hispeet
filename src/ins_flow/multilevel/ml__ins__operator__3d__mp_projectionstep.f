!> summary:  Multilevel Navier-Stokes projection step
!> author:   Joerg Stiller
!> date:     2026/06/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule (ML__INS__Operator__3D) MP_ProjectionStep
  use Constants
  use Array_Assignments
  use TPO__AAA__3D
  use TPO__Div__3D
  use TPO__Grad__3D
  use Trace_Operators__3D
  use Boundary_Variable__3D
  use Parent_To_Child_Interpolation__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Computes the pressure and minimizes the divergence of the given velocity

  module subroutine ProjectionStep(this, tau, bv, u, l_top)
    class(ML_INS_Operator_3D),     intent(in)    :: this
    real(RNP),                     intent(in)    :: tau   !< effective time step
    class(ML_BoundaryVariable_3D), intent(in)    :: bv    !< boundary values
    class(ML_MeshVariable_3D),     intent(inout) :: u     !< solution
    integer,             optional, intent(in)    :: l_top !< top level [auto]

    ! internal variables :::::::::::::::::::::::::::::::::::::::::::::::::::::::

    type(ML_MeshVariable_3D),     allocatable, save :: f, q, w
    type(ML_BoundaryVariable_3D), allocatable, save :: bv_q

    integer :: l_top_
    integer :: l
    logical :: mixed_order

    ! initialization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if (present(l_top)) then
      l_top_ = min(l_top, size(this%ins_op))
    else
      l_top_ = size(this%ins_op)
    end if

    !$omp master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(w, f, q, bv_q)
    call w    % Init(this%ml_op_u, nc = 3, l_top = l_top_)
    call f    % Init(this%ml_op_p, nc = 1, l_top = l_top_)
    call q    % Init(this%ml_op_p, nc = 1, l_top = l_top_)
    call bv_q % Init(this%ml_op_p, nc = 1, l_top = l_top_)
    !$omp end master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    mixed_order = this % ins_op(l_top_) % eop_u % po /=  &
                  this % ins_op(l_top_) % eop_p % po

    ! mass matrix and its inverse
    do l = 1, l_top_
      associate( sem_u  => this % ml_op_u % sem(l)       &
               , mm     => w % level(l) % val(:,:,:,:,1) &
               , mm_inv => w % level(l) % val(:,:,:,:,2) )
        block
          integer :: e

          call sem_u % Get_DG_DiagonalMassMatrix(mm)

          !$omp do
          do e = 1, sem_u % mesh % n_elem
            mm_inv(:,:,:,e) = 1 / mm(:,:,:,e)
          end do
          !$omp end do nowait

        end block
      end associate
    end do

    ! divergence :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    do l = 1, l_top_
      associate( eop_u => this % ins_op(l) % eop_u        &
               , sem_u => this % ml_op_u % sem(l)         &
               , v_    => u % level(l) % val(:,:,:,:,1:3) &
               , div_v => w % level(l) % val(:,:,:,:, 3 ) )
        block
          real(RNP), allocatable, save :: vp(:,:,:,:,:) ! v⁺
          integer :: ne, np

          np = size(v_, 1)
          ne = sem_u % mesh % n_elem

          !$omp master
          allocate(vp(np,np,6,ne,3))
          !$omp end master
          ! no barrier required ;)

          call GetOuterVectorTraces_3D(sem_u%mesh, v_, vp)
          call TPO_Div(eop_u, sem_u, v_, vp, div_v)

          !$omp master
          deallocate(vp)
          !$omp end master
        end block
      end associate
    end do

    ! pressure :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! start values .............................................................

    do l = 1, l_top_
      associate( iop_up => this % ins_op(l) % iop_up     &
               , p_     => u % level(l) % val(:,:,:,:,4) &
               , q_     => q % level(l) % val(:,:,:,:,1) )

        if (mixed_order) then
          call TPO_AAA(iop_up % A, p_, q_)
        else
          call SetArray(q_, p_)
        end if

      end associate
    end do

    ! boundary values ..........................................................

    do l = 1, l_top_
      associate( ins_op => this % ins_op(l)                &
               , mesh   => this % ins_op(l)% mesh          &
               , v      => u % level(l) % val(:,:,:,:,1:3) )
        block
          type(BoundaryVariable_3D), allocatable, save :: bv_u_(:)
          type(BoundaryVariable_3D), allocatable, save :: bv_p_(:), bv_q_(:)
          integer :: b

          !$omp master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          allocate(bv_u_(mesh % n_bound))
          allocate(bv_q_(mesh % n_bound))
          do b = 1, mesh % n_bound
            call bv   % level(l) % var(b) % GetSlice(bv_u_(b), first=1, last=4)
            call bv_q % level(l) % var(b) % GetSlice(bv_q_(b), first=1, last=1)
          end do
          if (mixed_order) then
            allocate(bv_p_(mesh % n_bound))
            do b = 1, mesh % n_bound
              call bv_p_(b) % Init(mesh%boundary(b), ins_op%eop_u%po, nc = 1)
            end do
          end if
          !$omp end master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !$omp barrier

          if (mixed_order) then
            call ins_op % GetPressureBoundaryValues(tau, v, bv_u_, bv_p_, bv_q_)
          else
            call ins_op % GetPressureBoundaryValues(tau, v, bv_u_, bv_q_)
          end if

          !$omp master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          deallocate(bv_u_, bv_q_)
          if (allocated(bv_p_)) deallocate(bv_p_)
          !$omp end master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        end block
      end associate
    end do

    ! RHS ......................................................................

    do l = 1, l_top_
      associate( iop_pu => this % ins_op(l) % iop_pu     &
               , mm     => w % level(l) % val(:,:,:,:,1) &
               , div_v  => w % level(l) % val(:,:,:,:,3) &
               , f_p    => w % level(l) % val(:,:,:,:,3) &
               , f_q    => f % level(l) % val(:,:,:,:,1) )

        f_p = -1/tau * mm * div_v

        if (mixed_order) then
          ! transfer sources to pressure space: w ≈ -1/τ ∫ 𝜑ᵖ f dx
          call TPO_AAA(transpose(iop_pu%A), f_p, f_q)
        else
          f_q = f_p
        end if

      end associate
    end do

    ! solution .................................................................

    call this % ml_solver_p % FAS_MG_Solver( bc     = this%problem%bc_p &
                                           , lambda = ZERO              &
                                           , nu     = ONE               &
                                           , bv     = bv_q              &
                                           , f      = f                 &
                                           , u      = q                 )

    ! transfer .................................................................

    do l = 1, l_top_
      associate( iop_pu => this % ins_op(l) % iop_pu     &
               , p_     => u % level(l) % val(:,:,:,:,4) &
               , q_     => q % level(l) % val(:,:,:,:,1) )

        if (mixed_order) then
          call TPO_AAA(iop_pu % A, q_, p_)
        else
          call SetArray(p_, q_)
        end if

      end associate
    end do

    ! correction :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    do l = 1, l_top_
      associate( eop_u  => this % ins_op(l) % eop_u        &
               , sem_u  => this % ml_op_u % sem(l)         &
               , v_     => u % level(l) % val(:,:,:,:,1:3) &
               , p_     => u % level(l) % val(:,:,:,:, 4 ) &
               , grad_p => w % level(l) % val(:,:,:,:,1:3) )
        block
          real(RNP), allocatable, save :: pp(:,:,:,:) ! p⁺
          integer :: c, na, ne, np

          np = size(v_, 1)
          na = sem_u % mesh % n_elem_active
          ne = sem_u % mesh % n_elem

          !$omp master
          allocate(pp(np,np,6,ne))
          !$omp end master
          ! no barrier required ;)

          call GetOuterTraces_3D(sem_u % mesh, p_, pp)
          call TPO_Grad(eop_u, sem_u, p_, pp, grad_p)

          ! correct velocity: v = v - τ∇p
          do c = 1, 3
            call MergeArrays(ONE, v_(:,:,:,:na,c), -tau, grad_p(:,:,:,:na,c))
          end do

          !$omp master
          deallocate(pp)
          !$omp end master
        end block
      end associate
    end do

    ! update frozen elements :::::::::::::::::::::::::::::::::::::::::::::::::::

    do l = 2, l_top_
      associate( iop_pl  => this % ml_op_u % iop_cf_x(l-1)    &
               , mesh_p  => this % ins_op(l-1)% mesh          &
               , mesh_l  => this % ins_op(l  )% mesh          &
               , v_p     => u % level(l-1) % val(:,:,:,:,1:3) &
               , v_l     => u % level(l  ) % val(:,:,:,:,1:3) &
               , w_l     => w % level(l  ) % val(:,:,:,:,1:3) )
        block
          integer :: c, e

          call ParentToChildInterpolation_3D(mesh_p, mesh_l, iop_pl, v_p, w_l)

          !$omp do collapse(3)
          do c = 1, 3
          do e = mesh_l%n_elem_active+1, mesh_l%n_elem
            v_l(:,:,:,e,c) = w_l(:,:,:,e,c)
          end do
          end do

        end block
      end associate
    end do

    ! finalization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    !$omp master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    deallocate(w, f, q, bv_q)
    !$omp end master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine ProjectionStep

  !=============================================================================

end submodule MP_ProjectionStep
