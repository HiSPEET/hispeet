!------------------------------------------------------------------------------!
! This file is part of HiSPEET: High-order Spectral Element Techniques         !
!                                                                              !
! Copyright (C) 2026 by the HiSPEET authors and the                            !
! Chair of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany                 !
!                                                                              !
! HiSPEET is free software: you can redistribute it and/or modify              !
! it under the terms of the GNU General Public License as published by         !
! the Free Software Foundation, either version 3 of the License, or            !
! (at your option) any later version.                                          !
!                                                                              !
! HiSPEET is distributed in the hope that it will be useful,                   !
! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                         !
! See the GNU General Public License for more details.                         !
!                                                                              !
! You should have received a copy of the GNU General Public License            !
! along with HiSPEET. If not, see <http://www.gnu.org/licenses/>.              !
!------------------------------------------------------------------------------!

!> summary:  FAS-MG residual
!> author:   Joerg Stiller
!> date:     2026/08/01
!===============================================================================

submodule(ML__DG__Elliptic_Solver__3D) MP_FAS_MG_Residual
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> FAS-MG residual with constant diffusivity
  !>
  !> Use `l_top` to specify a top level lower than `size(this%ml_op%sem)`

  module subroutine FAS_MG_Residual_C(this, bc, lambda, nu, bv, f, u, r, l_top)
    class(ML_DG_EllipticSolver_3D), intent(in)    :: this
    character,                      intent(in)    :: bc(:)
    real(RNP),                      intent(in)    :: lambda
    real(RNP),                      intent(in)    :: nu
    class(ML_BoundaryVariable_3D),  intent(in)    :: bv
    class(ML_MeshVariable_3D),      intent(in)    :: f
    class(ML_MeshVariable_3D),      intent(in)    :: u
    class(ML_MeshVariable_3D),      intent(inout) :: r
    integer,              optional, intent(in)    :: l_top

    call FAS_MG_Residual_X(this, bc, lambda, nu, null(), bv, f, u, r, l_top)

  end subroutine FAS_MG_Residual_C

  !-----------------------------------------------------------------------------
  !> FAS-MG residual with variable diffusivity
  !>
  !> Use `l_top` to specify a top level lower than `size(this%ml_op%sem)`

  module subroutine FAS_MG_Residual_V(this, bc, lambda, nu, bv, f, u, r, l_top)
    class(ML_DG_EllipticSolver_3D), intent(in)    :: this
    character,                      intent(in)    :: bc(:)
    real(RNP),                      intent(in)    :: lambda
    class(ML_MeshVariable_3D),      intent(in)    :: nu
    class(ML_BoundaryVariable_3D),  intent(in)    :: bv
    class(ML_MeshVariable_3D),      intent(in)    :: f
    class(ML_MeshVariable_3D),      intent(in)    :: u
    class(ML_MeshVariable_3D),      intent(inout) :: r
    integer,              optional, intent(in)    :: l_top

    call FAS_MG_Residual_X(this, bc, lambda, null(), nu, bv, f, u, r, l_top)

  end subroutine FAS_MG_Residual_V

  !-----------------------------------------------------------------------------
  !> Generic FAS-MG residual with constant or variable diffusivity
  !>
  !> If the boundary values `bv` and the sources `f` are absent, the result `r`
  !> equals the negative homogeneous FAS operator
  !> Use `l_top` to specify a top level lower than `size(this%ml_op%sem)`.

  module subroutine FAS_MG_Residual_X( this, bc, lambda, nu_0, nu_v, bv, f &
                                     , u, r, l_top)
    class(ML_DG_EllipticSolver_3D),                  intent(in)    :: this
    character,                                       intent(in)    :: bc(:)
    real(RNP),                                       intent(in)    :: lambda
    real(RNP),                             optional, intent(in)    :: nu_0
    class(ML_MeshVariable_3D),             optional, intent(in)    :: nu_v
    class(ML_BoundaryVariable_3D), target, optional, intent(in)    :: bv
    class(ML_MeshVariable_3D),             optional, intent(in)    :: f
    class(ML_MeshVariable_3D),                       intent(in)    :: u
    class(ML_MeshVariable_3D),                       intent(inout) :: r
    integer,                               optional, intent(in)    :: l_top

    type(ML_MeshVariable_3D), allocatable, save :: g, v
    type(BoundaryVariable_3D), pointer, save :: bv_l(:) => null()
    type(BoundaryVariable_3D), pointer, save :: bv_p(:) => null()

    integer :: e, l, l_top_

    character(len=:), allocatable :: prefix
    logical :: logging

    ! start logging ............................................................

    if (log_level > 0) then
      associate(proc => this % ml_op % sem(1) % mesh % proc)
        logging = proc == 0 .or. log_level > 1
        prefix  = LoggingPrefix('FAS_MG_Residual_X', proc)
      end associate
    else
      logging = .false.
    end if

    if (logging) then
      print '(2A)', prefix, 'start'
    end if

    ! initialization ...........................................................

    if (present(l_top)) then
      l_top_ = min(l_top, size(this%ml_op%sem))
    else
      l_top_ = size(this%ml_op%sem)
    end if

    !$omp master
    allocate(g, v)
    call g % Init(this%ml_op, nc = 1, l_top = l_top_)
    call v % Init(this%ml_op, nc = 1, l_top = l_top_)
    !$omp end master

    if (present(f)) then
      !$omp barrier !! g must be initialized !!
      call ML_SetArray_3D(g, f, l_top = l_top_)
    end if

    do l = l_top_, 1, -1
      associate( g_l => g % level(l) % val(:,:,:,:,1) &
               , u_l => u % level(l) % val(:,:,:,:,1) &
               , r_l => r % level(l) % val(:,:,:,:,1) )

        !$omp master
        if (present(bv)) then
          bv_l => bv % level(l) % var
        end if
        !$omp end master

        ! residual .............................................................

        if (present(nu_0)) then
          call this % Residual(l, bc, lambda, nu_0, g_l, bv_l, u_l, r_l)
        else
          associate(nu_l => nu_v % level(l) % val(:,:,:,:,1))
            call this % Residual(l, bc, lambda, nu_l, g_l, bv_l, u_l, r_l)
          end associate
        end if

        if (l == 1) exit

        ! parent FAS-RHS .......................................................

        associate( mesh_l => this % ml_op % sem(l  ) % mesh  &
                 , mesh_p => this % ml_op % sem(l-1) % mesh  &
                 , pop_lp => this % ml_op % pop_fc_x(l)      &
                 , iop_lp => this % ml_op % iop_fc_x(l)      &
                 , iop_pl => this % ml_op % iop_cf_x(l-1)    &
                 , ell_p  => this % elliptic_op(l-1)         &
                 , g_p    => g % level(l-1) % val(:,:,:,:,1) &
                 , r_p    => r % level(l-1) % val(:,:,:,:,1) &
                 , u_p    => u % level(l-1) % val(:,:,:,:,1) &
                 , v_p    => v % level(l-1) % val(:,:,:,:,1) )

          !$omp master
          if (present(bv)) then
            bv_p => bv % level(l-1) % var
          end if
          !$omp end master

          ! project solution to regularly refined parent elements
          select case(this % fc_projection)
          case('I')
            ! interpolation
            call ChildToParentProjection_3D(mesh_l, mesh_p, iop_lp, u_l, v_p)
          case('P')
            ! L²-projection
            call ChildToParentProjection_3D(mesh_l, mesh_p, pop_lp, u_l, v_p)
          end select

          ! restrict residual
          select case(this % fc_restriction)
          case('C')
            ! canonical restriction
            call ChildToParentRestriction_3D(mesh_l, mesh_p, iop_pl, r_l, r_p)
          case('P')
            ! L²-projection
            call ChildToParentProjection_3D(mesh_l, mesh_p, pop_lp, r_l, r_p)
          end select

          do e = 1, mesh_p % n_elem
            if (mesh_p % element(e) % adaptation % refinement >= 1000) then
              ! residual contribution to FAS-RHS in parent twigs
              g_p(:,:,:,e) = r_p(:,:,:,e)
            else
              ! set v_p to solution in leaves
              v_p(:,:,:,e) = u_p(:,:,:,e)
            end if
          end do

          ! apply parent operator to projected solution
          if (present(nu_0)) then
            call ell_p % Apply(bc, lambda, nu_0, bv_p, v_p, r_p)
          else
            associate(nu_p => nu_v % level(l-1) % val(:,:,:,:,1))
              call ell_p % Apply(bc, lambda, nu_p, bv_p, v_p, r_p)
            end associate
          end if

          ! add operator contribution to parent FAS-RHS
          do e = 1, mesh_p % n_elem
            if (mesh_p % element(e) % adaptation % refinement >= 1000) then
              g_p(:,:,:,e) = g_p(:,:,:,e) + r_p(:,:,:,e)
            end if
          end do

        end associate

      end associate
    end do

    ! finalization .............................................................

    !$omp master
    deallocate(g, v)
    bv_l => null()
    bv_p => null()
    !$omp end master

    ! exit logging .............................................................

    if (logging) then
      print '(2A)', prefix, 'exit'
    end if

  end subroutine FAS_MG_Residual_X

  !=============================================================================

end submodule MP_FAS_MG_Residual
