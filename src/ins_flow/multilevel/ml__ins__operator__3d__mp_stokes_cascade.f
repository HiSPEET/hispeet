!> summary:  Cascade start for Stokes part in semi-implicit INS solvers
!> author:   Joerg Stiller
!> date:     2025/05/19
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule (ML__INS__Operator__3D) MP_Stokes_Cascade
  use Parent_To_Child_Interpolation__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Cascade start procedure for the Stokes part

  module subroutine Stokes_Cascade(this, tau, mu, nu, bv, f, u)
    class(ML_INS_Operator_3D), intent(in) :: this
    real(RNP), intent(in) :: tau
      !< effective time step
    class(ML_MeshVariable_3D), optional, intent(in) :: mu
      !< bulk viscosity
    class(ML_MeshVariable_3D), optional, intent(in) :: nu
      !< shear viscosity
    class(ML_BoundaryVariable_3D), intent(in) :: bv
      !< boundary values
    class(ML_MeshVariable_3D), intent(inout) :: f
      !< unweighted RHS: f = v₀/τ + f_c + f_s + ...
    class(ML_MeshVariable_3D), intent(inout) :: u
      !< solution

    integer :: l

    do l = 1, size(u%level) - 1

      if (present(mu) .and. present(nu)) then
        call this % ins_op(l) % StokesSolver( tau                          &
                                            , f  % level(l)%val            &
                                            , bv % level(l)%var            &
                                            , mu % level(l)%val(:,:,:,:,1) &
                                            , nu % level(l)%val(:,:,:,:,1) &
                                            , u  % level(l)%val            )
      else
        call this % ins_op(l) % StokesSolver( tau                    &
                                            , f  = f  % level(l)%val &
                                            , bv = bv % level(l)%var &
                                            , u  = u  % level(l)%val )
      end if

      call ParentToChildInterpolation_3D                   &
               ( parent = this % ml_op_u % sem(l  ) % mesh &
               , child  = this % ml_op_u % sem(l+1) % mesh &
               , iop    = this % ml_op_u % iop_cf_x(l)     &
               , v_p    = u % level(l  ) % val             &
               , v_c    = u % level(l+1) % val             )

    end do

  end subroutine Stokes_Cascade

  !=============================================================================

end submodule MP_Stokes_Cascade
