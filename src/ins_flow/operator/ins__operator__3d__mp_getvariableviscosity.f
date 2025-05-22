!> summary:  Provision of variable viscosity coefficients
!> author:   Joerg Stiller
!> date:     2025/05/21
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(INS__Operator__3D) MP_GetVariableViscosity
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Variable viscosity coefficients

  module subroutine GetVariableViscosity(this, t, u, mu, nu)
    class(INS_Operator_3D), intent(in)  :: this
    real(RNP),              intent(in)  :: t
    real(RNP),  contiguous, intent(in)  :: u(:,:,:,:,:)
    real(RNP),  contiguous, intent(out) :: mu(:,:,:,:)
    real(RNP),  contiguous, intent(out) :: nu(:,:,:,:)

    associate(x => this % sem_u % metrics %x)

      ! physical shear viscosity ...............................................

      if (this % problem % HasVariableProperties()) then
        call this % problem % GetViscosity(x, t, u, nu)
      else
        call SetArray(nu, this % nu_0)
      end if

      ! subgrid scale viscosity ................................................

      ! TBD

      ! artificial bulk viscosity ..............................................

      ! TBD: element dependent μ

      call SetArray(mu, this % mu_0)

    end associate

  end subroutine GetVariableViscosity

  !=============================================================================

end submodule MP_GetVariableViscosity
