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

    associate(x => this % sem_u % metrics % x)

      ! physical shear viscosity ...............................................

      if (this % problem % HasVariableProperties()) then
        call this % problem % GetViscosity(x, t, u, nu)
      else
        call SetArray(nu, this % nu_0)
      end if

      ! subgrid scale viscosity ................................................

      ! TBD

      ! artificial bulk viscosity ..............................................

      if (this % c_mu > 0) then
        call GetVariableBulkViscosity(this, u, nu, mu)
      else
        call SetArray(mu, this % mu_0)
      end if

    end associate

  end subroutine GetVariableViscosity

  !-----------------------------------------------------------------------------
  !> Variable bulk viscosity

  subroutine GetVariableBulkViscosity(this, u, nu, mu)
    class(INS_Operator_3D), intent(in)  :: this
    real(RNP),  contiguous, intent(in)  :: u(:,:,:,:,:)
    real(RNP),  contiguous, intent(in)  :: nu(:,:,:,:)
    real(RNP),  contiguous, intent(out) :: mu(:,:,:,:)

    real(RNP) :: dx(3), hh, nn, vv
    integer :: e

    !$omp do
    do e = 1, size(mu,4)
      call this % mesh % element(e) % GetCuboidDimensions(dx)
      hh = (product(dx)**THIRD / this%eop_u%po)**2                     ! hh = h²
      vv = maxval(u(:,:,:,e,1)**2 + u(:,:,:,e,2)**2 + u(:,:,:,e,3)**2) ! vv = v²
      nn = maxval(nu(:,:,:,e))**2                                      ! nn = ν²
      mu(:,:,:,e) = this%c_mu * sqrt(nn + vv * hh)
    end do

  end subroutine GetVariableBulkViscosity

  !=============================================================================

end submodule MP_GetVariableViscosity
