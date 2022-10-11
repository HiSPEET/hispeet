!> summary:  3D gradient operator
!> author:   Joerg Stiller
!> date:     2022/10/09
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Grad__3D
  use Kind_Parameters, only: RDP
  use Standard_Operators__1D
  use Mesh__3D
  use Spectral_Element_Mesh__3D

  use TPO__Grad__3D_R ! TPO for regular meshes
  use TPO__Grad__3D_D ! TPO for deformed meshes

  private

  public :: TPO_Grad

  interface TPO_Grad
    module procedure TPO_Grad_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Automatic gradient TPO with double precision

  subroutine TPO_Grad_RDP(eop, sem, u, up, v)
    class(StandardOperators_1D),   intent(in) :: eop
    class(SpectralElementMesh_3D), intent(in) :: sem
    real(RDP),           intent(in)  :: u(:,:,:,:)   !< 3D scalar field u
    real(RDP), optional, intent(in)  :: up(:,:,:,:)  !< exterior traces u⁺
    real(RDP),           intent(out) :: v(:,:,:,:,:) !< gradient of u

    if (sem % mesh % regular) then
      call TPO_Grad( eop % w, eop % D, sem % mesh % dx, u, up, v )
    else
      call TPO_Grad( eop % w            &
                   , eop % D            &
                   , sem % metrics % Jd &
                   , sem % metrics % Ji &
                   , sem % metrics % a  &
                   , sem % metrics % n  &
                   , u, up, v           )
    end if

  end subroutine TPO_Grad_RDP

  !=============================================================================

end module TPO__Grad__3D
