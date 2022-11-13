!> summary:  3D divergence operator
!> author:   Joerg Stiller
!> date:     2022/10/09
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__Div__3D
  use Kind_Parameters, only: RDP
  use Standard_Operators__1D
  use Mesh__3D
  use Spectral_Element_Mesh__3D

  use TPO__Div__3D_R ! TPO for regular meshes
  use TPO__Div__3D_D ! TPO for deformed meshes

  private

  public :: TPO_Div

  interface TPO_Div
    module procedure TPO_Div_RDP
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Automatic divergence TPO with double precision

  subroutine TPO_Div_RDP(eop, sem, u, up, v)
    class(StandardOperators_1D),   intent(in) :: eop
    class(SpectralElementMesh_3D), intent(in) :: sem
    real(RDP),           intent(in)  :: u(:,:,:,:,:)  !< 3D vector field
    real(RDP), optional, intent(in)  :: up(:,:,:,:,:) !< exterior traces u⁺
    real(RDP),           intent(out) :: v(:,:,:,:)    !< divergence of u

    if (sem % mesh % regular) then
      call TPO_Div_R( eop % w, eop % D, sem % mesh % dx, u, up, v )
    else
      call TPO_Div_D( eop % w            &
                    , eop % D            &
                    , sem % metrics % Jd &
                    , sem % metrics % Ji &
                    , sem % metrics % a  &
                    , sem % metrics % n  &
                    , u, up, v           )
    end if

  end subroutine TPO_Div_RDP

  !=============================================================================

end module TPO__Div__3D
