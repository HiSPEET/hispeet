!> summary:  DG Schwarz operator for elliptic equations: base type
!> author:   Joerg Stiller
!> date:     2022/02/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module DG__Schwarz_Operator__3D__Base
! use Kind_Parameters

  implicit none
  private

  public :: DG_SCHWARZ_BC_3D
  public :: DG_SchwarzOperator_3D

  !-----------------------------------------------------------------------------
  !> Supported boundary conditions
  !>
  !> Periodic boundaries are treated as interior: 'P' → ' '

  character, parameter :: DG_SCHWARZ_BC_3D(3) = [ ' ', 'D', 'N']

  !-----------------------------------------------------------------------------
  !> Schwarz operator

  type, abstract :: DG_SchwarzOperator_3D

    integer :: po = -1     !< polynomial order
    integer :: no = -1     !< number of overlapped layers
    integer :: nc = -1     !< number of 1D configurations
    logical :: restrictive !< T/F for ex/including neighbor results

    integer, allocatable :: cfg(:,:) !< subdomain configurations (nc,*)

  contains

    procedure, nopass :: ConfigurationID

  end type DG_SchwarzOperator_3D

contains

  !=============================================================================
  ! Utilities

  !-----------------------------------------------------------------------------
  !> Returns the 1D subdomain configuration ID corresponding to the given BCs

  pure integer function ConfigurationID(bc) result(cfg)
    character, intent(in) :: bc(2) !< left/right boundary types {' ','D','N','P'}

    integer :: i1, i2

    select case(bc(1))
    case('D')
      i1 = 2
    case('N')
      i1 = 3
    case default ! ' ' and 'P'
      i1 = 1
    end select

    select case(bc(2))
    case('D')
      i2 = 2
    case('N')
      i2 = 3
    case default ! ' ' and 'P'
      i2 = 1
    end select

    cfg = i1 + (i2 - 1) * size(DG_SCHWARZ_BC_3D)

  end function ConfigurationID

  !=============================================================================

end module DG__Schwarz_Operator__3D__Base
