!> summary:  Definition of intrinsic type kind parameters
!> author:   Joerg Stiller
!> date:     2014/03/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Definition of kind parameters for intrinsic types
!===============================================================================

module Kind_Parameters
  use, intrinsic :: ISO_Fortran_Env
  implicit none
  public

  !-----------------------------------------------------------------------------
  ! integer kinds

  integer, parameter :: IXS = selected_int_kind( 2) !<  8 Bit, range > 10^2
  integer, parameter :: IXL = selected_int_kind(18) !< 64 Bit, range > 10^18

  !-----------------------------------------------------------------------------
  ! real kinds

  integer, parameter :: RSP = kind(1E0) !< real single precision
  integer, parameter :: RDP = kind(1D0) !< real double precision

  !> real high precision: 128 bit, if available
  integer, parameter :: RHP = merge(REAL128, RDP, REAL128 > 0)

  !> real normal precision
  integer, parameter :: RNP = RDP

  !=============================================================================

end module Kind_Parameters