!> summary:  Extended Fortran binding to MPI
!> author:   Joerg Stiller
!> date:     2014/11/06, revised 2017/04/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module XMPI

  use MPI_Binding

  use XMPI__Character
  use XMPI__Integer
  use XMPI__Integer_IXS
  use XMPI__Integer_IXL
  use XMPI__Logical
  use XMPI__Real_RSP
  use XMPI__Real_RDP

  logical, private :: initialized = .false.

contains

!-------------------------------------------------------------------------------
!> Initializes the extended MPI Fortran binding

subroutine XMPI_Init()

  if (initialized) then
    return
  else
    call Init_MPI_Binding
    initialized = .true.
  end if

end subroutine XMPI_Init

!===============================================================================

end module XMPI
