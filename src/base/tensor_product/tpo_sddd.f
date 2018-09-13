!> summary:  Evaluation of  v = s(DxDxD) u with scalar s and diagonal D
!> author:   Joerg Stiller
!> date:     2017/05/07
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Evaluation of  v = s(DxDxD) u with scalar s and diagonal D
!>
!> @note
!>   *  As the benefit of parametrization is uncertain, only a generic variant
!>      is provided, so far.
!>   *  For larger operators (np > 10?) it may be more efficient to collapse
!>      dimensions to apply the resulting 3D diagonal operator
!> @endnote
!===============================================================================

module TPO_sDDD
  use Kind_Parameters, only: RNP
  implicit none
  private

  public :: TPO_sDDD_Proc    ! procedure interface
  public :: TPO_sDDD_Assign  ! procedure assignment
  public :: TPO_sDDD_Eval    ! operator evaluation


  abstract interface
    subroutine TPO_sDDD_Proc(np, ne, s, D, u, v)
      use Kind_Parameters, only: RNP
      integer,   intent(in)  :: np             !< number of points per direction
      integer,   intent(in)  :: ne             !< number of elements
      real(RNP), intent(in)  :: s              !< scaling factor
      real(RNP), intent(in)  :: D(np)          !< diagonal 1D operator
      real(RNP), intent(in)  :: u(np,np,np,ne) !< 3D operand
      real(RNP), intent(out) :: v(np,np,np,ne) !< 3D result
    end subroutine TPO_sDDD_Proc
  end interface

contains

!-------------------------------------------------------------------------------
!> Assingment of a DDD operator procedure matching given dimensions

subroutine TPO_sDDD_Assign(np, Proc)
  integer, intent(in) :: np  !< operator dimension
  procedure(TPO_sDDD_Proc), pointer, intent(out) :: Proc

  ! external procedures ........................................................

  procedure(TPO_sDDD_Proc) :: TPO_sDDD__gen

  ! initialization .............................................................

  Proc => null()

  ! parametrized procedures ....................................................

  ! *** not implemented ***

  ! generic procedures .........................................................

  if (.not. associated(Proc)) then
    Proc => TPO_sDDD__gen
  end if

end subroutine TPO_sDDD_Assign

!-------------------------------------------------------------------------------
!> Evaluation of the sDDD operator

subroutine TPO_sDDD_Eval(np, ne, s, D, u, v)
  integer,   intent(in)  :: np             !< number of points per direction
  integer,   intent(in)  :: ne             !< number of elements
  real(RNP), intent(in)  :: s              !< scaling factor
  real(RNP), intent(in)  :: D(np)          !< diagonal 1D operator
  real(RNP), intent(in)  :: u(np,np,np,ne) !< 3D operand
  real(RNP), intent(out) :: v(np,np,np,ne) !< 3D result

  procedure(TPO_sDDD_Proc), pointer :: Proc

  call TPO_sDDD_Assign(np, Proc)
  call Proc(np, ne, s, D, u, v)

end subroutine TPO_sDDD_Eval

!===============================================================================

end module TPO_sDDD
