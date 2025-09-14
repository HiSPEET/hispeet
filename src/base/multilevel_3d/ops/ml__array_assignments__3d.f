!> summary:  Multilevel array assignments
!> author:   Joerg Stiller
!> date:     2025/09/13
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ML__Array_Assignments__3D
  use Kind_Parameters
  use Execution_Control
  use Array_Assignments
  use XMPI
  use ML__Mesh_Variable__3D
  implicit none
  private

  public :: ML_SetArray_3D
  public :: ML_MergeArrays_3D
  public :: ML_CalibrateArray_3D

  interface ML_SetArray_3D
    module procedure ML_SetArray_S
    module procedure ML_SetArray_A
  end interface

contains

  !=============================================================================
  ! SetArray

  !-----------------------------------------------------------------------------
  !> Assigment of a scalar

  subroutine ML_SetArray_S(a, s, l_top)
    class(ML_MeshVariable_3D), intent(inout) :: a
    real(RNP),         intent(in) :: s     !< assigned scalar
    integer, optional, intent(in) :: l_top !< top level  [auto]

    integer :: l, l_top_

    if (present(l_top)) then
      l_top_ = min(l_top, size(a % level))
    else
      l_top_ = size(a % level)
    end if

    do l = 1, l_top_
      call SetArray(a % level(l) % val, s, multi = .true.)
    end do

  end subroutine ML_SetArray_S

  !-----------------------------------------------------------------------------
  !> Assigment of a matching array variable

  subroutine ML_SetArray_A(a, b, l_top)
    class(ML_MeshVariable_3D), intent(inout) :: a
    class(ML_MeshVariable_3D), intent(in)    :: b
    integer,         optional, intent(in)    :: l_top !< top level  [auto]

    integer :: l, l_top_

    if (present(l_top)) then
      l_top_ = min(l_top, size(a % level))
    else
      l_top_ = size(a % level)
    end if

    do l = 1, l_top_
      call SetArray(a % level(l) % val, b % level(l) % val, multi = .true.)
    end do

  end subroutine ML_SetArray_A

  !=============================================================================
  ! MergeArrays

  !-----------------------------------------------------------------------------
  !> Performs  a = alpha * a + beta * b + s  with arrays a, b and scalar s

  subroutine ML_MergeArrays_3D(alpha, a, beta, b, s, l_top)
    real(RNP),                 intent(in)    :: alpha !< coefficient to a
    class(ML_MeshVariable_3D), intent(inout) :: a     !< array merged to
    real(RNP),                 intent(in)    :: beta  !< coefficient to b
    class(ML_MeshVariable_3D), intent(in)    :: b     !< array to merge
    real(RNP),       optional, intent(in)    :: s     !< scalar [0]
    integer,         optional, intent(in)    :: l_top !< top level  [auto]

    integer :: l, l_top_

    if (present(l_top)) then
      l_top_ = min(l_top, size(a % level))
    else
      l_top_ = size(a % level)
    end if

    do l = 1, l_top_
      call MergeArrays( alpha, a % level(l) % val                    &
                      , beta , b % level(l) % val, s, multi = .true. )
    end do

  end subroutine ML_MergeArrays_3D

  !=============================================================================
  ! CalibrateArray

  !-----------------------------------------------------------------------------
  !> Remove non-zero arithmetic mean over all levels, ignoring frozen elements

  subroutine ML_CalibrateArray_3D(a)
    class(ML_MeshVariable_3D), intent(inout) :: a

    real(RNP), save :: s_loc = 0, s_glob, a_mean
    integer,   save :: n_loc = 0, n_glob

    integer :: e, k, l, l_top

    l_top = size(a % level)

    ! local amounts ............................................................

    do l = 1, l_top
      associate( nc  => a % level(l) % nc                    &
               , ne  => a % level(l) % mesh % n_elem_active  &
               , a_l => a % level(l) % val                   )

        if (ne == 0) cycle

        !$omp master
        n_loc = n_loc + size(a_l(:,:,:,1:ne,1:nc))
        !$omp end master

        !$omp do collapse(2), reduction(+:s_loc) private(e,k)
        do k = 1, nc
        do e = 1, ne
          s_loc = s_loc + sum(a_l(:,:,:,e,k))
        end do
        end do

      end associate
    end do

    ! mean value ...............................................................

    !$omp master
    associate(comm => a % level(1) % mesh % comm_world)
      call XMPI_Allreduce(n_loc, n_glob, MPI_SUM, comm)
      call XMPI_Allreduce(s_loc, s_glob, MPI_SUM, comm)
      a_mean = s_glob / n_glob
      n_loc = 0
      s_loc = 0
    end associate
    !$omp end master
    !$omp barrier

    ! calibration to zero mean value ...........................................

    do l = 1, l_top
      associate( nc  => a % level(l) % nc                    &
               , ne  => a % level(l) % mesh % n_elem_active  &
               , a_l => a % level(l) % val                   )

        !$omp do collapse(2)
        do k = 1, nc
        do e = 1, ne
          a_l(:,:,:,e,k) = a_l(:,:,:,e,k) - a_mean
        end do
        end do

      end associate
    end do

  end subroutine ML_CalibrateArray_3D

  !=============================================================================

end module ML__Array_Assignments__3D
