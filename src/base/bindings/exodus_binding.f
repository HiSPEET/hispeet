!> summary:  Interface to subset of ExodusII routines
!> author:   Joerg Stiller
!> date:     2013/05/13
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Interface to subset of ExodusII routines
!===============================================================================

module Exodus_Binding
  use C_Binding
  implicit none
  private

  public :: ExOpen
  public :: ExClose
  public :: ExGetInit
  public :: ExGetCoord
  public :: ExGetElemBlock
  public :: ExGetElemBlkIds
  public :: ExGetElemConn
  public :: ExGetName
  public :: ExGetSideSetIds
  public :: ExGetSideSetParam
  public :: ExGetSideSet

  ! Exodus control parameters and dimensions ...................................
  ! NOTE: values adopted from exodusII.h

  ! Exodus control parameters
  integer(C_INT), public ::     &
    EX_WRITE = 1,               & ! open existing file for appending
    EX_READ  = 2                  ! open file for reading (default)

  ! Exodus string leghths (without null char)
  integer, parameter, public :: &
    EX_MAX_STR_LENGTH  =  32,   & ! = MAX_STR_LENGTH
    EX_MAX_NAME_LENGTH =  32,   & ! = MAX_NAME_LENGTH
    EX_MAX_LINE_LENGTH =  80,   & ! = MAX_LINE_LENGTH
    EX_MAX_ERR_LENGTH  = 256      ! = MAX_ERR_LENGTH

  ! Exodus entity types
  integer(C_INT), parameter, public :: &
    EX_ELEM_BLOCK  =  1,        &
    EX_SIDE_SET    =  3

  ! Interfaces to ExodusII C routines ..........................................

  interface

    function ExOpen(path, mode, comp_ws, io_ws, version) result(exoid) &
        bind(C, name='exf_open')
      import :: C_CHAR, C_INT, C_FLOAT
      character(C_CHAR),     intent(in)    :: path(*)
      integer(C_INT), value, intent(in)    :: mode
      integer(C_INT),        intent(inout) :: comp_ws
      integer(C_INT),        intent(inout) :: io_ws
      real(C_FLOAT),         intent(out)   :: version
      integer(C_INT)                       :: exoid
    end function ExOpen

    subroutine ExClose(exoid, err) bind(C, name='exf_close')
      import :: C_INT
      integer(C_INT), value, intent(in)  :: exoid
      integer(C_INT),        intent(out) :: err
    end subroutine ExClose

    subroutine ExGetInit(exoid, title, num_dim, num_nodes, num_elem, &
        num_elem_blk, num_node_sets, num_side_sets, err) &
        bind(C, name='exf_get_init')
      import :: C_INT, C_CHAR
      integer(C_INT), value, intent(in)  :: exoid
      character(C_CHAR),     intent(out) :: title(*)
      integer(C_INT),        intent(out) :: num_dim
      integer(C_INT),        intent(out) :: num_nodes
      integer(C_INT),        intent(out) :: num_elem
      integer(C_INT),        intent(out) :: num_elem_blk
      integer(C_INT),        intent(out) :: num_node_sets
      integer(C_INT),        intent(out) :: num_side_sets
      integer(C_INT),        intent(out) :: err
    end subroutine ExGetInit

    subroutine ExGetElemBlock(exoid, elem_blk_id, elem_type,  &
        num_elem_this_blk, num_nodes_per_elem, num_attr, err) &
        bind(C, name='exf_get_elem_block')
      import :: C_INT, C_CHAR
      integer(C_INT), value, intent(in)  :: exoid
      integer(C_INT), value, intent(in)  :: elem_blk_id
      character(C_CHAR),     intent(out) :: elem_type(*)
      integer(C_INT),        intent(out) :: num_elem_this_blk
      integer(C_INT),        intent(out) :: num_nodes_per_elem
      integer(C_INT),        intent(out) :: num_attr
      integer(C_INT),        intent(out) :: err
    end subroutine ExGetElemBlock

    subroutine ExGetElemBlkIds(exoid, elem_blk_ids, err) &
        bind(C, name='exf_get_elem_blk_ids')
      import :: C_INT
      integer(C_INT), value, intent(in)  :: exoid
      integer(C_INT),        intent(out) :: elem_blk_ids(*)
      integer(C_INT),        intent(out) :: err
    end subroutine ExGetElemBlkIds

    subroutine ExGetElemConn(exoid, elem_blk_id, connect, err) &
        bind(C, name='exf_get_elem_conn')
      import :: C_INT
      integer(C_INT), value, intent(in)    :: exoid
      integer(C_INT), value, intent(in)    :: elem_blk_id
      integer(C_INT),        intent(inout) :: connect(*)
      integer(C_INT),        intent(out)   :: err
    end subroutine ExGetElemConn

    subroutine ExGetName(exoid, obj_type, id, name, err) &
        bind(C, name='exf_get_name')
      import :: C_INT, C_CHAR
      integer(C_INT), value, intent(in)  :: exoid
      integer(C_INT), value, intent(in)  :: obj_type
      integer(C_INT), value, intent(in)  :: id
      character(C_CHAR),     intent(out) :: name(*)
      integer(C_INT),        intent(out) :: err
    end subroutine ExGetName

    subroutine ExGetSideSetIds(exoid, side_set_ids, err) &
        bind(C, name='exf_get_side_set_ids')
      import :: C_INT
      integer(C_INT), value, intent(in)  :: exoid
      integer(C_INT),        intent(out) :: side_set_ids(*)
      integer(C_INT),        intent(out) :: err
    end subroutine ExGetSideSetIds

    subroutine ExGetSideSetParam(exoid, side_set_id, num_side_in_set, &
        num_dist_fact_in_set, err) &
        bind(C, name='exf_get_side_set_param')
      import :: C_INT
      integer(C_INT), value, intent(in)  :: exoid
      integer(C_INT), value, intent(in)  :: side_set_id
      integer(C_INT),        intent(in)  :: num_side_in_set
      integer(C_INT),        intent(in)  :: num_dist_fact_in_set
      integer(C_INT),        intent(out) :: err
    end subroutine ExGetSideSetParam

    subroutine ExGetSideSet(exoid, side_set_id, side_set_elem_list, &
        side_set_side_list, err) &
        bind(C, name='exf_get_side_set')
      import :: C_INT
      integer(C_INT), value, intent(in)  :: exoid
      integer(C_INT), value, intent(in)  :: side_set_id
      integer(C_INT),        intent(out) :: side_set_elem_list(*)
      integer(C_INT),        intent(out) :: side_set_side_list(*)
      integer(C_INT),        intent(out) :: err
    end subroutine ExGetSideSet

  end interface

  ! generic interfaces to float and double versions ............................

  interface ExGetCoord
    module procedure ExGetCoord_Float
    module procedure ExGetCoord_Double
  end interface ExGetCoord

contains

!===============================================================================
! specific implementations of ExGetCoord

!-------------------------------------------------------------------------------
!> ExGetCoord -- float version

subroutine ExGetCoord_Float(exoid, x_coor, y_coor, z_coor, err)
  integer(C_INT), value, intent(in)  :: exoid
  real(C_FLOAT),         intent(out) :: x_coor(*)
  real(C_FLOAT),         intent(out) :: y_coor(*)
  real(C_FLOAT),         intent(out) :: z_coor(*)
  integer(C_INT),        intent(out) :: err

  interface
    subroutine exf_get_coord(exoid, x_coor, y_coor, z_coor, err) &
        bind(C, name='exf_get_coord')
      import :: C_INT, C_FLOAT
      integer(C_INT), value, intent(in)  :: exoid
      real(C_FLOAT),         intent(out) :: x_coor(*)
      real(C_FLOAT),         intent(out) :: y_coor(*)
      real(C_FLOAT),         intent(out) :: z_coor(*)
      integer(C_INT),        intent(out) :: err
    end subroutine exf_get_coord
  end interface

  call exf_get_coord(exoid, x_coor, y_coor, z_coor, err)

end subroutine ExGetCoord_Float

!-------------------------------------------------------------------------------
!> ExGetCoord -- double version

subroutine ExGetCoord_Double(exoid, x_coor, y_coor, z_coor, err)
  integer(C_INT), value, intent(in)  :: exoid
  real(C_DOUBLE),        intent(out) :: x_coor(*)
  real(C_DOUBLE),        intent(out) :: y_coor(*)
  real(C_DOUBLE),        intent(out) :: z_coor(*)
  integer(C_INT),        intent(out) :: err

  interface
    subroutine exf_get_coord(exoid, x_coor, y_coor, z_coor, err) &
        bind(C, name='exf_get_coord')
      import :: C_INT, C_DOUBLE
      integer(C_INT), value, intent(in)  :: exoid
      real(C_DOUBLE),        intent(out) :: x_coor(*)
      real(C_DOUBLE),        intent(out) :: y_coor(*)
      real(C_DOUBLE),        intent(out) :: z_coor(*)
      integer(C_INT),        intent(out) :: err
    end subroutine exf_get_coord
  end interface

  call exf_get_coord(exoid, x_coor, y_coor, z_coor, err)

end subroutine ExGetCoord_Double

!===============================================================================

end module Exodus_Binding
