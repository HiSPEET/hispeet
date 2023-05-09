!> summary:  Generate a generic 3d mesh from Gmsh meshfile
!> author:   Matthias Frey, Joerg Stiller
!> date:     2023/05/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @note
!> This version is restricted to hexahedral elements.
!===============================================================================

module Import__GMSH__3D
  use Kind_Parameters
  use Execution_Control   ! to throw warnings or errors
  use Generic_Mesh__3D
  implicit none
  private

  public :: Import_GMSH_3D

contains

  !-----------------------------------------------------------------------------
  !> Reads the mesh from file and converts it to a GenericMesh_3D object

  subroutine Import_GMSH_3D(file, mesh)
    character(len=*),      intent(in)  :: file !< Exodus file
    class(GenericMesh_3D), intent(out) :: mesh !< output mesh

    !---------------------------------------------------------------------------
    ! Declarations

    ! GMSH interface ...........................................................

    ! declare here variables for reading data from the Gmsh file
    ! - adopt GMSH names for clarity
    ! - think about defining a data type for accommodating the Gmsh data

    ! auxiliary variables ......................................................

    !---------------------------------------------------------------------------
    ! Read Gmsh mesh file

    ! dimensions ...............................................................

    ! or other meta data ...

    ! nodes ....................................................................

    ! elements .................................................................

    ! ...

    !---------------------------------------------------------------------------
    !> Process Gmsh mesh data

    ! perform necessary transformations ...
    !   - identification/numbering of vertices
    !   - identification of element points
    !   - transformation from equidistant to Lobatto points
    !   - identification of boundaries

    !---------------------------------------------------------------------------
    !> Create HiSPEET generic mesh

    ! ...

  end subroutine Import_GMSH_3D

  !=============================================================================

end module Import__GMSH__3D
