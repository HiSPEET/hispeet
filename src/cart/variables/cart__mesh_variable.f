!> summary:  Multi-component variable defined on mesh partition
!> author:   Joerg Stiller
!> date:     2016/03/18
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Multi-component variable defined on mesh partition
!>
!> @todo
!> This is just a sketch
!> @endtodo
!===============================================================================

module Mesh_Variable

  use Kind_Parameters, only: RNP
  use CART__Mesh_Partition

  implicit none
  private

  public :: MeshVariable

  !-----------------------------------------------------------------------------
  !> Type for handling mesh variables

  type MeshVariable

    type(MeshPartition),   pointer :: mesh => null() !< associated partition
    real(RNP), contiguous, pointer :: val(:,:,:,:,:) => null() !< values
    character(:),          pointer :: name(:) => null()! component names

    logical :: is_original = .false. !< indicates original mesh variable
    logical :: has_ghosts  = .false. !< indicates whether ghosts are included

  contains

  ! New
  ! Set
  ! Merge
  ! GetTrace:  extract trace variable
  ! GetHandle: return mesh variable linked to a section of the given one
  ! GetClone?  get full or partial copy
  ! Access:    get direct access to single or section of components

  end type MeshVariable

!===============================================================================

end module Mesh_Variable
