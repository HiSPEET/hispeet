
  type MeshVariable_3D
    class(MeshMetrics_3D),  pointer :: metrics
    real(RNP), contiguous,  pointer :: value(:,:,:,:,:)      ! (0:po,0:po,0:po,n_elem,n_comp)
    real(RNP), allocatable, private :: storage(:,:,:,:,:)
  end type MeshVariable_3D

! Remarks + issues:
! - original if storage is allocated
! - boundary conditions optional
! - extract slice (subset) as new variable
! – easy access to values, e.g: v(0:,0:,0:,1:,1:) => var % val(:,:,:,:,1:3)
! - extended types: scalar, vector
! - whether/how to include BC ?
    type(MeshBoundaryCondition_3D), allocatable :: bcond(:)

  type MeshBoundaryVariable_3D
    real(RNP), contiguous,  pointer :: value(:,:,:,:)      ! (0:po,0:po,n_bface,n_comp)
    real(RNP), allocatable, private :: storage(:,:,:,:)
  end type MeshBoundaryVariable_3D

  type, extends(MeshBoundaryVariable_3D) :: MeshBoundaryCondition_3D
    character :: typ
  end type MeshBoundaryCondition_3D

! Remarks
! - independent boundary variables required or always part of MeshVariable_3D ?
! - visible or private ?
! - default: same number of components as mesh variable

