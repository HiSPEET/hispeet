
  type Mesh3d_Variable
    class(Mesh3d_Metrics),  pointer :: metrics
    real(RNP), contiguous,  pointer :: value(:,:,:,:,:)      ! (0:po,0:po,0:po,n_elem,n_comp)
    real(RNP), allocatable, private :: storage(:,:,:,:,:)
  end type Mesh3d_Variable

! Remarks + issues:
! - original if storage is allocated
! - boundary conditions optional
! - extract slice (subset) as new variable
! – easy access to values, e.g: v(0:,0:,0:,1:,1:) => var % val(:,:,:,:,1:3)
! - extended types: scalar, vector
! - whether/how to include BC ?
    type(Mesh3d_BoundaryCondition), allocatable :: bcond(:)

  type Mesh3d_BoundaryVariable
    real(RNP), contiguous,  pointer :: value(:,:,:,:)      ! (0:po,0:po,n_bface,n_comp)
    real(RNP), allocatable, private :: storage(:,:,:,:)
  end type Mesh3d_BoundaryVariable

  type, extends(Mesh3d_BoundaryVariable) :: Mesh3d_BoundaryCondition
    character :: typ
  end type Mesh3d_BoundaryCondition

! Remarks
! - independent boundary variables required or always part of Mesh3d_Variable ?
! - visible or private ?
! - default: same number of components as mesh variable

