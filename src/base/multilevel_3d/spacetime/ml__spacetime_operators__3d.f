!> summary:  3D multilevel spectral element spacetime operators
!> author:   Joerg Stiller
!> date:     2024/11/24
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ML__Spacetime_Operators__3D
  use Spectral_Deferred_Correction
  use HP__Refinement_Operator__1D
  use HP__Coarsening_Operator__1D
  use ML__Mesh_Operators__3D
  implicit none
  private

  public :: ML_SpacetimeOperators_3D

  !-----------------------------------------------------------------------------
  !> 3D multilevel spectral element spacetime mesh operators

  type, extends(ML_MeshOperators_3D) :: ML_SpacetimeOperators_3D

    integer, allocatable :: n_time(:)
      !< number of time steps in one slice [1:l_top]
    type(SDC_Method), allocatable :: sdc(:)
      !< temporal SDC/collocation method [1:l_top]
    type(HP_RefinementOperator_1D), allocatable :: iop_cf_t(:)
      !< fine-to-coarse interpolation operators in time [1:l_top-1]
    type(HP_CoarseningOperator_1D), allocatable :: iop_fc_t(:)
      !< coarse-to-fine interpolation operators in time [2:l_top]
    type(HP_CoarseningOperator_1D), allocatable :: pop_fc_t(:)
      !< coarse-to-fine L2-projection operators in time [2:l_top]

! contains
!
!   procedure :: Init_ML_SpacetimeOperators_3D

  end type ML_SpacetimeOperators_3D

  ! constructor interface
! interface ML_SpacetimeOperators_3D
!   procedure New_ML_SpacetimeOperators_3D
! end interface

contains

  !-----------------------------------------------------------------------------
  !> Constructor of 3D multilevel spacetime operators



  !-----------------------------------------------------------------------------
  !> Initialization of 3D multilevel spacetime operators


  !=============================================================================

end module ML__Spacetime_Operators__3D
