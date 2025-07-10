!> summary:  3D multilevel spectral element space operators
!> author:   Joerg Stiller
!> date:     2024/06/26
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ML__Mesh_Operators__3D
  use Kind_Parameters, only: RNP
  use Execution_Control
  use HP__Refinement_Operator__1D
  use HP__Coarsening_Operator__1D
  use Spectral_Element_Mesh__3D
  use ML__Mesh__3D
  implicit none
  private

  public :: ML_MeshOperators_3D

  !-----------------------------------------------------------------------------
  !> 3D multilevel spectral element mesh operators

  type ML_MeshOperators_3D
    type(SpectralElementMesh_3D), allocatable :: sem(:)
      !< sequence of spectral element meshes [1:l_top]
    type(HP_RefinementOperator_1D), allocatable :: iop_cf_x(:)
      !< coarse-to-fine interpolation operators in space [1:l_top-1]
    type(HP_CoarseningOperator_1D), allocatable :: iop_fc_x(:)
      !< fine-to-coarse interpolation operators in space [2:l_top]
    type(HP_CoarseningOperator_1D), allocatable :: pop_fc_x(:)
      !< fine-to-coarse L2-projection operators in space [2:l_top]
  contains
    procedure :: Init_ML_MeshOperators_3D
    procedure :: Get_Volume
    procedure :: Get_SurfaceAreas
  end type ML_MeshOperators_3D

  ! constructor interface
  interface ML_MeshOperators_3D
    procedure New_ML_MeshOperators_3D
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Constructor of 3D multilevel mesh operators

  function New_ML_MeshOperators_3D(ml_mesh, po, nodes, smooth) result(this)
    class(ML_Mesh_3D),   intent(in) :: ml_mesh !< multilevel mesh partition
    integer,             intent(in) :: po(:)   !< sequence of polynomial orders
    character, optional, intent(in) :: nodes   !< 'G' or 'L' ['L']
    integer,   optional, intent(in) :: smooth  !< fine-to-coarse discontinuity
                                               !! smoothing {0,1,2} [0]
    type(ML_MeshOperators_3D) :: this

    call Init_ML_MeshOperators_3D(this, ml_mesh, po, nodes, smooth)

  end function New_ML_MeshOperators_3D

  !-----------------------------------------------------------------------------
  !> Initialization of 3D multilevel mesh operators

  subroutine Init_ML_MeshOperators_3D(this, ml_mesh, po, nodes, smooth)
    class(ML_MeshOperators_3D), intent(inout) :: this
    class(ML_Mesh_3D),   intent(in) :: ml_mesh !< multilevel mesh partition
    integer,             intent(in) :: po(:)   !< sequence of polynomial orders
    character, optional, intent(in) :: nodes   !< 'G' or 'L' ['L']
    integer,   optional, intent(in) :: smooth  !< fine-to-coarse discontinuity
                                               !! smoothing {0,1,2} [0]

    type(HP_RefinementOptions_1D) :: opt_cf
    type(HP_CoarseningOptions_1D) :: opt_fc

    integer :: l, l_top

    associate(mesh => ml_mesh % mesh)

      ! initialization .........................................................

      if (size(po) < size(mesh)) then
        call Error('Init_ML_MeshOperators_3D', 'size(po) < l_top', &
                   'ML__Mesh_Operators__3D')
      end if

      l_top = size(mesh)

      allocate(this % sem      (1:l_top  ))
      allocate(this % iop_cf_x (1:l_top-1))
      allocate(this % iop_fc_x (2:l_top  ))
      allocate(this % pop_fc_x (2:l_top  ))

      if (present(nodes)) then
        opt_cf % nodes = nodes
        opt_fc % nodes = nodes
      end if

      if (present(smooth)) then
        opt_fc % smooth = smooth
      end if

      ! creation of levels .....................................................

      do l = 1, l_top

        this % sem(l) = SpectralElementMesh_3D(mesh(l), po(l), nodes)

        ! coarse-to-fine transfers operators
        if (l < l_top) then
          ! interpolation
          opt_cf % po_c = po(l)
          opt_cf % po_f = po(l+1)
          opt_cf % mode = Mode(mesh(l)%refinement, po(l), po(l+1))
          this % iop_cf_x(l) = HP_RefinementOperator_1D(opt_cf)
        end if

        ! fine-to-coarse transfers operators
        if (l > 1) then
          ! interpolation
          opt_fc % po_f   = po(l)
          opt_fc % po_c   = po(l-1)
          opt_fc % method = 'I'
          opt_fc % mode   = Mode(mesh(l-1)%refinement, po(l-1), po(l))
          this % iop_fc_x(l) = HP_CoarseningOperator_1D(opt_fc)
          ! L2-projection
          opt_fc % method = 'P'
          this % pop_fc_x(l) = HP_CoarseningOperator_1D(opt_fc)
        end if

        ! skip levels above toplevel mesh
        if (mesh(l)%is_top .and. l < l_top) then
          call Warning('Init_ML_MeshOperators_3D', 'ML mesh id incomplete', &
                       'ML__Mesh_Operators__3D')
          exit
        end if

      end do

    end associate

  end subroutine Init_ML_MeshOperators_3D

  !-----------------------------------------------------------------------------
  !> Refinement mode from given refinement type and polynomial orders

  pure integer function Mode(refinement, po_c, po_f)
    character, intent(in) :: refinement !< 'c' cloning or 's' subdivision
    integer,   intent(in) :: po_c       !< polynomial order of coarse level
    integer,   intent(in) :: po_f       !< polynomial order of fine level

    select case(refinement)
    case('c')
      ! cloning
      if (po_c /= po_f) then
        ! p-refinement
        Mode = 1
      else
        ! identity
        Mode = 0
      end if
    case('s')
      ! subdivision
      Mode = 2
    case default
      ! no or unknown refinement type
      Mode = -1
    end select

  end function Mode

  !-----------------------------------------------------------------------------
  !> TBP for computing the volume of the computational domain

  subroutine Get_Volume(this, vol)
    class(ML_MeshOperators_3D), intent(in) :: this
    real(RNP), intent(out) :: vol !< volume, should be PRIVATE with OpenMP

    real(RNP) :: vol_l
    integer :: l

    vol = 0
    do l = 1, size(this % sem)
      call this % sem(l) % Get_Volume(vol_l, leaf = .true.)
      vol = vol + vol_l
    end do

  end subroutine Get_Volume

  !-----------------------------------------------------------------------------
  !> TBP for computing the areas of the boundary surfaces

  subroutine Get_SurfaceAreas(this, area)
    class(ML_MeshOperators_3D), intent(in) :: this
    real(RNP), intent(out) :: area(:) !< areas, should be PRIVATE with OpenMP

    real(RNP) :: area_l(size(area))
    integer :: l

    area_l = 0
    do l = 1, size(this % sem)
      call this % sem(l) % Get_SurfaceAreas(area_l, leaf = .true.)
      area = area + area_l
    end do

  end subroutine Get_SurfaceAreas

  !=============================================================================

end module ML__Mesh_Operators__3D
