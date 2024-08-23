!> summary:  3D multilevel spectral element mesh
!> author:   Joerg Stiller
!> date:     2024/06/26
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ML__Mesh_Operators__3D
  use Execution_Control
  use HP__Refinement_Operator__1D
  use HP__Coarsening_Operator__1D
  use Spectral_Element_Mesh__3D
  use ML__Mesh__3D
  implicit none
  private

  public :: ML_MeshOperators_3D

  !-----------------------------------------------------------------------------
  !> 3D multilevel spectral element mesh

  type ML_MeshOperators_3D
    type(SpectralElementMesh_3D), allocatable :: sem(:)
      !< sequence of spectral element meshes [1:l_top]
    type(HP_RefinementOperator_1D), allocatable :: iop_cf(:)
      !< fine-to-coarse interpolation operators [1:l_top-1]
    type(HP_CoarseningOperator_1D), allocatable :: iop_fc(:)
      !< coarse-to-fine interpolation operators [2:l_top]
    type(HP_CoarseningOperator_1D), allocatable :: pop_fc(:)
      !< coarse-to-fine L2-projection operators [2:l_top]
  contains
    procedure :: Init_ML_MeshOperators_3D
  end type ML_MeshOperators_3D

 ! constructor interface
 interface ML_MeshOperators_3D
   procedure New_ML_MeshOperators_3D
 end interface

contains

  !-----------------------------------------------------------------------------
  !> Constructor of 3D multilevel mesh operators

  function New_ML_MeshOperators_3D(ml_mesh, po, basis, smooth) result(this)
    class(ML_Mesh_3D),   intent(in) :: ml_mesh !< multilevel mesh partition
    integer,             intent(in) :: po(:)   !< sequence of polynomial orders
    character, optional, intent(in) :: basis   !< 'G' or 'L' ['L']
    integer,   optional, intent(in) :: smooth  !< fine-to-coarse discontinuity
                                               !! smoothing {0,1,2} [0]
    type(ML_MeshOperators_3D) :: this

    call Init_ML_MeshOperators_3D(this, ml_mesh, po, basis, smooth)

  end function New_ML_MeshOperators_3D

  !-----------------------------------------------------------------------------
  !> Initialization of 3D multilevel mesh operators

  subroutine Init_ML_MeshOperators_3D(this, ml_mesh, po, basis, smooth)
    class(ML_MeshOperators_3D), intent(inout) :: this
    class(ML_Mesh_3D),   intent(in) :: ml_mesh !< multilevel mesh partition
    integer,             intent(in) :: po(:)   !< sequence of polynomial orders
    character, optional, intent(in) :: basis   !< 'G' or 'L' ['L']
    integer,   optional, intent(in) :: smooth  !< fine-to-coarse discontinuity
                                               !! smoothing {0,1,2} [0]

    type(HP_RefinementOptions_1D) :: cf_opt
    type(HP_CoarseningOptions_1D) :: fc_opt

    integer :: l, l_top

    associate(mesh => ml_mesh % mesh)

      ! initialization .........................................................

      if (size(po) < size(mesh)) then
        call Error('Init_ML_MeshOperators_3D', 'size(po) < l_top', &
                   'ML__Mesh_Operators__3D')
      end if

      l_top = size(mesh)

      allocate(this % sem    (1:l_top  ))
      allocate(this % iop_cf (1:l_top-1))
      allocate(this % iop_fc (2:l_top  ))
      allocate(this % pop_fc (2:l_top  ))

      if (present(basis)) then
        cf_opt % basis = basis
        fc_opt % basis = basis
      end if

      if (present(smooth)) then
        fc_opt % smooth = smooth
      end if

      ! creation of levels .....................................................

        do l = 1, l_top

        this % sem(l) = SpectralElementMesh_3D(mesh(l), po(l), basis)

        ! coarse-to-fine transfers operators
        if (l < l_top) then
          ! interpolation
          cf_opt % po_c    = po(l)
          cf_opt % po_f    = po(l+1)
          cf_opt % mode    = Mode(mesh(l)%refinement, po(l), po(l+1))
          this % iop_cf(l) = HP_RefinementOperator_1D(cf_opt)
        end if

        ! fine-to-coarse transfers operators
        if (l > 1) then
          ! interpolation
          fc_opt % po_f    = po(l)
          fc_opt % po_c    = po(l-1)
          fc_opt % method  = 'I'
          fc_opt % mode    = Mode(mesh(l-1)%refinement, po(l-1), po(l))
          this % iop_fc(l) = HP_CoarseningOperator_1D(fc_opt)
          ! L2-projection
          fc_opt % method  = 'P'
          this % pop_fc(l) = HP_CoarseningOperator_1D(fc_opt)
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

  !=============================================================================

end module ML__Mesh_Operators__3D
