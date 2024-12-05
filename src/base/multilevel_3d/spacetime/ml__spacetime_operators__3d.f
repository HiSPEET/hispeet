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
  use ML__Mesh__3D
  implicit none
  private

  public :: ML_SpacetimeOperators_3D

  !-----------------------------------------------------------------------------
  !> 3D multilevel spectral element spacetime mesh operators

  type, extends(ML_MeshOperators_3D) :: ML_SpacetimeOperators_3D

    integer, allocatable :: po_t(:)
      !< polynomial degree in time [1:l_top]
    integer, allocatable :: ne_t(:)
      !< number of time elements in one slice [1:l_top]
    type(SDC_Method), allocatable :: sdc(:)
      !< temporal SDC/collocation method [1:l_top]
    type(HP_RefinementOperator_1D), allocatable :: iop_cf_t(:)
      !< coarse-to-fine interpolation operators in time [1:l_top-1]
    type(HP_CoarseningOperator_1D), allocatable :: iop_fc_t(:)
      !< fine-to-coarse interpolation operators in time [2:l_top]
    type(HP_CoarseningOperator_1D), allocatable :: pop_fc_t(:)
      !< fine-to-coarse L2-projection operators in time [2:l_top]

  contains

    procedure :: Init_ML_SpacetimeOperators_3D

  end type ML_SpacetimeOperators_3D

  ! constructor interface
  interface ML_SpacetimeOperators_3D
    procedure New_ML_SpacetimeOperators_3D
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Constructor of 3D multilevel spacetime operators

  function New_ML_SpacetimeOperators_3D &
      (ml_mesh, po_x, po_t, ne_t, nodes_x, nodes_t, smooth) result(this)

    class(ML_Mesh_3D),      intent(in) :: ml_mesh !< multilevel mesh partition
    integer, contiguous,    intent(in) :: po_x(:) !< polynomial orders in x
    integer, contiguous,    intent(in) :: po_t(:) !< polynomial orders in t
    integer, contiguous,    intent(in) :: ne_t(:) !< elements per slice in t
    character(2), optional, intent(in) :: nodes_x !< x nodes {G,L}     ['L']
    character(2), optional, intent(in) :: nodes_t !< t nodes {E,L,RR}  ['RR']
    integer     , optional, intent(in) :: smooth  !< fine-to-coarse discontinuity
                                                  !! smoothing {0,1,2} [0]
    type(ML_SpacetimeOperators_3D) :: this

    call Init_ML_SpacetimeOperators_3D( this, ml_mesh, po_x, po_t, ne_t &
                                      , nodes_x, nodes_t, smooth  )

  end function New_ML_SpacetimeOperators_3D

  !-----------------------------------------------------------------------------
  !> Initialization of 3D multilevel spacetime operators

  subroutine Init_ML_SpacetimeOperators_3D &
      (this, ml_mesh, po_x, po_t, ne_t, nodes_x, nodes_t, smooth)

    class(ML_SpacetimeOperators_3D), intent(inout) :: this
    class(ML_Mesh_3D),      intent(in) :: ml_mesh !< multilevel mesh partition
    integer, contiguous,    intent(in) :: po_x(:) !< polynomial orders in x
    integer, contiguous,    intent(in) :: po_t(:) !< polynomial orders in t
    integer, contiguous,    intent(in) :: ne_t(:) !< number of elements in t
    character(2), optional, intent(in) :: nodes_x !< x nodes {G,L}     ['L']
    character(2), optional, intent(in) :: nodes_t !< t nodes {E,L,RR}  ['RR']
    integer     , optional, intent(in) :: smooth  !< fine-to-coarse discontinuity
                                                  !! smoothing {0,1,2} [0]

    type(SDC_Options)             :: opt_sdc
    type(HP_RefinementOptions_1D) :: opt_cf
    type(HP_CoarseningOptions_1D) :: opt_fc

    integer :: l, l_top

    ! spatial operators ........................................................

    call this % Init_ML_MeshOperators_3D(ml_mesh, po_x, nodes_x, smooth)

    ! temporal operators .......................................................

    l_top = size(ml_mesh % mesh)

    this % po_t = po_t
    this % ne_t = ne_t

    allocate(this % sdc      (1:l_top  ))
    allocate(this % iop_cf_t (1:l_top-1))
    allocate(this % iop_fc_t (2:l_top  ))
    allocate(this % pop_fc_t (2:l_top  ))

    if (present(nodes_t)) then
      opt_sdc % nodes = nodes_t
      opt_cf  % nodes = nodes_t
      opt_fc  % nodes = nodes_t
    end if

    if (present(smooth)) then
      opt_fc % smooth = smooth
    end if

    do l = 1, l_top

      ! SDC method
      opt_sdc % n_col = po_t(l) + 1
      this % sdc(l) = SDC_Method(opt_sdc)

      ! coarse-to-fine transfers operators
      if (l < l_top) then
        ! interpolation
        opt_cf % po_c = po_t(l)
        opt_cf % po_f = po_t(l+1)
        opt_cf % mode = Mode(po_t(l), ne_t(l), po_t(l+1), ne_t(l+1))
        this % iop_cf_t(l) = HP_RefinementOperator_1D(opt_cf)
      end if

      ! fine-to-coarse transfers operators
      if (l > 1) then
        ! interpolation
        opt_fc % po_f   = po_t(l)
        opt_fc % po_c   = po_t(l-1)
        opt_fc % method = 'I'
        opt_fc % mode   = Mode(po_t(l-1), ne_t(l-1), po_t(l), ne_t(l))
        this % iop_fc_t(l) = HP_CoarseningOperator_1D(opt_fc)
        ! L2-projection
        opt_fc % method = 'P'
        this % pop_fc_t(l) = HP_CoarseningOperator_1D(opt_fc)
      end if

    end do

  end subroutine Init_ML_SpacetimeOperators_3D

  !-----------------------------------------------------------------------------
  !> Refinement mode from given element orders and numbers

  pure integer function Mode(po_c, ne_c, po_f, ne_f)
    integer, intent(in) :: po_c !< order  of coarse elements
    integer, intent(in) :: ne_c !< number of coarse elements
    integer, intent(in) :: po_f !< order  of fine elements
    integer, intent(in) :: ne_f !< number of fine elements

    if (ne_f == ne_c) then
      ! cloning
      if (po_c /= po_f) then
        ! p-refinement
        Mode = 1
      else
        ! identity
        Mode = 0
      end if
    else if (ne_f == 2*ne_c) then
      ! subdivision
      Mode = 2
    else
      ! unknown refinement type
      Mode = -1
    end if

  end function Mode

  !=============================================================================

end module ML__Spacetime_Operators__3D
