!> summary:  3D multilevel spectral element space operators
!> author:   Joerg Stiller
!> date:     2024/06/26
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ML__Mesh_Operators__3D
  use Kind_Parameters
  use Constants
  use XMPI
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
    procedure :: Get_MeshCharacteristics
    procedure :: Print_MeshCharacteristics
  end type ML_MeshOperators_3D

  ! constructor interface
  interface ML_MeshOperators_3D
    procedure New_ML_MeshOperators_3D
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Constructor of 3D multilevel mesh operators

  function New_ML_MeshOperators_3D(ml_mesh, po, nodes, smooth, filter) &
      result(this)
    class(ML_Mesh_3D),   intent(in) :: ml_mesh !< multilevel mesh partition
    integer,             intent(in) :: po(:)   !< sequence of polynomial orders
    character, optional, intent(in) :: nodes   !< 'G' or 'L' ['L']
    integer,   optional, intent(in) :: smooth  !< fine-to-coarse discontinuity
                                               !! smoothing {0,1,2} [0]
    integer,   optional, intent(in) :: filter  !< fine-to-coarse filter order,
                                               !! default: no filtering
    type(ML_MeshOperators_3D) :: this

    call Init_ML_MeshOperators_3D(this, ml_mesh, po, nodes, smooth, filter)

  end function New_ML_MeshOperators_3D

  !-----------------------------------------------------------------------------
  !> Initialization of 3D multilevel mesh operators

  subroutine Init_ML_MeshOperators_3D(this, ml_mesh, po, nodes, smooth, filter)
    class(ML_MeshOperators_3D), intent(inout) :: this
    class(ML_Mesh_3D),   intent(in) :: ml_mesh !< multilevel mesh partition
    integer,             intent(in) :: po(:)   !< sequence of polynomial orders
    character, optional, intent(in) :: nodes   !< 'G' or 'L' ['L']
    integer,   optional, intent(in) :: smooth  !< fine-to-coarse discontinuity
                                               !! smoothing {0,1,2} [0]
    integer,   optional, intent(in) :: filter  !< fine-to-coarse filter order,
                                               !! default: no filtering

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

      if (present(filter)) then
        opt_fc % filter = filter
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

  !-----------------------------------------------------------------------------
  !> Get characteristics of the multilevel spectral element mesh

  subroutine Get_MeshCharacteristics(this, n_elem, n_active, n_leaf, emq)

    class(ML_MeshOperators_3D), intent(in) :: this

    integer, allocatable, optional, intent(inout) :: n_elem(:,:)
      !< n_elem(l_top,4), element counts:
      !! - n_elem(l,1) = loc num in level l
      !! - n_elem(l,2) = min num per partition in level l
      !! - n_elem(l,3) = max num per partition in level l
      !! - n_elem(l,4) = tot num in level l

    integer, allocatable, optional, intent(inout) :: n_active(:,:)
      !< n_active(l_top,4), active element counts:
      !! - n_active(l,1) = loc num in level l
      !! - n_active(l,2) = min num per partition in level l
      !! - n_active(l,3) = max num per partition in level l
      !! - n_active(l,4) = tot num in level l

    integer, allocatable, optional, intent(inout) :: n_leaf(:,:)
      !< n_leaf(l_top,4), leaf element counts:
      !! - n_leaf(l,1) = loc num in level l
      !! - n_leaf(l,2) = min num per partition in level l
      !! - n_leaf(l,3) = max num per partition in level l
      !! - n_leaf(l,4) = tot num in level l

    real(RNP), allocatable, optional, intent(inout) :: emq(:,:)
      !< emq(l_top,6), active element metrics and quality
      !! - emq(l,1) = min size in level l (cube root of V)
      !! - emq(l,2) = max size in level l (cube root of V)
      !! - emq(l,3) = min aspect ratio in level l
      !! - emq(l,4) = max aspect ratio in level l
      !! - emq(l,5) = min scaled Jacobian level l
      !! - emq(l,6) = max scaled Jacobian level l

    real(RNP), save :: dx_min, dx_max, ar_min, ar_max, qj_min, qj_max

    real(RNP) :: dx_x(3)
    real(RNP) :: dx_e, ar_e, qj_e
    integer   :: e, l, l_top

    associate(comm => this % sem(1) % mesh % comm_world)

      l_top = size(this % sem)

      ! number of elements .....................................................

      if (present(n_elem)) then

        !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        !$omp master

        if (allocated(n_elem)) then
          if (any(shape(n_elem) /= [l_top,4])) deallocate(n_elem)
        end if
        if (.not. allocated(n_elem)) then
          allocate(n_elem(l_top,4))
        end if

        do l = 1, l_top
          n_elem(l,1) = max(0, this % sem(l) % mesh % n_elem)
        end do

        call XMPI_Allreduce(n_elem(:,1), n_elem(:,2), MPI_MIN, comm)
        call XMPI_Allreduce(n_elem(:,1), n_elem(:,3), MPI_MAX, comm)
        call XMPI_Allreduce(n_elem(:,1), n_elem(:,4), MPI_SUM, comm)

        !$omp end master
        !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      end if

      ! number of active elements ..............................................

      if (present(n_active)) then

        !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        !$omp master

        if (allocated(n_active)) then
          if (any(shape(n_active) /= [l_top,4])) deallocate(n_active)
        end if
        if (.not. allocated(n_active)) then
          allocate(n_active(l_top,4))
        end if

        do l = 1, l_top
          n_active(l,1) = max(0, this % sem(l) % mesh % n_elem_active)
        end do

        call XMPI_Allreduce(n_active(:,1), n_active(:,2), MPI_MIN, comm)
        call XMPI_Allreduce(n_active(:,1), n_active(:,3), MPI_MAX, comm)
        call XMPI_Allreduce(n_active(:,1), n_active(:,4), MPI_SUM, comm)

        !$omp end master
        !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      end if

      ! number of leaf elements ................................................

      if (present(n_leaf)) then

        !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        !$omp master

        if (allocated(n_leaf)) then
          if (any(shape(n_leaf) /= [l_top,4])) deallocate(n_leaf)
        end if
        if (.not. allocated(n_leaf)) then
          allocate(n_leaf(l_top,4))
        end if

        do l = 1, l_top
          associate(mesh => this % sem(l) % mesh)
            if (mesh % n_elem_active > 0) then
              n_leaf(l,1) = count(mesh % element % IsLeaf())
            else
              n_leaf(l,1) = 0
            end if
          end associate
        end do

        call XMPI_Allreduce(n_leaf(:,1), n_leaf(:,2), MPI_MIN, comm)
        call XMPI_Allreduce(n_leaf(:,1), n_leaf(:,3), MPI_MAX, comm)
        call XMPI_Allreduce(n_leaf(:,1), n_leaf(:,4), MPI_SUM, comm)

        !$omp end master
        !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      end if

      ! element metrics and quality ............................................

      if (present(emq)) then

        !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        !$omp master

        if (allocated(emq)) then
          if (any(shape(emq) /= [l_top,5])) deallocate(emq)
        end if
        if (.not. allocated(emq)) then
          allocate(emq(l_top,6))
        end if

        !$omp end master
        !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

        do l = 1, l_top
          associate( mesh => this % sem(l) % mesh         &
                   , Jd   => this % sem(l) % metrics % Jd )

            !$omp master
            dx_min =  huge(ONE)
            dx_max = -huge(ONE)
            ar_min =  huge(ONE)
            ar_max = -huge(ONE)
            qj_min =  huge(ONE)
            qj_max = -huge(ONE)
            !$omp end master

            !$omp do reduction(min:dx_min,ar_min,qj_min) &
            !$omp &  reduction(max:dx_max,ar_max,qj_max)
            do e = 1, mesh % n_elem_active

              call mesh % element(e) % GetCuboidDimensions(dx_x)

              dx_e = product(dx_x) ** THIRD
              ar_e = maxval(dx_x) / minval(dx_x)
              qj_e = minval(Jd(:,:,:,e)) / maxval(Jd(:,:,:,e))

              dx_min = min(dx_min, dx_e)
              dx_max = max(dx_max, dx_e)
              ar_min = min(ar_min, ar_e)
              ar_max = max(ar_max, ar_e)
              qj_min = min(qj_min, qj_e)
              qj_max = max(qj_max, qj_e)

            end do

            !$omp master
            call XMPI_Allreduce(dx_min, emq(l,1), MPI_MIN, comm)
            call XMPI_Allreduce(dx_max, emq(l,2), MPI_MAX, comm)
            call XMPI_Allreduce(ar_min, emq(l,3), MPI_MAX, comm)
            call XMPI_Allreduce(ar_max, emq(l,4), MPI_MAX, comm)
            call XMPI_Allreduce(qj_min, emq(l,5), MPI_MIN, comm)
            call XMPI_Allreduce(qj_max, emq(l,6), MPI_MAX, comm)
            !$omp end master

          end associate
        end do

      end if
    end associate

  end subroutine Get_MeshCharacteristics

  !-----------------------------------------------------------------------------
  !> Print characteristics of the multilevel spectral element mesh

  subroutine Print_MeshCharacteristics(this)
    class(ML_MeshOperators_3D), intent(in) :: this

    integer  , allocatable, save :: n_elem(:,:), n_active(:,:), n_leaf(:,:)
    real(RNP), allocatable, save :: emq(:,:)

    integer(IXL) :: np_leaf, np_tot
    integer :: l

    call this % Get_MeshCharacteristics(n_elem, n_active, n_leaf, emq)

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !$omp master

    if (this % sem(1) % mesh % proc == 0) then

      write(*,'(/,A4,X,6A9,A4,2(3X,A9),2(X,A9))') &
          '   l'      , &
          '  n_parts' , &
          '   na_min' , &
          '   na_max' , &
          '   na_tot' , &
          '   ne_tot' , &
          '  ne_leaf' , &
          '  po'      , &
          '   dx_min' , &
          '   dx_max' , &
          '   ar_max' , &
          '   qj_min'

      np_tot  = 0
      np_leaf = 0

      do l = 1, size(this % sem)
        associate( mesh => this % sem(l) % mesh        &
                 , po   => this % sem(l) % std_op % po )

        np_tot  = np_tot  + n_elem(l,4) * (po + 1)**3
        np_leaf = np_leaf + n_leaf(l,4) * (po + 1)**3

        write(*,'(I4,X,6I9,I4,2(2X,ES10.3),2F10.3)') &
            l             , &
            mesh%n_parts  , &
            n_active(l,2) , &
            n_active(l,3) , &
            n_active(l,4) , &
            n_elem(l,4)   , &
            n_leaf(l,4)   , &
            po            , &
            emq(l,1:2)    , &
            emq(l,4:5)

        end associate
      end do

      write(*,*)
      write(*,'(2X,A)') 'global'
      write(*,'(T5,A,I0)') 'ne_tot  = ', sum(n_elem(:,4))
      write(*,'(T5,A,I0)') 'ne_leaf = ', sum(n_leaf(:,4))
      write(*,'(T5,A,I0)') 'np_tot  = ', np_tot
      write(*,'(T5,A,I0)') 'np_leaf = ', np_leaf

    end if

    if (allocated( n_elem   )) deallocate( n_elem   )
    if (allocated( n_active )) deallocate( n_active )
    if (allocated( n_leaf   )) deallocate( n_leaf   )
    if (allocated( emq      )) deallocate( emq      )

    !$omp end master
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  end subroutine Print_MeshCharacteristics

  !=============================================================================

end module ML__Mesh_Operators__3D
