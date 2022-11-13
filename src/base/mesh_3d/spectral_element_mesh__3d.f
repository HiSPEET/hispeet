!> summary:  3D spectral element mesh
!> author:   Joerg Stiller
!> date:     2021/06/29
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Spectral_Element_Mesh__3D
  use Kind_Parameters, only: RNP
  use Constants      , only: ZERO
  use XMPI
  use Standard_Operators__1D
  use Mesh__3D
  use Mesh_Metrics__3D
  implicit none
  private

  public :: SpectralElementMesh_3D

  !-----------------------------------------------------------------------------
  !> 3D spectral element mesh

  type SpectralElementMesh_3D
    class(Mesh_3D), pointer    :: mesh    !< mesh partition
    type(StandardOperators_1D) :: std_op  !< standard element operators
    type(MeshMetrics_3D)       :: metrics !< mesh points and metrics
  contains
    procedure :: Init_SpectralElementMesh_3D
    procedure :: Get_Volume
    procedure :: Get_SurfaceAreas
    procedure :: Get_DG_DiagonalMassMatrix
  end type SpectralElementMesh_3D

  ! constructor interface
  interface SpectralElementMesh_3D
    module procedure New_SpectralElementMesh_3D
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D spectral element mesh constructor

  function New_SpectralElementMesh_3D(mesh, po, basis) result(this)
    class(Mesh_3D), target, intent(in) :: mesh  !< mesh partition
    integer,                intent(in) :: po    !< polynomial order
    character,    optional, intent(in) :: basis !< 'G' or 'L' ['L']

    type(SpectralElementMesh_3D) :: this

    call Init_SpectralElementMesh_3D(this, mesh, po, basis)

  end function New_SpectralElementMesh_3D

  !-----------------------------------------------------------------------------
  !> 3D spectral element mesh initialization

  subroutine Init_SpectralElementMesh_3D(this, mesh, po, basis)
    class(SpectralElementMesh_3D), intent(inout) :: this
    class(Mesh_3D),        target, intent(in)    :: mesh  !< mesh partition
    integer,                       intent(in)    :: po    !< polynomial order
    character,           optional, intent(in)    :: basis !< 'G' or 'L' ['L']

    this % mesh    => mesh
    this % std_op  =  StandardOperators_1D(po, basis)
    this % metrics =  MeshMetrics_3D(mesh, this % std_op)

  end subroutine Init_SpectralElementMesh_3D

  !-----------------------------------------------------------------------------
  !> TBP for computing the volume of the computational domain

  subroutine Get_Volume(this, vol)
    class(SpectralElementMesh_3D), intent(in) :: this
    real(RNP), intent(out) :: vol !< volume, should be PRIVATE with OpenMP

    real(RNP), save :: v_loc, v_glob

    real(RNP), allocatable :: www(:,:,:)
    integer :: e, i, j, k

    if (this % mesh % part < 0) then
      vol = 0
      return
    end if

    associate( mesh   => this % mesh         &
             , std_op => this % std_op       &
             , Jd     => this % metrics % Jd )

      if (mesh % regular) then

        !$omp master
        v_loc = product(mesh % dx)  *  mesh % n_elem
        !$omp end master

      else

        ! precompute 3D quadrature weights
        allocate(www(0:std_op%po, 0:std_op%po, 0:std_op%po))
        do k = 0, std_op%po
        do j = 0, std_op%po
        do i = 0, std_op%po
          www(i,j,k) = std_op % w(i) * std_op % w(j) * std_op % w(k)
        end do
        end do
        end do

        !$omp single
        v_loc = ZERO
        !$omp end single

        !$omp do reduction(+:v_loc) schedule(static)
        do e = 1, mesh % n_elem
          v_loc = v_loc + sum(www * Jd(:,:,:,e))
        end do

      end if

      !$omp master
      call XMPI_Allreduce(v_loc, v_glob, MPI_SUM, mesh % comm_parts)
      !$omp end master
      !$omp barrier

      vol = v_glob

    end associate

  end subroutine Get_Volume

  !-----------------------------------------------------------------------------
  !> TBP for computing the areas of the boundary surfaces

  subroutine Get_SurfaceAreas(this, area)
    class(SpectralElementMesh_3D), intent(in) :: this
    real(RNP), intent(out) :: area(:) !< areas, should be PRIVATE with OpenMP

    real(RNP), allocatable, save :: a_loc(:), a_glob(:)

    if (this % mesh % part < 0) then
      area = 0
      return
    end if

    !$omp master
    allocate(a_loc (this % mesh % n_bound), source = ZERO)
    allocate(a_glob(this % mesh % n_bound))
    !$omp end master
    !$omp barrier

    if (this % mesh % regular) then
      call Get_SurfaceAreas_R(this, a_loc)
    else
      call Get_SurfaceAreas_G(this, a_loc)
    end if

    !$omp master
    call XMPI_Allreduce(a_loc, a_glob, MPI_SUM, this % mesh % comm_parts)
    !$omp end master
    !$omp barrier

    area = a_glob

  end subroutine Get_SurfaceAreas

  !-----------------------------------------------------------------------------
  !> Computation of the surface areas for a regular mesh

  subroutine Get_SurfaceAreas_R(this, area)
    class(SpectralElementMesh_3D), intent(in) :: this
    real(RNP), intent(inout) :: area(this%mesh%n_bound)

    real(RNP), save :: a_face(6)
    real(RNP) :: a_priv
    integer :: b, f

    associate(mesh => this % mesh)

      if (mesh % part < 0) return

      a_face(1:2) = mesh % dx(2) * mesh % dx(3)
      a_face(3:4) = mesh % dx(3) * mesh % dx(1)
      a_face(5:6) = mesh % dx(1) * mesh % dx(2)

      do b = 1, mesh % n_bound
        a_priv = ZERO
        !$omp do
        do f = 1, mesh % boundary(b) % n_face
          a_priv = a_priv + a_face(mesh % boundary(b) % face(f) % element_face)
        end do
        !$omp end do nowait
        !$omp atomic
        area(b) = area(b) + a_priv
        !$omp barrier
      end do

    end associate

  end subroutine Get_SurfaceAreas_R

  !-----------------------------------------------------------------------------
  !> Computation of the surface areas for a general mesh

  subroutine Get_SurfaceAreas_G(this, area)
    class(SpectralElementMesh_3D), intent(in) :: this
    real(RNP), intent(inout) :: area(this%mesh%n_bound)

    real(RNP), allocatable :: ww(:,:)
    real(RNP) :: a_priv
    integer :: b, e, f, i, j, s

    associate( mesh   => this % mesh        &
             , std_op => this % std_op      &
             , a      => this % metrics % a )

      if (mesh % part < 0) return

      ! precompute 2D quadrature weights
      allocate(ww(0:std_op%po, 0:std_op%po))
      do j = 0, std_op%po
      do i = 0, std_op%po
        ww(i,j) = std_op % w(i) * std_op % w(j)
      end do
      end do

      do b = 1, mesh % n_bound
        a_priv = ZERO
        !$omp do
        do f = 1, mesh % boundary(b) % n_face
          e = mesh % boundary(b) % face(f) % element_id   ! element ID
          s = mesh % boundary(b) % face(f) % element_face ! element side
          a_priv = a_priv + sum(ww * a(:,:,s,e))
        end do
        !$omp end do nowait
        !$omp atomic
        area(b) = area(b) + a_priv
        !$omp barrier
      end do

    end associate

  end subroutine Get_SurfaceAreas_G

  !-----------------------------------------------------------------------------
  !> TBP to compute the quadrature-based diagonal mass matrix for DG-SEM

  subroutine Get_DG_DiagonalMassMatrix(this, mm)
    class(SpectralElementMesh_3D), intent(in) :: this
    real(RNP), intent(out) :: mm(:,:,:,:) !< diagonal entries of mass matrix

    real(RNP), allocatable :: www(:,:,:)
    integer :: e, i, j, k

    ! skip computation for empty partitions
    if (size(mm) == 0) return

    associate( mesh   => this % mesh         &
             , std_op => this % std_op       &
             , Jd     => this % metrics % Jd )

      allocate(www(0:std_op%po, 0:std_op%po, 0:std_op%po))
      do k = 0, std_op%po
      do j = 0, std_op%po
      do i = 0, std_op%po
        www(i,j,k) = std_op % w(i) * std_op % w(j) * std_op % w(k)
      end do
      end do
      end do

      if (mesh % regular) then

        !$omp do
        do e = 1, mesh % n_elem
          mm(:,:,:,e) = www * Jd(:,:,:,1)
        end do

      else

        !$omp do
        do e = 1, mesh % n_elem
          mm(:,:,:,e) = www * Jd(:,:,:,e)
        end do

      end if
    end associate

  end subroutine Get_DG_DiagonalMassMatrix

  !=============================================================================

end module Spectral_Element_Mesh__3D
