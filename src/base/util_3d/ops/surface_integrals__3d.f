!> summary:  Surface integrals of 3D spectral element boundary variables
!> author:   Joerg Stiller
!> date:     2022/11/11
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Surface_Integrals__3D
  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO
  use Execution_Control, only: Error
  use XMPI
  use Spectral_Element_Mesh__3D
  use Boundary_Variable__3D
  implicit none

  private

  public :: GetSurfaceIntegral
  public :: GetSurfaceIntegrals

contains

  !-----------------------------------------------------------------------------
  !> Compute the integral of variables on a single boundary surface

  subroutine GetSurfaceIntegral(sem, bv_u, int_ub, leaf)
    class(SpectralElementMesh_3D), intent(in)  :: sem
    class(BoundaryVariable_3D),    intent(in)  :: bv_u
    real(RNP),                     intent(out) :: int_ub(:)
    logical,             optional, intent(in)  :: leaf !< constrain to leaves

    real(RNP), allocatable, save :: int_ub_glob(:), int_ub_loc(:)

    real(RNP), allocatable :: ww(:,:)
    integer :: i, j

    if (sem % mesh % part < 0) then
      int_ub = 0
      return
    end if

    associate(po => bv_u % po, nc => bv_u % nc, w => sem % std_op % w)

      !$omp master
      if (sem % std_op % po /= po .or. size(int_ub) /= nc) then
        call Error( 'GetSurfaceIntegral'     &
                  , 'arguments do not match' &
                  , 'Surface_Integrals__3D'  )
      end if

      if (allocated(int_ub_glob)) deallocate(int_ub_glob)
      if (allocated(int_ub_loc )) deallocate(int_ub_loc)

      allocate(int_ub_loc (nc), source = ZERO)
      allocate(int_ub_glob(nc), source = ZERO)
      !$omp end master
      !$omp barrier

      ! 2D quadrature weights
      allocate(ww(0:po,0:po))
      do j = 0, po
      do i = 0, po
        ww(i,j) = w(i) * w(j)
      end do
      end do

      call GetLocalSurfaceIntegral(sem, bv_u, ww, int_ub_loc, leaf)

      !$omp master
      call XMPI_Allreduce(int_ub_loc, int_ub_glob, MPI_SUM, sem%mesh%comm_parts)
      !$omp end master
      !$omp barrier

      int_ub = int_ub_glob

    end associate

  end subroutine GetSurfaceIntegral

  !-----------------------------------------------------------------------------
  !> Compute the integral of a variable on multiple boundary surfaces

  subroutine GetSurfaceIntegrals(sem, bv_u, int_ub, leaf)
    class(SpectralElementMesh_3D), intent(in)  :: sem
    class(BoundaryVariable_3D),    intent(in)  :: bv_u(:)
    real(RNP),                     intent(out) :: int_ub(:,:)
    logical,             optional, intent(in)  :: leaf !< constrain to leaves

    real(RNP), allocatable, save :: int_ub_glob(:,:), int_ub_loc(:,:)

    real(RNP), allocatable :: ww(:,:)
    integer :: b, i, j, po, nc, nb

    if (sem % mesh % part < 0 .or. size(bv_u) == 0) then
      int_ub = 0
      return
    end if

    po = sem % std_op % po
    nc = size(int_ub,1)
    nb = size(int_ub,2)

    associate(w => sem % std_op % w)

      !$omp master
      if ( size(bv_u) /= nb .or. any(bv_u%po /= po) .or. any(bv_u%nc /= nc)) then
        call Error( 'GetSurfaceIntegrals'    &
                  , 'arguments do not match' &
                  , 'Surface_Integrals__3D'  )
      end if

      if (allocated(int_ub_glob)) deallocate(int_ub_glob)
      if (allocated(int_ub_loc )) deallocate(int_ub_loc)

      allocate(int_ub_loc (nc,nb), source = ZERO)
      allocate(int_ub_glob(nc,nb), source = ZERO)
      !$omp end master
      !$omp barrier

      ! 2D quadrature weights
      allocate(ww(0:po,0:po))
      do j = 0, po
      do i = 0, po
        ww(i,j) = w(i) * w(j)
      end do
      end do

      do b = 1, nb
        call GetLocalSurfaceIntegral(sem, bv_u(b), ww, int_ub_loc(:,b), leaf)
      end do

      !$omp master
      call XMPI_Allreduce(int_ub_loc, int_ub_glob, MPI_SUM, sem%mesh%comm_parts)
      !$omp end master
      !$omp barrier

      int_ub = int_ub_glob

    end associate

  end subroutine GetSurfaceIntegrals

  !-----------------------------------------------------------------------------
  !> Computation of local surface integrals

  subroutine GetLocalSurfaceIntegral(sem, bv_u, ww, int_ub_loc, leaf)
    class(SpectralElementMesh_3D), intent(in)    :: sem
    class(BoundaryVariable_3D),    intent(in)    :: bv_u
    real(RNP),                     intent(in)    :: ww(:,:)
    real(RNP),                     intent(inout) :: int_ub_loc(:)
    logical,             optional, intent(in)    :: leaf !< constrain to leaves

    real(RNP) :: int_ub_priv(bv_u%nc)
    logical   :: complete
    integer   :: c, e, f, s

    if (present(leaf)) then
      complete = .not. leaf
    else
      complete = .true.
    end if

    associate( a => sem % metrics % a, boundary => bv_u % boundary)

      int_ub_priv = ZERO

      !$omp do
      do f = 1, boundary % n_face
        e = boundary % face(f) % element_id   ! element ID
        s = boundary % face(f) % element_face ! element side

        if (complete .or. sem % mesh % element(e) % IsLeaf()) then
          do c = 1, bv_u % nc
            int_ub_priv(c) = int_ub_priv(c) &
                           + sum(ww * a(:,:,s,e) * bv_u % val(:,:,f,c))
          end do
        end if
      end do
      !$omp end do nowait

      do c = 1, bv_u%nc
        !$omp atomic
        int_ub_loc(c) = int_ub_loc(c) + int_ub_priv(c)
      end do
      !$omp barrier

    end associate

  end subroutine GetLocalSurfaceIntegral

  !=============================================================================

end module Surface_Integrals__3D
