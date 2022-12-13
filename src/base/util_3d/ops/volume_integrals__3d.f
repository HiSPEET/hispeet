!> summary:  Volume integrals of 3D spectral element variables
!> author:   Joerg Stiller
!> date:     2022/11/11
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Volume_Integrals__3D
  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO
  use XMPI
  use Spectral_Element_Mesh__3D
  implicit none

  private

  public :: GetVolumeIntegral

  interface GetVolumeIntegral
    module procedure GetVolumeIntegral_S
    module procedure GetVolumeIntegral_A
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Compute the volume integral of a scalar variable

  subroutine GetVolumeIntegral_S(sem, u, int_u)
    class(SpectralElementMesh_3D), intent(in) :: sem
    real(RNP), contiguous, intent(in) :: u(:,:,:,:)
    real(RNP), intent(out) :: int_u

    real(RNP) :: int_u_(1)

    if (sem % mesh % part >= 0) then
      call GetVolumeIntegral_X( sem                 &
                              , sem % std_op % po   &
                              , sem % mesh % n_elem &
                              , 1, u, int_u_        )
      int_u = int_u_(1)
    else
      int_u = 0
    end if

  end subroutine GetVolumeIntegral_S

  !-----------------------------------------------------------------------------
  !> Compute the volume integral of an array variable

  subroutine GetVolumeIntegral_A(sem, u, int_u)
    class(SpectralElementMesh_3D), intent(in) :: sem
    real(RNP), contiguous, intent(in) :: u(:,:,:,:,:)
    real(RNP), intent(out) :: int_u(:)

    if (sem % mesh % part >= 0) then
      call GetVolumeIntegral_X( sem                 &
                              , sem % std_op % po   &
                              , sem % mesh % n_elem &
                              , size(u,5), u, int_u )
    else
      int_u = 0
    end if

  end subroutine GetVolumeIntegral_A

  !-----------------------------------------------------------------------------
  !> Computation of volume integrals, explicit shape

  subroutine GetVolumeIntegral_X(sem, po, ne, nc, u, int_u)
    class(SpectralElementMesh_3D), intent(in) :: sem
    integer,   intent(in)  :: po
    integer,   intent(in)  :: ne
    integer,   intent(in)  :: nc
    real(RNP), intent(in)  :: u(0:po,0:po,0:po,ne,nc)
    real(RNP), intent(out) :: int_u(nc)

    real(RNP), allocatable, save :: int_u_loc(:), int_u_glob(:)
    real(RNP), allocatable :: int_u_priv(:), www(:,:,:)
    real(RNP) :: Jd0
    integer   :: c, e, i, j, k

    associate( mesh   => sem % mesh          &
             , std_op => sem % std_op        &
             , Jd     => sem % metrics % Jd  )

      allocate(int_u_priv(nc), source = ZERO)

      !$omp master
      if (allocated(int_u_glob)) deallocate(int_u_glob)
      if (allocated(int_u_loc )) deallocate(int_u_loc)

      allocate(int_u_glob(nc), source = ZERO)
      allocate(int_u_loc (nc), source = ZERO)
      !$omp end master
      !$omp barrier

      ! precompute 3D quadrature weights
      allocate(www(0:po, 0:po, 0:po))
      do k = 0, po
      do j = 0, po
      do i = 0, po
        www(i,j,k) = std_op % w(i) * std_op % w(j) * std_op % w(k)
      end do
      end do
      end do

      if (mesh % regular) then

        Jd0 = product(mesh % dx) / 8

        !$omp do schedule(static)
        do e = 1, mesh % n_elem
          do c = 1, nc
            int_u_priv(c) = int_u_priv(c) + Jd0 * sum(www * u(:,:,:,e,c))
          end do
        end do
        !$omp end do nowait

      else

        !$omp do schedule(static)
        do e = 1, mesh % n_elem
          do c = 1, nc
            int_u_priv(c) = int_u_priv(c) &
                          + sum(www * Jd(:,:,:,e) * u(:,:,:,e,c))
          end do
        end do
        !$omp end do nowait

      end if

      !$omp atomic
      int_u_loc = int_u_loc + int_u_priv
      !$omp barrier

      !$omp master
      call XMPI_Allreduce(int_u_loc, int_u_glob, MPI_SUM, mesh % comm_parts)
      !$omp end master
      !$omp barrier

      int_u = int_u_glob

    end associate

  end subroutine GetVolumeIntegral_X

  !=============================================================================

end module Volume_Integrals__3D
