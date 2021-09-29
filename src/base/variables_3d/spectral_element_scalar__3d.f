!> summary:  3D spectral element scalar variable
!> author:   Joerg Stiller
!> date:     2021/7/01
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Spectral_Element_Scalar__3D
  use Kind_Parameters, only: RNP
  use Spectral_Element_Mesh__3D
  use Spectral_Element_Variable__3D
  implicit none
  private

  public :: SpectralElementScalar_3D

  !-----------------------------------------------------------------------------
  !> 3D spectral element scalar variable

  type, extends(SpectralElementVariable_3D) :: SpectralElementScalar_3D
  contains
    ! wrappers eliminating the obsolete component dimension of the result
    procedure :: GetVolumeIntegral
    procedure :: GetTrace
  end type SpectralElementScalar_3D

  ! constructor
  interface SpectralElementScalar_3D
    module procedure New_SpectralElementScalar_3D
    module procedure WrapScalar
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D spectral element scalar constructor

  function New_SpectralElementScalar_3D(sem) result(this)
    class(SpectralElementMesh_3D), target, intent(in) :: sem
    type(SpectralElementScalar_3D) :: this

    call this % Init_SpectralElementVariable_3D(sem, nc = 1)

  end function New_SpectralElementScalar_3D

  !-----------------------------------------------------------------------------
  !> Wrap scalar mesh variable into spectral element scalar

  function WrapScalar(sem, u) result(this)
    class(SpectralElementMesh_3D), target, intent(in) :: sem
    real(RNP), contiguous, intent(in) :: u(0:,0:,0:,:)
    type(SpectralElementScalar_3D) :: this

    this % sem => sem
    call MapVariable(sem % std_op % po, sem % mesh % n_elem, u, this % val)

  end function WrapScalar

  !-----------------------------------------------------------------------------
  !> TBP for computing the volume integral of the scalar

  subroutine GetVolumeIntegral(this, vi)
    class(SpectralElementScalar_3D), intent(in) :: this
    real(RNP), intent(out) :: vi !< volume integral

    real(RNP) :: vi_all(1)

    call this % GetVolumeIntegrals(vi_all)
    vi = vi_all(1)

  end subroutine GetVolumeIntegral

  !-----------------------------------------------------------------------------
  !> Extract the trace of the scalar

  subroutine GetTrace(this, tr_val, align)
    class(SpectralElementScalar_3D), intent(in) :: this
    real(RNP), intent(inout) :: tr_val(:,:,:,:) !< trace of the scalar
    logical, optional, intent(in) :: align !< align traces with mesh face [F]

    real(RNP), pointer, save :: tr_all(:,:,:,:,:)

    !$omp single
    call MapTrace(size(tr_val,1), size(tr_val,4), tr_s=tr_val, tr_a=tr_all)
    !$omp end single

    call this % GetTraces(tr_all, align)

  end subroutine GetTrace

  !=============================================================================
  ! Helpers

  !-----------------------------------------------------------------------------
  !> Map scalar to array variable

  subroutine MapVariable(po, ne, s, a)
    integer, intent(in) :: po !< polynomial order
    integer, intent(in) :: ne !< number of elements
    real(RNP), target , intent(in)  :: s((po+1)**3 * ne) !< scalar variable
    real(RNP), pointer, intent(out) :: a(:,:,:,:,:)      !< array variable

    a(0:po,0:po,0:po,1:ne,1:1) => s

  end subroutine MapVariable

  !-----------------------------------------------------------------------------
  !> Map scalar to array trace

  subroutine MapTrace(np, ne, tr_s, tr_a)
    integer, intent(in) :: np, ne
    real(RNP), target , intent(in)  :: tr_s(np*np*6*ne)
    real(RNP), pointer, intent(out) :: tr_a(:,:,:,:,:)

    tr_a(1:np,1:np,1:6,1:ne,1:1) => tr_s

  end subroutine MapTrace

  !=============================================================================

end module Spectral_Element_Scalar__3D
