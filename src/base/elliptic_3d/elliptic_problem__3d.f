!> summary:  Interface to elliptic problems
!> author:   Joerg Stiller
!> date:     2019/01/27
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Elliptic_Problem__3D
  use Kind_Parameters, only: RNP

  implicit none
  private

  public :: EllipticProblem_3D

  !-----------------------------------------------------------------------------
  !> Abstract type for defining and handling an elliptic problem

  type, abstract :: EllipticProblem_3D

    real(RNP) :: lambda = 0  !< Helmholtz parameter
    real(RNP) :: nu_0   = 1  !< diffusivity mean value ν₀
    real(RNP) :: nu_1   = 0  !< diffusivity fluctuation amplitude ν₁
    real(RNP) :: d_nu   = 0  !< diffusivity fluctuation phase shift
    integer   :: k_nu   = 1  !< diffusivity fluctuation wave number

  contains

    procedure :: GetSource

    procedure(GetExactSolution),       deferred :: GetExactSolution
    procedure(GetExactGradient),       deferred :: GetExactGradient
    procedure(GetExactLaplacian),      deferred :: GetExactLaplacian
    procedure(GetDiffusivity),         deferred :: GetDiffusivity
    procedure(GetDiffusivityGradient), deferred :: GetDiffusivityGradient

  end type EllipticProblem_3D

  abstract interface

    !---------------------------------------------------------------------------
    !> Exact solution

    subroutine GetExactSolution(this, x, u)
      import
      class(EllipticProblem_3D), intent(in) :: this
      real(RNP), intent(in)  :: x(:,:,:,:,:) !< mesh points
      real(RNP), intent(out) :: u(:,:,:,:)   !< solution, u(x)
    end subroutine GetExactSolution

    !---------------------------------------------------------------------------
    !> Exact solution gradient

    subroutine GetExactGradient(this, x, grad_u)
      import
      class(EllipticProblem_3D), intent(in) :: this
      real(RNP), intent(in)  :: x(:,:,:,:,:)      !< mesh points
      real(RNP), intent(out) :: grad_u(:,:,:,:,:) !< ∇u(x)
    end subroutine GetExactGradient

    !---------------------------------------------------------------------------
    !> Exact laplacian

    subroutine GetExactLaplacian(this, x, laplace_u)
      import
      class(EllipticProblem_3D), intent(in) :: this
      real(RNP), intent(in)  :: x(:,:,:,:,:)       !< mesh points
      real(RNP), intent(out) :: laplace_u(:,:,:,:) !< ∇²u(x)
    end subroutine GetExactLaplacian

    !---------------------------------------------------------------------------
    !> Diffusivity

    subroutine GetDiffusivity(this, x, nu)
      import
      class(EllipticProblem_3D), intent(in) :: this
      real(RNP), intent(in)  :: x(:,:,:,:,:)  !< mesh points
      real(RNP), intent(out) :: nu(:,:,:,:)   !< ν(x)
    end subroutine GetDiffusivity

    !---------------------------------------------------------------------------
    !> Diffusivity gradient

    subroutine GetDiffusivityGradient(this, x, grad_nu)
      import
      class(EllipticProblem_3D), intent(in) :: this
      real(RNP), intent(in)  :: x(:,:,:,:,:)       !< mesh points
      real(RNP), intent(out) :: grad_nu(:,:,:,:,:) !< ∇ν(x)
    end subroutine GetDiffusivityGradient

  end interface

contains

  !-----------------------------------------------------------------------------
  !> Exact source

  subroutine GetSource(this, x, f)
    class(EllipticProblem_3D), intent(in) :: this
    real(RNP), intent(in)  :: x(:,:,:,:,:)  !< mesh points
    real(RNP), intent(out) :: f(:,:,:,:)    !< source, f(x)

    real(RNP), dimension(:,:,:,:),   allocatable, save :: nu, u, laplace_u
    real(RNP), dimension(:,:,:,:,:), allocatable, save :: grad_nu, grad_u

    integer :: np, ne
    integer :: e, i, j, k

    allocate(nu        , mold=f)
    allocate(u         , mold=f)
    allocate(laplace_u , mold=f)
    allocate(grad_nu   , mold=x)
    allocate(grad_u    , mold=x)

    call this % GetExactSolution       (x, u)
    call this % GetExactGradient       (x, grad_u)
    call this % GetExactLaplacian      (x, laplace_u)
    call this % GetDiffusivity         (x, nu)
    call this % GetDiffusivityGradient (x, grad_nu)

    np = size(x,1)
    ne = size(x,4)

    associate(lambda => this % lambda)

      !$omp do
      do e = 1, ne
        do k = 1, np
        do j = 1, np
        do i = 1, np

          f(i,j,k,e) = lambda * u(i,j,k,e)                        &
                      - ( nu(i,j,k,e) * laplace_u(i,j,k,e)        &
                        + grad_nu(i,j,k,e,1) * grad_u(i,j,k,e,1)  &
                        + grad_nu(i,j,k,e,2) * grad_u(i,j,k,e,2)  &
                        + grad_nu(i,j,k,e,3) * grad_u(i,j,k,e,3) )
        end do
        end do
        end do
      end do

    end associate

    deallocate(nu, u, laplace_u, grad_nu, grad_u)

  end subroutine GetSource

  !=============================================================================

end module Elliptic_Problem__3D
