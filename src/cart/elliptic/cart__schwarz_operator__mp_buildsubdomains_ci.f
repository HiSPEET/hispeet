!> summary:  Build Schwarz subdomains for constant isotropic coefficients
!> author:   Joerg Stiller
!> date:     2018/11/08
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Build Schwarz subdomains for constant isotropic coefficients
!===============================================================================

submodule(CART__Schwarz_Operator) MP_BuildSubdomains_CI
  use Constants,     only: ZERO, ONE, HALF
  use Eigenproblems, only: SolveGeneralizedEigenproblem
  use Schwarz_Weighting
  use IP_Element_Operators_1D
  implicit none

contains

!-------------------------------------------------------------------------------
!> Set subdomain configurations and inverse 3D eigenvalues: const isotropic

module subroutine BuildSubdomains_CI(this, mesh, lambda, nu, bc)
  class(SchwarzOperator3D),   intent(inout) :: this   !< Schwarz operator
  class(MeshPartition),       intent(in)    :: mesh   !< mesh partition
  real(RNP),                  intent(in)    :: lambda !< Helmholtz parameter
  real(RNP),                  intent(in)    :: nu     !< diffusivity
  character,                  intent(in)    :: bc(:)  !< BC {'D','N'}

  character :: bc_face(size(bc))
  real(RNP) :: g0, g1, g2, g3
  integer   :: e, i, j, k, n1, n2, n3

  ! preliminaries ..............................................................

  n1 = size(this % V1, 1)
  n2 = size(this % V2, 1)
  n3 = size(this % V3, 1)

  if (allocated(this%cfg)) then
    if (any(shape(this%cfg) /= [ 3, mesh%ne ])) then
      deallocate(this%cfg)
    end if
  end if

  if (allocated(this%D_inv)) then
    if (any(shape(this%D_inv) /= [ n1, n2, n3, mesh%ne ])) then
      deallocate(this%D_inv)
    end if
  end if

  if (.not. allocated(this%cfg)) then
    allocate(this%cfg(3, mesh%ne))
  end if

  if (.not. allocated(this%D_inv)) then
    allocate(this%D_inv(n1, n2, n3, mesh%ne))
  end if

  ! cfg and D_inv ..............................................................

  associate( V1 => this % V1,  cfg   => this % cfg   &
           , V2 => this % V2,  D_inv => this % D_inv &
           , V3 => this % V3,  dx    => mesh % dx    )

    do e = 1, mesh%ne

      where(mesh % element(e) % face % boundary > 0)
        bc_face = bc( mesh % element(e) % face % boundary )
      elsewhere
        bc_face = ''
      end where

      cfg(1,e) = ConfigurationID(bc_face(1:2))
      cfg(2,e) = ConfigurationID(bc_face(3:4))
      cfg(3,e) = ConfigurationID(bc_face(5:6))

      g0 = dx(1) * dx(2) * dx(3)
      g1 = dx(2) * dx(3) / dx(1)
      g2 = dx(3) * dx(1) / dx(2)
      g3 = dx(1) * dx(2) / dx(3)

      do k = 1, n3
      do j = 1, n2
      do i = 1, n1

        D_inv(i,j,k,e) = ONE / ( lambda * g0                  &
                               + nu * ( g1 * V1(i, cfg(1,e))  &
                                      + g2 * V2(j, cfg(2,e))  &
                                      + g3 * V3(k, cfg(3,e))  &
                                      )                       &
                               )
      end do
      end do
      end do

    end do

  end associate

end subroutine BuildSubdomains_CI

!===============================================================================

end submodule MP_BuildSubdomains_CI
