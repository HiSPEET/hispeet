!> summary:  Build Schwarz subdomains for variable isotropic coefficients
!> author:   Joerg Stiller
!> date:     2018/12/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Build Schwarz subdomains for variable isotropic coefficients
!===============================================================================

submodule(CART__Schwarz_Operator) MP_BuildSubdomains_VI
  use Constants, only: ONE
  implicit none

contains

!-------------------------------------------------------------------------------
!> Set subdomain configurations and inverse 3D eigenvalues: variable isotropic

module subroutine BuildSubdomains_VI(this, eop, mesh, lambda, nu, bc)
  class(SchwarzOperator3D),   intent(inout) :: this        !< Schwarz operator
  class(StandardOperators1D), intent(in)    :: eop         !< 1D SE operators
  class(MeshPartition),       intent(in)    :: mesh        !< mesh partition
  real(RNP),                  intent(in)    :: lambda      !< Helmholtz parameter
  real(RNP),                  intent(in)    :: nu(:,:,:,:) !< diffusivity
  character,                  intent(in)    :: bc(:)       !< BC {'D','N'}

  character :: bc_face(size(bc))
  integer   :: e, i, j, k, n1, n2, n3, np
  real(RNP) :: g0, g1, g2, g3, nu_0
  real(RNP), allocatable :: nu_3(:), nu_23(:,:)

  ! preliminaries ..............................................................

  n1 = size(this % V1, 1)
  n2 = size(this % V2, 1)
  n3 = size(this % V3, 1)
  np = eop % po + 1

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

  allocate(nu_3(np), nu_23(np,np))

  ! cfg and D_inv ..............................................................

  associate( V1    => this % V1,     &
             V2    => this % V2,     &
             V3    => this % V3,     &
             cfg   => this % cfg,    &
             D_inv => this % D_inv,  &
             dx    => mesh % dx      )

    do e = 1, mesh%ne

      associate(face => mesh % element(e) % face(1:6))
        where(face(:)%boundary > 0)
          bc_face = bc(face%boundary)
        elsewhere
          bc_face = ''
        end where
      end associate

      cfg(1,e) = ConfigurationID(bc_face(1:2))
      cfg(2,e) = ConfigurationID(bc_face(3:4))
      cfg(3,e) = ConfigurationID(bc_face(5:6))

      g0 = dx(1) * dx(2) * dx(3)
      g1 = dx(2) * dx(3) / dx(1)
      g2 = dx(3) * dx(1) / dx(2)
      g3 = dx(1) * dx(2) / dx(3)

      ! mean diffusivity: ν₀ = (Δξ Δη Δζ)⁻¹ ∫∫∫ ν dξ dη dζ
      do k = 1, np
        do j = 1, np
          nu_23(j,k) = dot_product(eop%w, nu(:,j,k,e))
        end do
        nu_3(k) = dot_product(eop%w, nu_23(:,k))
      end do
      nu_0 = dot_product(eop%w, nu_3) / 8

      do k = 1, n3
      do j = 1, n2
      do i = 1, n1

        D_inv(i,j,k,e) = ONE / ( lambda * g0                    &
                               + nu_0 * ( g1 * V1(i, cfg(1,e))  &
                                        + g2 * V2(j, cfg(2,e))  &
                                        + g3 * V3(k, cfg(3,e))  &
                                        )                       &
                               )
      end do
      end do
      end do

    end do

  end associate

end subroutine BuildSubdomains_VI

!===============================================================================

end submodule MP_BuildSubdomains_VI
