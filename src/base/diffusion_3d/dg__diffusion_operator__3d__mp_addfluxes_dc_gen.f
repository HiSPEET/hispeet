!> summary:  3D DG diffusion operator: application of fluxes with deformed mesh
!>           and constant diffusivity, generic version
!> author:   Joerg Stiller
!> date:     2021/11/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(DG__Diffusion_Operator__3D:MP_Apply_DC) MP_AddFluxes_DC_Gen
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Compute & add fluxes: generic version

  module subroutine AddFluxes_DC_Gen(mesh, eop, a, Ji_n, nu, tr, r, f)

    ! arguments ................................................................

    class(Mesh_3D),                intent(in) :: mesh !< mesh partition
    class(DG_ElementOperators_1D), intent(in) :: eop  !< ID-DG element operators

    real(RNP), intent(in)    :: a(0:,0:,:,:)      !< area coeff @ element faces
    real(RNP), intent(in)    :: Ji_n(0:,0:,:,:,:) !< J⁻¹⋅n      @ element faces
    real(RNP), intent(in)    :: nu                !< diffusivity
    real(RNP), intent(in)    :: tr(0:,0:,:,:,:)   !< traces of u, q_n
    real(RNP), intent(inout) :: r(0:,0:,0:,:)     !< result
    real(RNP), intent(in), optional :: f(0:,0:,0:,:) !< RHS

    contiguous :: a, Ji_n, tr, r, f

    ! local variables ..........................................................

    real(RNP), dimension(0:eop%po, 0:eop%po)    :: Mf, Ds_t, jmp_u, avg_q
    real(RNP), dimension(0:eop%po, 0:eop%po, 2) :: Ct, Ds_Ct
    real(RNP), dimension(0:eop%po, 0:eop%po, 6) :: Cn
    real(RNP) :: mu_nu, tmp
    integer   :: e, i, j, k, m
    logical   :: present_f, struct

    associate( P  => eop % po &
             , Ms => eop % w  &
             , Ds => eop % D  )

      ! auxiliaries ............................................................

      ! face standard mass matrix
      do j = 0, P
      do i = 0, P
        Mf(i,j) = Ms(i) * Ms(j)
      end do
      end do

      Ds_t = transpose(eop % D)

      present_f = present(f)
      struct    = mesh % structured

      ! add fluxes .............................................................

      !$omp do private(e)
      do e = 1, mesh % n_elem
        associate(element => mesh % element(e))

          ! faces 1+2 (west + east)  . . . . . . . . . . . . . . . . . . . . . .

          do m = 1, 2

            ! jmp_u = n·[u], avg_q = n·{ν∇u}
            call GetBoundaryFluxes(element, struct, e, m, tr, jmp_u, avg_q)
            mu_nu = eop % PenaltyFactor(mesh % dx_mean(m,e)) * nu
            i = (m-1) * P

            ! normal and tangential contributions of jumps
            do k = 0, P
            do j = 0, P
              tmp = Mf(j,k) * HALF * nu * a(j,k,m,e) * jmp_u(j,k)
              Cn(j,k,m) = tmp * Ji_n(j,k,m,e,1)
              Ct(j,k,1) = tmp * Ji_n(j,k,m,e,2) ! Ct(:,:,2) is transposed for
              Ct(k,j,2) = tmp * Ji_n(j,k,m,e,3) ! fused evaluation of derivative
            end do
            end do

            ! derivatives of tangential jump contributions
            select case(P)
            case(1) ! allow compiler to optimize inlined code for P = 1
              call TangentialDerivatives(1, Ds, Ct, Ds_Ct)
            case(2) ! dito for P = 2
              call TangentialDerivatives(2, Ds, Ct, Ds_Ct)
            case default
              call TangentialDerivatives(P, Ds, Ct, Ds_Ct)
            end select

            ! add flux contributions, accounting for transposition of Ct(:,:,2)
            do k = 0, P
            do j = 0, P

              r(i,j,k,e) = r(i,j,k,e)                           &
                         - ( Mf(j,k) * a(j,k,m,e)               &
                           * (avg_q(j,k) - mu_nu * jmp_u(j,k))  &
                           + Ds_Ct(j,k,1)                       &
                           + Ds_Ct(k,j,2)                       &
                           )
            end do
            end do

          end do

          ! faces 3+4 (south + north)  . . . . . . . . . . . . . . . . . . . . .

          do m = 3, 4

            call GetBoundaryFluxes(element, struct, e, m, tr, jmp_u, avg_q)
            mu_nu = eop % PenaltyFactor(mesh % dx_mean(m,e)) * nu
            j = (m-3) * P

            do k = 0, P
            do i = 0, P
              tmp = Mf(i,k) * HALF * nu * a(i,k,m,e) * jmp_u(i,k)
              Cn(i,k,m) = tmp * Ji_n(i,k,m,e,2)
              Ct(i,k,1) = tmp * Ji_n(i,k,m,e,1)
              Ct(k,i,2) = tmp * Ji_n(i,k,m,e,3)
            end do
            end do

            select case(P)
            case(1)
              call TangentialDerivatives(1, Ds, Ct, Ds_Ct)
            case(2)
              call TangentialDerivatives(2, Ds, Ct, Ds_Ct)
            case default
              call TangentialDerivatives(P, Ds, Ct, Ds_Ct)
            end select

            do k = 0, P
            do i = 0, P

              r(i,j,k,e) = r(i,j,k,e)                           &
                         - ( Mf(i,k) * a(i,k,m,e)               &
                           * (avg_q(i,k) - mu_nu * jmp_u(i,k))  &
                           + Ds_Ct(i,k,1)                       &
                           + Ds_Ct(k,i,2)                       &
                           )
            end do
            end do

          end do

          ! faces 5+6 (bottom + top) . . . . . . . . . . . . . . . . . . . . . .

          do m = 5, 6

            call GetBoundaryFluxes(element, struct, e, m, tr, jmp_u, avg_q)
            mu_nu = eop % PenaltyFactor(mesh % dx_mean(m,e)) * nu
            k = (m-5) * P

            do j = 0, P
            do i = 0, P
              tmp = Mf(i,j) * HALF * nu * a(i,j,m,e) * jmp_u(i,j)
              Cn(i,j,m) = tmp * Ji_n(i,j,m,e,3)
              Ct(i,j,1) = tmp * Ji_n(i,j,m,e,1)
              Ct(j,i,2) = tmp * Ji_n(i,j,m,e,2)
            end do
            end do

            select case(P)
            case(1)
              call TangentialDerivatives(1, Ds, Ct, Ds_Ct)
            case(2)
              call TangentialDerivatives(2, Ds, Ct, Ds_Ct)
            case default
              call TangentialDerivatives(P, Ds, Ct, Ds_Ct)
            end select

            do j = 0, P
            do i = 0, P

              r(i,j,k,e) = r(i,j,k,e)                           &
                         - ( Mf(i,j) * a(i,j,m,e)               &
                           * (avg_q(i,j) - mu_nu * jmp_u(i,j))  &
                           + Ds_Ct(i,j,1)                       &
                           + Ds_Ct(j,i,2)                       &
                           )
            end do
            end do

          end do

          if (present_f) then
            do k = 0, P
            do j = 0, P
            do i = 0, P
              r(i,j,k,e) = r(i,j,k,e)               &
                         - ( Ds_t(i,0) * Cn(j,k,1)  &
                           + Ds_t(i,P) * Cn(j,k,2)  &
                           + Ds_t(j,0) * Cn(i,k,3)  &
                           + Ds_t(j,P) * Cn(i,k,4)  &
                           + Ds_t(k,0) * Cn(i,j,5)  &
                           + Ds_t(k,P) * Cn(i,j,6)  &
                           + f(i,j,k,e)             &
                           )
            end do
            end do
            end do
          else
            do k = 0, P
            do j = 0, P
            do i = 0, P
              r(i,j,k,e) = r(i,j,k,e)               &
                         - ( Ds_t(i,0) * Cn(j,k,1)  &
                           + Ds_t(i,P) * Cn(j,k,2)  &
                           + Ds_t(j,0) * Cn(i,k,3)  &
                           + Ds_t(j,P) * Cn(i,k,4)  &
                           + Ds_t(k,0) * Cn(i,j,5)  &
                           + Ds_t(k,P) * Cn(i,j,6)  &
                           )
            end do
            end do
            end do
          end if

        end associate

      end do
    end associate

  end subroutine AddFluxes_DC_Gen

  !-----------------------------------------------------------------------------
  !> Computation of tangential derivatives

  pure subroutine TangentialDerivatives(P, D, C, DC)
    integer,   intent(in)  :: P              !< polynomial order
    real(RNP), intent(in)  :: D  (0:P,0:P)   !< diff matrix
    real(RNP), intent(in)  :: C  (0:P,0:P,2) !< tangential jump contributions
    real(RNP), intent(out) :: DC (0:P,0:P,2) !< derivatives of C

    integer :: i, j, m

    do j = 0, P
    do i = 0, P
      DC(i,j,1) = 0
      DC(i,j,2) = 0
      do m = 0, P
        DC(i,j,1) = DC(i,j,1) + D(m,i) * C(m,j,1)
        DC(i,j,2) = DC(i,j,2) + D(m,i) * C(m,j,2)
      end do
    end do
    end do

  end subroutine TangentialDerivatives

  !=============================================================================

end submodule MP_AddFluxes_DC_Gen
