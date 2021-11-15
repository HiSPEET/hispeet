!> summary:  3D DG diffusion operator: application with deformed mesh and
!>           constant diffusivity
!> author:   Joerg Stiller
!> date:     2021/10/6
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(DG__Diffusion_Operator__3D:MP_Apply) MP_Apply_DC
  use TPO__Diffusion__3D
  use Mesh_Metrics__3D
  use Element_Face_Transfer_Buffer__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Application with irregular (deformed) mesh and constant diffusivity

  module subroutine Apply_DC(this, u, v, f)
    class(DG_DiffusionOperator_3D), intent(in) :: this
    real(RNP), intent(in)  :: u(:,:,:,:) !< operand
    real(RNP), intent(out) :: v(:,:,:,:) !< result
    real(RNP), intent(in), optional :: f(:,:,:,:) !< RHS

    ! local variables ..........................................................

    ! trace transfer buffers operators
    type(ElementFaceTransferBuffer_3D), asynchronous, allocatable, save :: tr_buf

    ! trace variables
    real(RNP), allocatable, save :: tr(:,:,:,:,:) ! traces of u and q_n

    integer :: po, ne, ng, np

    associate( mesh   => this % sem % mesh           &
             , Jd     => this % sem % metrics % Jd   &
             , G      => this % sem % metrics % G    &
             , a      => this % sem % metrics % a    &
             , Ji_n   => this % sem % metrics % Ji_n &
             , lambda => this % lambda               &
             , nu     => this % nu_pc                &
             , eop    => this % eop                  )

      ! initialization .........................................................

      po = eop  % po
      ne = mesh % n_elem
      ng = mesh % n_ghost
      np = po + 1

      ! workspace and operators
      !$omp master
      allocate(tr(0:po, 0:po, 6, ne+ng,2))
      tr_buf = ElementFaceTransferBuffer_3D(mesh, tr)
      !$omp end master
      !$omp barrier

      ! apply element diffusion operator .......................................

      call TPO_Diffusion( eop%w, eop%D, Jd, G, lambda, nu, u, v &
                        , Ji_n                                  &
                        , ub = tr(:,:,:,:,1)                    &
                        , qb = tr(:,:,:,:,2)                    )

      ! transfer traces and apply boundary conditions ..........................

      call tr_buf % Transfer(mesh, tr, tag=1000)

      call ApplyBoundaryConditions( mesh, this%bc          &
                                  , tr_u  = tr(:,:,:,:,1)  &
                                  , tr_qn = tr(:,:,:,:,2)  )
      call tr_buf % Merge(tr)

      ! add fluxes .............................................................

      call AddFluxes( mesh, eop, a, Ji_n, nu &
                    , tr_u  = tr(:,:,:,:,1)  &
                    , tr_qn = tr(:,:,:,:,2)  &
                    , f     = f              &
                    , v     = v              )

      ! clean-up ...............................................................

      !$omp master
      deallocate(tr, tr_buf)
      !$omp end master

    end associate

  end subroutine Apply_DC

  !-----------------------------------------------------------------------------
  !> Compute & add fluxes through element boundaries and, optionally, apply RHS

  subroutine AddFluxes(mesh, eop, a, Ji_n, nu, tr_u, tr_qn, f, v)

    ! arguments ................................................................

    class(Mesh_3D), intent(in) :: mesh !< mesh partition
    class(DG_ElementOperators_1D), intent(in) :: eop !< ID-DG element operators

    real(RNP), intent(in)    :: a(0:,0:,:,:)      !< area coeff @ element faces
    real(RNP), intent(in)    :: Ji_n(0:,0:,:,:,:) !< J⁻¹⋅n      @ element faces
    real(RNP), intent(in)    :: nu                !< diffusivity
    real(RNP), intent(in)    :: tr_u (0:,0:,:,:)  !< u nᵢ       @ element faces
    real(RNP), intent(in)    :: tr_qn(0:,0:,:,:)  !< q_n        @ element faces
    real(RNP), intent(inout) :: v(0:,0:,0:,:)     !< result

    real(RNP), optional, intent(in) :: f(0:,0:,0:,:) !< RHS

    ! local variables ..........................................................

    real(RNP), dimension(0:eop%po, 0:eop%po, 3) :: Cu
    real(RNP), dimension(0:eop%po, 0:eop%po)    :: Mf, Aq, Ju
    real(RNP) :: mu_nu, tmp
    integer   :: e, i, j, k, m, n
    logical   :: present_f, struct

    associate( P  => eop % po, &
               Ms => eop % w,  &
               Ds => eop % D   )

      ! auxiliaries ............................................................

      ! face standard mass matrix
      do j = 0, P
      do i = 0, P
        Mf(i,j) = Ms(i) * Ms(j)
      end do
      end do

      present_f = present(f)
      struct    = mesh % structured

      ! add fluxes .............................................................

      !$omp do private(e)
      do e = 1, mesh % n_elem
        associate(element => mesh % element(e))

          ! faces 1+2 (west + east)  . . . . . . . . . . . . . . . . . . . . . .

          do m = 1, 2

            ! Ju = n·[u], Aq = n·{ν∇u}
            call GetElementBoundaryFluxes(element,struct,e,m,tr_u,tr_qn,Ju,Aq)

            mu_nu = eop % PenaltyFactor(mesh % dx_mean(m,e)) * nu

            i = (m-1) * P

            do k = 0, P
            do j = 0, P
              tmp = Mf(j,k) * HALF * nu * a(j,k,m,e) * Ju(j,k)
              Cu(j,k,1) = tmp * Ji_n(j,k,m,e,1)
              Cu(j,k,2) = tmp * Ji_n(j,k,m,e,2)
              Cu(j,k,3) = tmp * Ji_n(j,k,m,e,3)
            end do
            end do

            do k = 0, P
            do j = 0, P

              do n = 0, P
                v(n,j,k,e) = v(n,j,k,e) - Ds(i,n) * Cu(j,k,1)
              end do

              tmp = Mf(j,k) * a(j,k,m,e) * (Aq(j,k) - mu_nu * Ju(j,k))
              do n = 0, P
                tmp = tmp + Ds(n,j) * Cu(n,k,2) + Ds(n,k) * Cu(j,n,3)
              end do
              v(i,j,k,e) = v(i,j,k,e) - tmp

            end do
            end do

          end do

          ! faces 3+4 (south + north)  . . . . . . . . . . . . . . . . . . . . .

          do m = 3, 4

            call GetElementBoundaryFluxes(element,struct,e,m,tr_u,tr_qn,Ju,Aq)
            mu_nu = eop % PenaltyFactor(mesh % dx_mean(m,e)) * nu
            j = (m-3) * P

            do k = 0, P
            do i = 0, P
              tmp = Mf(i,k) * HALF * nu * a(i,k,m,e) * Ju(i,k)
              Cu(i,k,1) = tmp * Ji_n(i,k,m,e,1)
              Cu(i,k,2) = tmp * Ji_n(i,k,m,e,2)
              Cu(i,k,3) = tmp * Ji_n(i,k,m,e,3)
            end do
            end do

            do k = 0, P
            do i = 0, P

              do n = 0, P
                v(i,n,k,e) = v(i,n,k,e) - Ds(j,n) * Cu(i,k,2)
              end do

              tmp = Mf(i,k) * a(i,k,m,e) * (Aq(i,k) - mu_nu * Ju(i,k))
              do n = 0, P
                tmp = tmp + Ds(n,i) * Cu(n,k,1) + Ds(n,k) * Cu(i,n,3)
              end do
              v(i,j,k,e) = v(i,j,k,e) - tmp

            end do
            end do

          end do

          ! faces 5+6 (bottom + top) . . . . . . . . . . . . . . . . . . . . . .

          do m = 5, 6

            call GetElementBoundaryFluxes(element,struct,e,m,tr_u,tr_qn,Ju,Aq)
            mu_nu = eop % PenaltyFactor(mesh % dx_mean(m,e)) * nu
            k = (m-5) * P

            do j = 0, P
            do i = 0, P
              tmp = Mf(i,j) * HALF * nu * a(i,j,m,e) * Ju(i,j)
              Cu(i,j,1) = tmp * Ji_n(i,j,m,e,1)
              Cu(i,j,2) = tmp * Ji_n(i,j,m,e,2)
              Cu(i,j,3) = tmp * Ji_n(i,j,m,e,3)
            end do
            end do

            do j = 0, P
            do i = 0, P

              do n = 0, P
                v(i,j,n,e) = v(i,j,n,e) - Ds(k,n) * Cu(i,j,3)
              end do

              tmp = Mf(i,j) * a(i,j,m,e) * (Aq(i,j) - mu_nu * Ju(i,j))
              do n = 0, P
                tmp = tmp + Ds(n,i) * Cu(n,j,1) + Ds(n,j) * Cu(i,n,2)
              end do
              v(i,j,k,e) = v(i,j,k,e) - tmp

            end do
            end do

          end do

          ! RHS  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

          if (present_f) then
            do k = 0, P
            do j = 0, P
            do i = 0, P
              v(i,j,k,e) = v(i,j,k,e) - f(i,j,k,e)
            end do
            end do
            end do
          end if

        end associate

      end do
    end associate

  end subroutine AddFluxes

  !=============================================================================

end submodule MP_Apply_DC
