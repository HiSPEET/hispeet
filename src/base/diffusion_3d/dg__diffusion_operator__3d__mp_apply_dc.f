!> summary:  3D DG diffusion operator: application with deformed mesh and
!>           constant diffusivity
!> author:   Joerg Stiller
!> date:     2021/10/6
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(DG__Diffusion_Operator__3D:MP_Apply) MP_Apply_DC
  use TPO__Diffusion__3D
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

      ! transfer traces ........................................................

      call tr_buf % Transfer(mesh, tr, tag=1000)
      call tr_buf % Merge(tr)

      call ApplyBoundaryConditions( mesh, this%bc          &
                                  , tr_u  = tr(:,:,:,:,1)  &
                                  , tr_qn = tr(:,:,:,:,2)  )

      ! add fluxes .............................................................

      call AddFluxes( mesh, eop, nu         &
                    , tr_u  = tr(:,:,:,:,1) &
                    , tr_qn = tr(:,:,:,:,2) &
                    , f     = f             &
                    , v     = v             )

      ! clean-up ...............................................................

      !$omp master
      deallocate(tr, tr_buf)
      !$omp end master

    end associate

  end subroutine Apply_DC

  !-----------------------------------------------------------------------------
  !> Compute & add fluxes through element boundaries and, optionally, apply RHS

  subroutine AddFluxes(mesh, metrics, eop, nu, tr_u, tr_qn, f, v)

    ! arguments ................................................................

    class(Mesh_3D), intent(in) :: mesh !< mesh partition
    class(MeshMetrics_3D), intent(in) :: metrics !< metric coefficients
    class(DG_ElementOperators_1D), intent(in) :: eop  !< ID-DG element operators

    real(RNP), intent(in) :: nu               !< diffusivity
    real(RNP), intent(in) :: tr_u (0:,0:,:,:) !< u nᵢ @ element faces
    real(RNP), intent(in) :: tr_qn(0:,0:,:,:) !< q_n  @ element faces
    real(RNP), optional, intent(in) :: f(0:,0:,0:,:) !< RHS
    real(RNP), intent(inout) :: v(0:,0:,0:,:) !< result

    real(RNP), dimension(0:eop%po, 0:eop%po, 3) :: B
    real(RNP), dimension(0:eop%po, 0:eop%po)    :: Mf, Aq, Ju
    real(RNP) :: tmp
    integer   :: en(6), fn(6)
    integer   :: e, f, i, j, k, l

    associate( P    => eop % po,      &
               Ms   => eop % w,       &
               Ds   => eop % D,       &
               a    => metrics % a,   &
               Ji_n => metrics % Ji_n )

      ! face standard mass matrix
      do j = 0, P
      do i = 0, P
        Mf(i,j) = Ms(i) * Ms(j)
      end do
      end do

      !$omp do private(e)
      do e = 1, mesh % n_elem
        associate(element => mesh % element(e))

          ! identify the neighbor elements and their adjoining faces
          do k = 1, 6
            i = element % face(k) % i_neighbor
            if (i > 0) then
              ! assumes i_neighbor ≥ 0, iff neighbor exists (including ghosts)
              en(k) = element % neighbor(i) % id
              fn(k) = element % neighbor(i) % component
            else
              en(k) = 0
              fn(k) = 0
            end if
          end do

          ! face 1 (west) ......................................................

          ! n·[u] and n·{ν∇u}
          call GetElementBoundaryFluxes( e, 1, en(1), fn(1), P+1 &
                                       , 1, tr_u, tr_qn, Ju, Aq  )


          ! face 2 (east) ......................................................

          f = 2

### dx based on adjoinig cuboids (precomputed)
### mu_nu(f) = eop % PenaltyFactor(dx(f))

          ! Ju = n·[u], Aq = n·{ν∇u}
          call GetElementBoundaryFluxes( e, f, en(f), fn(f), P+1 &
                                       , 1, tr_u, tr_qn, Ju, Aq  )

          do k = 0, P
          do j = 0, P
            B(j,k,1) = HALF * nu * a(j,k,f,e) * Ji_n(j,k,f,e,1) * Ju(j,k)
            B(j,k,2) = HALF * nu * a(j,k,f,e) * Ji_n(j,k,f,e,2) * Ju(j,k)
            B(j,k,3) = HALF * nu * a(j,k,f,e) * Ji_n(j,k,f,e,3) * Ju(j,k)
          end do
          end do

          do k = 0, P
          do j = 0, P

            do i = 0, P
              v(i,j,k,e) = v(i,j,k,e) + Mf(j,k) * Ds(P,i) * B(j,k,1)
            end do

            tmp = a(j,k,f,e) * (Aq(j,k) - mu_nu(f) * Ju(j,k))
            do l = 0, P
              tmp = tmp + Ds(l,j) * B(l,k,2) + Ds(l,k) * B(j,l,3)
            end do
            v(P,j,k,e) = v(P,j,k,e) + Mf(j,k) * tmp

          end do
          end do

        end associate

      end do
    end associate

  end subroutine AddFluxes

  !=============================================================================

end submodule MP_Apply_DC
