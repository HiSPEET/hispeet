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


  !=============================================================================
  ! Interfaces to separate module procedures

  interface

    !---------------------------------------------------------------------------
    !> Compute & add fluxes: generic version

    module subroutine AddFluxes_DC_Gen(mesh, eop, a, Ji_n, nu, tr_u, tr_qn, f, v)
      class(Mesh_3D),                  intent(in)    :: mesh
      class(DG_ElementOperators_1D),   intent(in)    :: eop
      real(RNP), contiguous,           intent(in)    :: a(0:,0:,:,:)
      real(RNP), contiguous,           intent(in)    :: Ji_n(0:,0:,:,:,:)
      real(RNP),                       intent(in)    :: nu
      real(RNP), contiguous,           intent(in)    :: tr_u (0:,0:,:,:)
      real(RNP), contiguous,           intent(in)    :: tr_qn(0:,0:,:,:)
      real(RNP), contiguous,           intent(inout) :: v(0:,0:,0:,:)
      real(RNP), contiguous, optional, intent(in)    :: f(0:,0:,0:,:)
    end subroutine AddFluxes_DC_Gen

    !---------------------------------------------------------------------------
    !> Compute & add fluxes: LIBXSMM version

    module subroutine AddFluxes_DC_XSMM(mesh, eop, a, Ji_n, nu, tr_u, tr_qn, f, v)
      class(Mesh_3D),                  intent(in)    :: mesh
      class(DG_ElementOperators_1D),   intent(in)    :: eop
      real(RNP), contiguous,           intent(in)    :: a(0:,0:,:,:)
      real(RNP), contiguous,           intent(in)    :: Ji_n(0:,0:,:,:,:)
      real(RNP),                       intent(in)    :: nu
      real(RNP), contiguous,           intent(in)    :: tr_u (0:,0:,:,:)
      real(RNP), contiguous,           intent(in)    :: tr_qn(0:,0:,:,:)
      real(RNP), contiguous,           intent(inout) :: v(0:,0:,0:,:)
      real(RNP), contiguous, optional, intent(in)    :: f(0:,0:,0:,:)
    end subroutine AddFluxes_DC_XSMM

  end interface !!

contains

  !-----------------------------------------------------------------------------
  !> Application with irregular (deformed) mesh and constant diffusivity

  module subroutine Apply_DC(this, u, v, f)
    class(DG_DiffusionOperator_3D), intent(in) :: this
    real(RNP), contiguous, intent(in)  :: u(:,:,:,:) !< operand
    real(RNP), contiguous, intent(out) :: v(:,:,:,:) !< result
    real(RNP), contiguous, intent(in), optional :: f(:,:,:,:) !< RHS

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

      select case(po)
      case(1:2)
        call AddFluxes_DC_Gen ( mesh, eop, a, Ji_n, nu &
                              , tr_u  = tr(:,:,:,:,1)  &
                              , tr_qn = tr(:,:,:,:,2)  &
                              , f     = f              &
                              , v     = v              )
      case default

#ifdef __LIBXSMM__
        call AddFluxes_DC_XSMM( mesh, eop, a, Ji_n, nu &
                              , tr_u  = tr(:,:,:,:,1)  &
                              , tr_qn = tr(:,:,:,:,2)  &
                              , f     = f              &
                              , v     = v              )
#else
         call AddFluxes_DC_Gen ( mesh, eop, a, Ji_n, nu &
                               , tr_u  = tr(:,:,:,:,1)  &
                               , tr_qn = tr(:,:,:,:,2)  &
                               , f     = f              &
                               , v     = v              )
#endif

      end select

      ! clean-up ...............................................................

      !$omp master
      deallocate(tr, tr_buf)
      !$omp end master

    end associate

  end subroutine Apply_DC

  !=============================================================================

end submodule MP_Apply_DC
