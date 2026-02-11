!> summary:  3D DG elliptic operator: evaluation with deformed mesh and
!>           variable isotropic diffusivity
!> author:   Joerg Stiller
!> date:     2024/08/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(DG__Elliptic_Operator__3D) MP_Eval_DV
  use TPO__Elliptic__3D
  use Mesh_Metrics__3D
  use Mesh__3D
  use Element_Face_Transfer_Buffer__3D
  implicit none

  !=============================================================================
  ! Interfaces to separate module procedures

  interface

    !---------------------------------------------------------------------------
    !> Compute & add fluxes: generic version

    module subroutine AddFluxes_DV_Gen(this, hom_bc, a, Ji_n, tr, r, f)
      class(DG_EllipticOperator_3D),   intent(in)    :: this
      logical,                         intent(in)    :: hom_bc
      real(RNP), contiguous,           intent(in)    :: a(0:,0:,:,:)
      real(RNP), contiguous,           intent(in)    :: Ji_n(0:,0:,:,:,:)
      real(RNP), contiguous,           intent(in)    :: tr(0:,0:,:,:,:)
      real(RNP), contiguous,           intent(inout) :: r(0:,0:,0:,:)
      real(RNP), contiguous, optional, intent(in)    :: f(0:,0:,0:,:)
    end subroutine AddFluxes_DV_Gen

    !---------------------------------------------------------------------------
    !> Compute & add fluxes: LIBXSMM version

    module subroutine AddFluxes_DV_XSMM(this, hom_bc, a, Ji_n, tr, r, f)
      class(DG_EllipticOperator_3D),   intent(in)    :: this
      logical,                         intent(in)    :: hom_bc
      real(RNP), contiguous,           intent(in)    :: a(0:,0:,:,:)
      real(RNP), contiguous,           intent(in)    :: Ji_n(0:,0:,:,:,:)
      real(RNP), contiguous,           intent(in)    :: tr(0:,0:,:,:,:)
      real(RNP), contiguous,           intent(inout) :: r(0:,0:,0:,:)
      real(RNP), contiguous, optional, intent(in)    :: f(0:,0:,0:,:)
    end subroutine AddFluxes_DV_XSMM

  end interface

contains

  !-----------------------------------------------------------------------------
  !> Evaluation with deformed mesh and variable isotropic diffusivity
  !>
  !> Usage
  !>   1) `f` and `bv` given:  computation of the residual, `r = f - Au`
  !>   2) `f` and `bv` absent: evaluation of the homogeneous operator, `r = Au`

  module subroutine Eval_DV(this, bc, lambda, nu, u, r, f, bv)
    class(DG_EllipticOperator_3D),   intent(in)  :: this
    character,                       intent(in)  :: bc(:)       !< {P,D,N}
    real(RNP),                       intent(in)  :: lambda      !< λ
    real(RNP), contiguous,           intent(in)  :: nu(:,:,:,:) !< ν
    real(RNP), contiguous,           intent(in)  :: u(:,:,:,:)  !< operand
    real(RNP), contiguous,           intent(out) :: r(:,:,:,:)  !< result
    real(RNP), contiguous, optional, intent(in)  :: f(:,:,:,:)  !< RHS
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:)
    !< boundary values

    ! local variables ..........................................................

    ! trace transfer buffers operators
    type(ElementFaceTransferBuffer_3D), asynchronous, allocatable, save :: buf_tr

    ! trace variables
    real(RNP), allocatable, save :: tr(:,:,:,:,:) ! traces of ν, u and q_n

    logical :: hom_bc
    integer :: po, ne, ng, np

    if (this % sem % mesh % part < 0) return

    associate( eop  => this % eop                  &
             , sem  => this % sem                  &
             , mesh => this % sem % mesh           &
             , a    => this % sem % metrics % a    &
             , Ji_n => this % sem % metrics % Ji_n )

      ! initialization .........................................................

      hom_bc = .not. present(bv)

      po = eop  % po
      ne = mesh % n_elem
      ng = mesh % n_ghost
      np = po + 1

      ! workspace and operators
      !$omp master
      allocate(tr(0:po, 0:po, 6, ne+ng,3))
      buf_tr = ElementFaceTransferBuffer_3D(mesh, tr)
      !$omp end master
      !$omp barrier

      ! apply element diffusion operator .......................................

      call TPO_Elliptic( eop, sem            & ! As r may contain ghost entries
                       , lambda, nu, u       & ! it is explicitly restricted to
                       , v   = r(:,:,:,:ne)  & ! local elements since no proper
                       , nub = tr(:,:,:,:,1) & ! bound checking is performed in
                       , ub  = tr(:,:,:,:,2) & ! bound checking is performed in
                       , qb  = tr(:,:,:,:,3) ) ! TPO_Elliptic.

      ! transfer traces and apply boundary conditions ..........................

      call buf_tr % Transfer(mesh, tr, tag=1000)

      call EnforceBoundaryConditions( this, bc, bv          &
                                    , jmp_u = tr(:,:,:,:,2) &
                                    , avg_q = tr(:,:,:,:,3) )

      call buf_tr % Merge(tr)

      ! add fluxes .............................................................

      select case(po)
      case(1:2)
        call AddFluxes_DV_Gen  (this, hom_bc, a, Ji_n, tr, r, f)
      case default
#ifdef __LIBXSMM__
        call AddFluxes_DV_XSMM (this, hom_bc, a, Ji_n, tr, r, f)
#else
        call AddFluxes_DV_Gen  (this, hom_bc, a, Ji_n, tr, r, f)
#endif
      end select

      ! frozen and ghost elements are set to zero ..............................

      call SetArray(r(:,:,:,mesh%n_elem_active+1:), ZERO)

      ! clean-up ...............................................................

      !$omp master
      deallocate(tr, buf_tr)
      !$omp end master

    end associate

  end subroutine Eval_DV

  !=============================================================================

end submodule MP_Eval_DV
