!> summary:  3D DG diffusion operator: TBP for setting diffusivity
!> author:   Joerg Stiller
!> date:     2021/08/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(DG__Diffusion_Operator__3D) MP_SetDiffusivity
  use Spectral_Element_Scalar__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Set constant physical and spectral diffusivities

  module subroutine SetDiffusivity_C(this, nu_p, nu_s)
    class(DG_DiffusionOperator_3D), intent(inout) :: this
    real(RNP),           intent(in) :: nu_p   !< physical diffusivity
    real(RNP), optional, intent(in) :: nu_s   !< spectral diffusivity [0]

    this % nu_pc  = nu_p

    if (present(nu_s)) then
      this % nu_sc = nu_s
    else
      this % nu_sc = 0
    end if

    if (allocated(this % nu_pv)) deallocate(this % nu_pv)
    if (allocated(this % nu_mx)) deallocate(this % nu_mx)

  end subroutine SetDiffusivity_C

  !-----------------------------------------------------------------------------
  !> Set variable physical diffusivity

  module subroutine SetDiffusivity_V(this, nu_p)
    class(DG_DiffusionOperator_3D), intent(inout) :: this
    real(RNP), intent(in) :: nu_p(:,:,:,:)  !< physical diffusivity

    type(SpectralElementScalar_3D), allocatable, save :: se_nu
    real(RNP), allocatable, save :: tr_nu(:,:,:,:)
    integer :: f, ne, nf, po

    this % nu_pc  = 0
    this % nu_sc  = 0

    po = this % sem % std_op % po
    ne = this % sem % mesh % n_elem
    nf = this % sem % mesh % n_face

    !$omp master

    if (allocated(this % nu_pv)) then
      if (any(shape(this % nu_pv) /= shape(nu_p))) deallocate(this % nu_pv)
    end if
    if (.not. allocated(this % nu_pv)) then
      allocate(this % nu_pv(0:po,0:po,0:po,ne), source = nu_p)
    else
      this % nu_pv = nu_p
    end if

    if (allocated(this % nu_mx)) then
      if (any(shape(this % nu_mx) /= [po+1,po+1,nf])) deallocate(this % nu_mx)
    end if
    if (.not. allocated(this % nu_mx)) then
      allocate(this % nu_mx(0:po,0:po,nf))
    end if

    allocate(tr_nu(0:po,0:po,2,nf))
    se_nu = SpectralElementScalar_3D(this % sem, this % nu_pv)

    !$omp end master

    ! apply tracing procedure of scalar SE variable
    call se_nu % GetTrace(tr_nu)

    ! determine maximum face diffusivitiy
    !$omp do
    do f = 1, nf
      this % nu_mx(:,:,f) = max(tr_nu(:,:,1,f), tr_nu(:,:,2,f))
    end do

    !$omp master
    deallocate(se_nu, tr_nu)
    !$omp end master

  end subroutine SetDiffusivity_V

  !=============================================================================

end submodule MP_SetDiffusivity
