!> summary:  Base type of SDC correctors for incompressible flows
!> author:   Joerg Stiller
!> date:     2020/05/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CART__ISP_Flow__SDC_Corrector
  use Kind_Parameters, only: RNP
  use ISP_Flow_Problem
  use CART__ISP_Flow__Operators
  implicit none
  private

  public :: SDC_Corrector

  !-----------------------------------------------------------------------------
  !> Abstract type of a SDC corrector for incompressible flow

  type, abstract :: SDC_Corrector
    class(FlowProblem),   pointer :: problem => null() !< flow problem
    class(FlowOperators), pointer :: flow_op => null() !< flow operators
  contains
    procedure :: Init_SDC_Corrector
    procedure(GetCorrectionRHS), deferred :: GetCorrectionRHS
    procedure(CorrectionStep),   deferred :: CorrectionStep
  end type SDC_Corrector

  abstract interface

    !---------------------------------------------------------------------------
    !> Computes F_ex and F_im as defined in the corrector

    subroutine GetCorrectionRHS(this, t, nu, u_ex, u_im, F_ex, F_im,)
      import

      class(SDC_Corrector_EulerDS), intent(in) :: this
      real(RNP), intent(in)  :: t                !< time
      real(RNP), intent(in)  :: nu   (:,:,:,:,:) !< diffusivity
      real(RNP), intent(in)  :: u_ex (:,:,:,:,:) !< u applied in F_ex
      real(RNP), intent(in)  :: u_im (:,:,:,:,:) !< u applied in F_im and F
      real(RNP), intent(out) :: F_ex (:,:,:,:,:) !< explicit RHS for corrector
      real(RNP), intent(out) :: F_im (:,:,:,:,:) !< implicit RHS for corrector

      optional :: nu

    end subroutine GetCorrectionRHS

    !---------------------------------------------------------------------------
    !> Execution of a single correction step

    subroutine CorrectionStep( this, t, dt               &
                             , F_ex_0_old, F_ex_0, F_ex  &
                             , F_im_0_old, F_im_0, F_im  &
                             , S, u_0, u, nu             )
      import

      class(SDC_Corrector_EulerDS), intent(inout) :: this
      real(RNP), intent(inout) :: t                      !< time t₀ → t
      real(RNP), intent(in)    :: dt                     !< step size ∆t = t-t₀
      real(RNP), intent(in)    :: F_ex_0_old (:,:,:,:,:) !< F^ex (t₀)ᵏ⁻¹
      real(RNP), intent(in)    :: F_ex_0     (:,:,:,:,:) !< F^ex (t₀)ᵏ
      real(RNP), intent(inout) :: F_ex       (:,:,:,:,:) !< F^ex (t )ᵏ⁻¹
      real(RNP), intent(in)    :: F_im_0_old (:,:,:,:,:) !< F^im (t₀)ᵏ⁻¹
      real(RNP), intent(in)    :: F_im_0     (:,:,:,:,:) !< F^im (t₀)ᵏ
      real(RNP), intent(inout) :: F_im       (:,:,:,:,:) !< F^im (t )ᵏ⁻¹
      real(RNP), intent(in)    :: S          (:,:,:,:,:) !< S    (t₀)ᵏ⁻¹
      real(RNP), intent(in)    :: u_0        (:,:,:,:,:) !< u    (t₀)ᵏ
      real(RNP), intent(inout) :: u          (:,:,:,:,:) !< u    (t )ᵏ⁻¹ → (t)ᵏ
      real(RNP), intent(inout) :: nu         (:,:,:,:,:) !< v    (t )ᵏ⁻¹

      optional :: nu

    end subroutine CorrectionStep

  end interface

contains

  !-----------------------------------------------------------------------------
  !> Initialization of SDC_Corrector object

  subroutine SDC_Corrector(this, problem, flow_op)
    class(SDC_Corrector),         intent(inout) :: this
    class(FlowProblem),   target, intent(in)    :: problem !< flow problem
    class(FlowOperators), target, intent(in)    :: flow_op !< flow operators

    this % problem => problem
    this % flow_op => flow_op

  end subroutine Init_SDC_Corrector

  !=============================================================================

end module CART__ISP_Flow__SDC_Corrector
