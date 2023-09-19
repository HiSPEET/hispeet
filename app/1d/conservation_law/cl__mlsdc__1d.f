!> summary:  MultiLevel SDC method for 1D conservation laws
!> author:   Joerg Stiller
!> date:     2023/09/17
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CL__MLSDC__1D
  use Kind_Parameters
  use Standard_Operators__1D
  use DG__Element_Operators__1D
  use DG__Schwarz_Operator__1D
  use Coarse_To_Fine_Interpolation__1D
  use Fine_To_Coarse_Projection__1D
  use CL__Problem__1D
  use CL__Operator__1D
  use CL__Time_Integrator__1D
  use CL__SDC__Method__1D
  use CL__MLSDC__Level__1D

  implicit none
  private

  !-----------------------------------------------------------------------------
  !> MLSDC for 1D conservation laws

  type, public :: CL_MLSDC_1D
    type(CL_MLSDC_Level_1D), allocatable :: level(:)
  end type CL_MLSDC_1D

  ! constructor
  interface CL_MLSDC_1D
    procedure New_CL_MLSDC_1D
  end interface

  !-----------------------------------------------------------------------------
  !> Options for initializing CL_MLSDC_1D objects

  type, public :: CL_MLSDC_Options_1D(n_level)

    integer, len :: n_level = 1  !< number of space-time levels

    ! level parameters .........................................................

    integer :: n_elem (n_level)  !< number of elements in space
    integer :: p_elem (n_level)  !< polynomial degree of elements space
    integer :: q_conv (n_level)  !< quadrature degree for convection
    integer :: n_step (n_level)  !< number of time steps in one slice
    integer :: n_sub  (n_level)  !< number of subintervals in time step

    ! space options ............................................................

    real(RNP) :: penalty = 2  !< penalty parameter > 1

    type(DG_SchwarzOptions_1D) :: schwarz_root !< Schwarz opts for root level
    type(DG_SchwarzOptions_1D) :: schwarz_fine !< Schwarz opts for finer levels

    ! transfer options .........................................................

    character :: projection_method = 'P' !< fine-to-coarse projection method:
                                         !! 'P'  L² projection,
                                         !! 'I'  interpolation

    integer :: projection_smoothing = 1 !< discontinuities in 2:1 projection:
                                        !!  0   no smoothing
                                        !!  1   removal by linear blending
                                        !!  0   removal by coefficient averaging

  end type CL_MLSDC_Options_1D

contains

  !-----------------------------------------------------------------------------
  !> New CL_MLSDC_1D object

  function New_CL_MLSDC_1D(opt, opt_pre, opt_sdc, cl_problem) result(this)
    class(CL_MLSDC_Options_1D),          intent(in) :: opt
    class(CL_TimeIntegrator_Options_1D), intent(in) :: opt_pre
    class(CL_SDC_Options_1D),            intent(in) :: opt_sdc
    class(CL_Problem_1D),                intent(in) :: cl_problem
    type(CL_MLSDC_1D) :: this

    call Init_CL_MLSDC_1D(this, opt, opt_pre, opt_sdc, cl_problem)

  end function New_CL_MLSDC_1D

  !-----------------------------------------------------------------------------
  !> Initialization of a CL_MLSDC_1D object

  subroutine Init_CL_MLSDC_1D(this, opt, opt_pre, opt_sdc, cl_problem)
    class(CL_MLSDC_1D),                  intent(inout) :: this
    class(CL_MLSDC_Options_1D),          intent(in)    :: opt
    class(CL_TimeIntegrator_Options_1D), intent(in)    :: opt_pre
    class(CL_SDC_Options_1D),            intent(in)    :: opt_sdc
    class(CL_Problem_1D),                intent(in)    :: cl_problem

    type(CL_MLSDC_Level_Options_1D) :: opt_level
    integer :: l

    ! preliminaries ............................................................

    if (allocated(this % level)) then
      deallocate(this % level)
    end if

    allocate(this % level(opt % n_level))

    ! consistency check ........................................................

    ! TBD

    ! generate levels ..........................................................

    do l = 1, opt % n_level

      ! space options  . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

      associate(cl_operator => opt_level % cl_operator)

        cl_operator % ne  = opt % n_elem(l)
        cl_operator % nc  = cl_problem % nc
        cl_operator % xb1 = cl_problem % xb1
        cl_operator % xb2 = cl_problem % xb2

        cl_operator % eop = &
            DG_ElementOptions_1D(po = opt % p_elem(l), penalty = opt % penalty)

        cl_operator % qop = StandardOperatorOptions_1D(po = opt % q_conv(l))

        select case(l)
        case(1)
          cl_operator % schwarz = opt % schwarz_root
        case default
          cl_operator % schwarz = opt % schwarz_fine
        end select

      end associate

      ! time options . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

      opt_level % n_step = opt % n_step(l)

      opt_level % cl_pre = opt_pre
      opt_level % cl_sdc = opt_sdc
      opt_level % cl_sdc % n_sub = opt % n_sub(l)

      ! transfer options . . . . . . . . . . . . . . . . . . . . . . . . . . . .

      if (l < opt % n_level) then

        opt_level % iop_cf_x =                              &
            CoarseToFineInterpolationOptions_1D(            &
                po_c = opt % p_elem(l),                     &
                po_f = opt % p_elem(l+1),                   &
                mode = opt % n_elem(l+1) / opt % n_elem(l)  )

        opt_level % iop_cf_t =                              &
            CoarseToFineInterpolationOptions_1D(            &
                po_c = opt % n_sub(l),                      &
                po_f = opt % n_sub(l+1),                    &
                mode = opt % n_step(l+1) / opt % n_step(l)  )

      end if

      if (l > 1) then

        opt_level % pop_fc_x =                                    &
            FineToCoarseProjectionOptions_1D(                     &
                po_f      = opt % p_elem(l),                      &
                po_c      = opt % p_elem(l-1),                    &
                mode      = opt % n_elem(l) / opt % n_elem(l-1),  &
                method    = opt % projection_method,              &
                smoothing = opt % projection_smoothing            )

        opt_level % pop_fc_t =                                    &
            FineToCoarseProjectionOptions_1D(                     &
                po_f      = opt % n_sub(l),                       &
                po_c      = opt % n_sub(l-1),                     &
                mode      = opt % n_step(l) / opt % n_step(l-1),  &
                method    = opt % projection_method               )

      end if

      ! create level  . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

      this % level(l) = CL_MLSDC_Level_1D(opt_level, cl_problem)

    end do

  end subroutine Init_CL_MLSDC_1D

  !=============================================================================

end module CL__MLSDC__1D
