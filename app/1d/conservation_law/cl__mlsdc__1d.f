!> summary:  MultiLevel SDC method for 1D conservation laws
!> author:   Joerg Stiller
!> date:     2023/09/17
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CL__MLSDC__1D
  use Kind_Parameters
  use Standard_Element_Operators__1D
  use DG__Element_Operators__1D
  use DG__Schwarz_Operator__1D
  use HP__Refinement_Operator__1D
  use HP__Coarsening_Operator__1D
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

  type, public :: CL_MLSDC_Options_1D

    integer :: n_level = -1  !< number of space-time levels

    ! level parameters
    integer, allocatable :: p_space(:) !< polynomial degree of elements space
    integer, allocatable :: n_space(:) !< number of elements in space
    integer, allocatable :: p_time(:)  !< polynomial degree of time step
    integer, allocatable :: n_time(:)  !< number of time steps in one slice
    integer, allocatable :: q_conv(:)  !< quadrature degree for convection

    ! space options
    real(RNP) :: penalty = 2                   !< IP-DG penalty parameter > 1
    type(DG_SchwarzOptions_1D) :: schwarz_root !< Schwarz opts for root level
    type(DG_SchwarzOptions_1D) :: schwarz_fine !< Schwarz opts for finer levels

    character :: projection_method = 'I' !< fine-to-coarse projection method:
                                         !! 'P'  L² projection,
                                         !! 'I'  interpolation

    integer :: projection_smooth = 1     !< discontinuities in 2:1 projection:
                                         !!  0   no smoothing
                                         !!  1   remove by linear blending
                                         !!  2   remove by coefficient averaging

  end type CL_MLSDC_Options_1D

  ! constructor
  interface CL_MLSDC_Options_1D
    procedure New_CL_MLSDC_Options_1D
  end interface

contains

  !=============================================================================
  ! TBP for CL_MLSDC_1D

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

      ! dimensions . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

      opt_level % p_space = opt % p_space(l)
      opt_level % n_space = opt % n_space(l)
      opt_level % p_time  = opt % p_time (l)
      opt_level % n_time  = opt % n_time (l)

      ! space options  . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

      associate(cl_operator => opt_level % cl_operator)

        cl_operator % ne  = opt % n_space(l)
        cl_operator % xb1 = cl_problem % xb1
        cl_operator % xb2 = cl_problem % xb2

        cl_operator % eop = &
            DG_ElementOptions_1D(po = opt % p_space(l), penalty = opt % penalty)

        cl_operator % qop = StandardElementOptions_1D(po = opt % q_conv(l))

        select case(l)
        case(1)
          cl_operator % schwarz = opt % schwarz_root
        case default
          cl_operator % schwarz = opt % schwarz_fine
        end select

      end associate

      ! time options . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

      opt_level % cl_pre = opt_pre
      opt_level % cl_sdc = opt_sdc

      ! number of collocation points = polynomial degree + 1
      opt_level % cl_sdc % n_col = opt % p_time(l) + 1

      ! transfer options . . . . . . . . . . . . . . . . . . . . . . . . . . . .

      if (l < opt % n_level) then

        opt_level % iop_cf_x =                                &
            HP_RefinementOptions_1D(                          &
                po_c  = opt % p_space(l),                     &
                po_f  = opt % p_space(l+1),                   &
                mode  = opt % n_space(l+1) / opt % n_space(l) )

        opt_level % iop_cf_t =                                &
            HP_RefinementOptions_1D(                          &
                nodes = opt_sdc % nodes,                      &
                po_c  = opt % p_time(l),                      &
                po_f  = opt % p_time(l+1),                    &
                mode  = opt % n_time(l+1) / opt % n_time(l)   )

      else

        ! reset interpolation options
        opt_level % iop_cf_x = HP_RefinementOptions_1D()
        opt_level % iop_cf_t = HP_RefinementOptions_1D()

      end if

      if (l > 1) then

        opt_level % pop_fc_x =                                  &
            HP_CoarseningOptions_1D(                            &
                po_f   = opt % p_space(l),                      &
                po_c   = opt % p_space(l-1),                    &
                mode   = opt % n_space(l) / opt % n_space(l-1), &
                method = opt % projection_method,               &
                smooth = opt % projection_smooth                )

        opt_level % pop_fc_t =                                  &
            HP_CoarseningOptions_1D(                            &
                nodes  = opt_sdc % nodes,                       &
                po_f   = opt % p_time(l),                       &
                po_c   = opt % p_time(l-1),                     &
                mode   = opt % n_time(l) / opt % n_time(l-1),   &
                method = opt % projection_method                )

      else

        ! reset projection options
        opt_level % pop_fc_x = HP_CoarseningOptions_1D()
        opt_level % pop_fc_t = HP_CoarseningOptions_1D()

      end if

      ! create level  . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

      this % level(l) = CL_MLSDC_Level_1D(opt_level, cl_problem)

    end do

  end subroutine Init_CL_MLSDC_1D

  !=============================================================================
  ! TBP for CL_MLSDC_Options_1D

  !-----------------------------------------------------------------------------
  !> New CL_MLSDC_1D options

  function New_CL_MLSDC_Options_1D(n_level) result(opt)
    integer, intent(in) :: n_level !< number of space-time levels
    type(CL_MLSDC_Options_1D) :: opt

    opt % n_level = n_level

    allocate(opt % p_space(n_level))
    allocate(opt % n_space(n_level))
    allocate(opt % p_time (n_level))
    allocate(opt % n_time (n_level))
    allocate(opt % q_conv (n_level))

  end function New_CL_MLSDC_Options_1D

  !=============================================================================

end module CL__MLSDC__1D
