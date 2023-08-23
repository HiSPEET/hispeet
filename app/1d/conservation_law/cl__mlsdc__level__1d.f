!> summary:  MultiLevel Spectral Deferred Correction of 1D conservation laws
!> author:   Joerg Stiller
!> date:     2023/08/22
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CL__MLSDC__Level__1D
  use Kind_Parameters
! use Constants
! use Execution_Control
  use Coarse_To_Fine_Interpolation__1D
  use Fine_To_Coarse_Interpolation__1D
  use CL__Problem__1D
  use CL__Operator__1D
  use CL__SDC__Method__1D
  use CL__SDC__Method__Euler__1D
  use CL__SDC__Method__ISD__1D

  implicit none
  private

  !-----------------------------------------------------------------------------
  !> Type for keeping one level of the MLSDC data structure

  type, public :: CL_MLSDC_Level_1D

    logical :: is_root = .true.   !< T if root (bottom) level mesh
    logical :: is_top  = .true.   !< T if top level mesh

    ! discretization and solvers ...............................................

    class(CL_Problem_1D),    pointer     :: cl_problem   !< same for all levels
    class(CL_Operator_1D),   allocatable :: cl_operator
    class(CL_SDC_Method_1D), allocatable :: cl_sdc

    integer :: nt = 1  !< number of time steps in one space-time slice

    ! interpolation operators ..................................................

    type(CoarseToFineInterpolation_1D) :: iop_cf_x !< C-F interpolation in x
    type(CoarseToFineInterpolation_1D) :: iop_cf_t !< C-F interpolation in t

    type(FineToCoarseInterpolation_1D) :: iop_fc_x !< F-C interpolation in x
    type(FineToCoarseInterpolation_1D) :: iop_fc_t !< F-C interpolation in t

    ! data .....................................................................

    !! real(RNP), allocatable :: u(:,:,:,:,:)
    !! !< approximate solution dimensioned as `u(0:po,1:ne,1:nc,0:ns,1:nt)`, where
    !!  - 'po'  is the polynomial order (degree) in space
    !!  - 'ne'  is the number of elements in space
    !!  - 'nc'  is the number of solution components
    !!  - 'ns'  is the number of subintervals in each time step (degree in time)
    !!  - 'nt'  is the number of time steps

    ! possibly some more/auxiliary data

  contains

    procedure :: Interpolate_CF !< solution interpolation to next finer   level
    procedure :: Interpolate_FC !< solution interpolation to next coarser level
    procedure :: Restrict_FC    !< residual restriction   to next coarser level

  end type CL_MLSDC_Level_1D

contains

  !-----------------------------------------------------------------------------
  !> Coarse-to-fine space-time interpolation of solution-like variables

  subroutine Interpolate_CF(this, u_c, u_f, complete)
    class(CL_MLSDC_Level_1D), intent(in) :: this !< coarse level
    real(RNP), intent(in)    :: u_c(0:,:,:,0:,:) !< coarse solution variable
    real(RNP), intent(inout) :: u_f(0:,:,:,0:,:) !< fine solution variable
    logical,   intent(in)    :: complete !< T/F for all/active coarse elements

    real(RNP), allocatable :: u_i(:,:,:,:,:) ! intermediate interpolant

    integer :: po_c, po_f ! polynomial order in space
    integer :: ne_c, ne_f ! number of elements in space
    integer :: ns_c, ns_f ! number of subintervals in each time step
    integer :: nt_c, nt_f ! number of time steps
    integer :: nc         ! number of components, must not change

    logical :: is_consistent
    integer :: c, e, l, m, n, n1, n2

    associate( iop_x => this % iop_fc_x &
             , iop_t => this % iop_fc_t &
             , mask  => this % cl_operator % mask )

      ! initialization .........................................................

      ! coarse dimensions
      po_c = ubound(u_c, 1)
      ne_c = ubound(u_c, 2)
      ns_c = ubound(u_c, 4)
      nt_c = ubound(u_c, 5)

      ! fine dimensions
      po_f = ubound(u_f, 1)
      ne_f = ubound(u_f, 2)
      ns_f = ubound(u_f, 4)
      nt_f = ubound(u_f, 5)

      ! number of components
      nc = this % cl_problem % nc

      ! consistency of components
      is_consistent = ubound(u_c, 3) == nc .and. &
                      ubound(u_f, 3) == nc

      ! consistency of spatial dimensions
      is_consistent = is_consistent        .and. &
                      po_c == iop_x % po_c .and. &
                      po_f == iop_x % po_f

      ! consistency with spatial interpolation mode
      select case(iop_x % mode)
      case(1)
        is_consistent = is_consistent .and. ne_f == ne_c
      case(2)
        is_consistent = is_consistent .and. ne_f == ne_c * 2
      end select

      ! consistency of temporal dimensions
      is_consistent = is_consistent        .and. &
                      ns_c == iop_t % po_c .and. &
                      ns_f == iop_t % po_f

      ! consistency with spatial interpolation mode
      select case(iop_t % mode)
      case(1)
        is_consistent = is_consistent .and. nt_f == nt_c
      case(2)
        is_consistent = is_consistent .and. nt_f == nt_c * 2
      end select

      ! evaluate result of checks
      if (.not. is_consistent) then
        call Error( 'Interpolate_CF'           &
                  , 'failed consistency check' &
                  , 'CL__MLSDC__Level__1D'     )
      end if

      ! workspace
      allocate(u_i(0:po_f,1:ne_f,nc,0:ns_c,1:nt_c))

      ! spatial interpolation ..................................................

      select case(iop_x % mode)
      case(0)
        u_i = u_c
      case(1)
        ! ne_f = ne_c
        do n = 1, nt_c
        do m = 0, ns_c
        do c = 1, nc
        do e = 1, ne_c
          if (complete .or. mask(e)) then
            u_i(:,e,c,m,n) = matmul(iop_x % A(:,:,1), u_c(:,e,c,m,n))
          end if
        end do
        end do
        end do
        end do
      case(2)
        ! ne_f = 2 * ne_c
        do n = 1, nt_c
        do m = 0, ns_c
        do c = 1, nc
        do e = 1, ne_c
          if (complete .or. mask(e)) then
            u_i(:,2*e-1,c,m,n) = matmul(iop_x % A(:,:,1), u_c(:,e,c,m,n))
            u_i(:,2*e  ,c,m,n) = matmul(iop_x % A(:,:,2), u_c(:,e,c,m,n))
          end if
        end do
        end do
        end do
        end do
      end select

      ! temporal interpolation .................................................

      select case(iop_t % mode)
      case(0)
        u_f = u_i
      case(1)
        ! nt_f = nt_c
        do n = 1, nt_c
        do m = 0, ns_f
        do c = 1, nc
        do e = 1, ne_f
          if (complete .or. mask(e)) then
            u_f(:,e,c,m,n) = 0
            do l = 0, ns_c
              u_f(:,e,c,m,n) = u_f(:,e,c,m,n) + iop_t % A(m,l,1) * u_i(:,e,c,l,n)
            end do
          end if
        end do
        end do
        end do
        end do
      case(2)
        ! nt_f = 2 * nt_c
        do n = 1, nt_c
        do m = 0, ns_f
        do c = 1, nc
        do e = 1, ne_f
          if (complete .or. mask(e)) then
            n1 = 2 * n - 1
            n2 = 2 * n
            u_f(:,e,c,m,n1) = 0
            u_f(:,e,c,m,n2) = 0
            do l = 0, ns_c
              u_f(:,e,c,m,n1) = u_f(:,e,c,m,n1) + iop_t % A(m,l,1) * u_i(:,e,c,l,n)
              u_f(:,e,c,m,n2) = u_f(:,e,c,m,n2) + iop_t % A(m,l,2) * u_i(:,e,c,l,n)
            end do
          end if
        end do
        end do
        end do
        end do
      end select

    end associate

  end subroutine Interpolate_CF

  !=============================================================================

end module CL__MLSDC__Level__1D
