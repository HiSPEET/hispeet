!> summary:  Polynomial multigrid initialization routines
!> author:   Joerg Stiller
!> date:     2019/02/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(CART__Elliptic_PMG) MP_Init
  implicit none

contains

!-------------------------------------------------------------------------------
!> PMG_Method3D initialization: IP, constant isotropic diffusivity

module subroutine Init_IP( this, mesh, ip_opt, pmg_opt )
  class(PMG_Method3D),          intent(inout) :: this
  class(MeshPartition), target, intent(in)    :: mesh    !< mesh partition
  class(IP_ElementOptions1D),   intent(in)    :: ip_opt  !< IP/DG options
  class(PMG_Options3D),         intent(in)    :: pmg_opt !< PMG options

  integer, allocatable :: po(:)
  integer :: l, l_top, ns1, ns2

  if (ip_opt % po /= pmg_opt % po_top) then
    call Error('CART__Elliptic_PMG: Init_IP', 'ip_opt % po /= pmg_opt % po_top')
  end if

  this % mesh => mesh
  call Assign_PMG_Options(pmg_opt, this)
  call CreatePolynomialLevels(pmg_opt, po)

  l_top = ubound(po,1)
  allocate(this % level(0:l_top))

  ns1 = pmg_opt % ns1
  ns2 = pmg_opt % ns2

  associate(level => this % level)

    if (l_top > 0) then
      call level(l_top) % Init_TopLevel_IP( ns1, ns2            &
                                          , pmg_opt % smoother  &
                                          , mesh, ip_opt        &
                                          , pmg_opt % schwarz   &
                                          , po(l_top-1)         )
    else
      call level(l_top) % Init_TopLevel_IP( ns1, ns2            &
                                          , pmg_opt % smoother  &
                                          , mesh, ip_opt        &
                                          , pmg_opt % schwarz   )
    end if

    do l = l_top-1, 1, -1
      ns1 = ns1 * pmg_opt % mvs
      ns2 = ns2 * pmg_opt % mvs
      call level(l) % Init_CoarseLevel( po(l), ns1, ns2, this%level(l+1)  &
                                      , pmg_opt % schwarz, po(l-1)        )
    end do

    if (l_top > 0) then
      call level(0) % Init_CoarseLevel( po(0), ns1, ns2, this%level(1)  &
                                      , pmg_opt % schwarz               )
    end if

  end associate

end subroutine Init_IP

!-------------------------------------------------------------------------------
!> Assignment of options to PMG_Method3D object

subroutine Assign_PMG_Options(opt, pmg)
  class(PMG_Options3D), intent(in)    :: opt
  class(PMG_Method3D),  intent(inout) :: pmg

  pmg % i_max   = opt % i_max
  pmg % r_red   = opt % r_red
  pmg % r_max   = opt % r_max
  pmg % dr_min  = opt % dr_min

  pmg % solver  = opt % solver
  pmg % i0_max  = opt % i0_max
  pmg % r0_red  = opt % r0_red
  pmg % r0_max  = opt % r0_max
  pmg % monitor = opt % monitor

end subroutine Assign_PMG_Options

!-------------------------------------------------------------------------------
!> Creates a set of ascending polynomial orders

subroutine CreatePolynomialLevels(opt, po)
  class(PMG_Options3D), intent(in)  :: opt
  integer, allocatable, intent(out) :: po(:) !< polynomial levels

  integer, allocatable :: q(:)
  integer :: lb, lt

  associate( po_top => opt % po_top, cr     => opt % cr,    &
             po_bot => opt % po_bot, cr_max => opt % cr_max )

    if (po_top > po_bot) then

      allocate(q(po_bot:po_top), source = po_top)
      lt = po_top
      lb = lt
      do
        lb = lb -1
        q(lb) = max( nint(q(lb+1)/cr), ceiling(q(lb+1)/cr_max), 1 )
        q(lb) = min( q(lb), q(lb+1) - 1 )
        if (q(lb) <= max(po_bot, 1)) then
          q(lb) = max(q(lb), po_bot)
          exit
        end if
      end do
      allocate(po(0:lt-lb), source=q(lb:lt))

    else if (po_top == po_bot) then

      allocate(po(0:0), source=po_top)

    else

      allocate(po(0:-1))

    end if

  end associate

end subroutine CreatePolynomialLevels

!===============================================================================

end submodule MP_Init
