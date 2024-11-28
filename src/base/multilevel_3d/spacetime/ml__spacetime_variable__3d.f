!> summary:  3D multilevel spacetime mesh variable
!> author:   Joerg Stiller
!> date:     2024/11/24
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ML__Spacetime_Variable__3D
  use Mesh_Variable__3D
  use Spacetime_Variable__3D
  use ML__Spacetime_Operators__3D
  implicit none
  private

  public :: ML_SpacetimeVariable_3D

  !-----------------------------------------------------------------------------
  !> 3D multilevel spacetime mesh variable
  !>
  !> The values of a multilevel spacetime variable `u` are accessed via
  !>
  !>       u % level(l) % var(m,n) % val(i,j,k,e,c)
  !>
  !>  where
  !>
  !>    - `1 ≤   l   ≤ l_top`:  level
  !>    - `0 ≤   m   ≤ pt(l)`:  temporal collocation point ID
  !>    - `1 ≤   n   ≤ nt(l)`:  temporal element ID in given time slice
  !>    - `0 ≤ i,j,k ≤ po(l)`:  spatial collocation point triple index
  !>    - `1 ≤   e   ≤ ne(l)`:  spatial element ID
  !>    - `1 ≤   c   ≤ nc   `:  component ID

  type ML_SpacetimeVariable_3D
    type(SpacetimeVariable_3D), allocatable :: level(:) !< variable per level
    character(len=:),           allocatable :: name(:)  !< component names
  contains
    procedure :: Init_SpacetimeVariable_3D
  end type ML_SpacetimeVariable_3D

  ! constructor
  interface ML_SpacetimeVariable_3D
    procedure New_ML_SpacetimeVariable_3D
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Constructor of 3D multilevel spacetime variable

  function New_ML_SpacetimeVariable_3D(ml_op, nc, name) result(this)
    class(ML_SpacetimeOperators_3D), intent(in) :: ml_op
    integer,                         intent(in) :: nc
    character(len=*),      optional, intent(in) :: name(nc)
    type(ML_SpacetimeVariable_3D) :: this

    call Init_SpacetimeVariable_3D(this, ml_op, nc, name)

  end function New_ML_SpacetimeVariable_3D

  !-----------------------------------------------------------------------------
  !> Initialization of 3D multilevel spacetime variable

  subroutine Init_SpacetimeVariable_3D(this, ml_op, nc, name)
    class(ML_SpacetimeVariable_3D),  intent(inout) :: this
    class(ML_SpacetimeOperators_3D), intent(in)    :: ml_op
    integer,                         intent(in)    :: nc
    character(len=*),      optional, intent(in)    :: name(nc)

    integer :: c, l, m, n

    allocate(this % level( size(ml_op%sem) ))

    do l = 1, size(this%level)
      associate( po_t => this % level(l) % po_t &
               , ne_t => this % level(l) % ne_t )

        po_t = ml_op % sdc(l) % p_col
        ne_t = ml_op % ne_t(l)

        allocate(this % level(l) % var(0:po_t,1:ne_t))

        do n = 1, ne_t
        do m = 0, po_t
          this % level(l) % var(m,n) =                             &
              MeshVariable_3D( mesh = ml_op % sem(l) % mesh        &
                             , po   = ml_op % sem(l) % std_op % po &
                             , nc   = nc                           )
        end do
        end do

      end associate
    end do

    if (present(name)) then
      this % name = name
    else
      l = 2 + int(log10(dble(nc)))
      allocate(character(len=l) :: this % name(nc))
      do c = 1, nc
        write(this%name(c), '(A,I0)') 'v', c
      end do
    end if

  end subroutine Init_SpacetimeVariable_3D

  !=============================================================================

end module ML__Spacetime_Variable__3D
