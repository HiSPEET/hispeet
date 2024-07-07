!> summary:  Exporting a multilevel mesh variable to VTK
!> author:   Joerg Stiller
!> date:     2024/07/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(ML__Mesh_Variable__3D) MP_EportVTK
  use Export_VTK_Volume_Data__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Export multilevel mesh variable to VTK
  !>
  !> Supported modes
  !>   - `1`  all elements
  !>   - `2`  only active elements
  !>   - `3`  only leaf elements

  module subroutine ExportVTK(this, ml_op, file, mode)
    class(ML_MeshVariable_3D),  intent(in) :: this  !< multilevel variable
    class(ML_MeshOperators_3D), intent(in) :: ml_op !< spectral element OPs
    character(len=*),           intent(in) :: file  !< export file base name
    integer,                    intent(in) :: mode  !< export mode {1,2,3}

    character(len=:), allocatable, save :: tag
    logical, allocatable, save :: mask(:)
    integer :: l, m

    do l = 1, size(this%level)
      associate(sem => ml_op%sem(l))

        if (sem%mesh%part < 0) cycle

        ! element export mask
        allocate(mask(sem%mesh%n_elem))
        select case(mode)
        case(1)
          mask = .true.
        case(2)
          mask = .not. sem % mesh % element % frozen
        case(3)
          mask = (.not. sem % mesh % element % frozen) .and. &
                 (sem % mesh % element % adaptation % refinement < 1000)
        end select

        ! level tag added to file name
        m = 3 + int(log10(dble(l)))
        allocate(character(len=m) :: tag)
        write(tag,'(A,I0)') '_l', l

        call ExportVTK_VolumeData( x       = sem%metrics%x     &
                                 , s       = this%level(l)%val &
                                 , sname   = this%name         &
                                 , file    = trim(file)//tag   &
                                 , part    = sem%mesh%part     &
                                 , n_parts = sem%mesh%n_parts  &
                                 , mask    = mask              )

        deallocate(tag, mask)

      end associate
    end do

  end subroutine ExportVTK

  !=============================================================================

end submodule MP_EportVTK
