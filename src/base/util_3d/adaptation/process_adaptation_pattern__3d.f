module Process_Adaptation_Pattern__3D
  use Mesh__3D
  use Mesh_Element_Indexing__3D
  use Element_Transfer_Buffer__3D

  implicit none
  private

  public :: ProcessAdaptationPattern_3D

contains

  !-----------------------------------------------------------------------------
  !> Convert adaptation flags to refinement marks and set number of sublevels

  subroutine ProcessAdaptationPattern_3D(mesh)
    class(Mesh_3D), intent(inout) :: mesh

    integer, allocatable, target :: mark(:)
    integer, contiguous, pointer :: mark_val(:,:,:,:)
    type(ElementTransferBuffer_3D), allocatable, asynchronous :: mark_buf

    logical :: comp_refined(26)
    integer :: e, i, i1, i2, k, l, nn

    ! adaptation marks .........................................................

    allocate(mark(mesh%n_elem + mesh%n_ghost), source = -1)

    ! extract adaptation marks
    do e = 1, mesh % n_elem
      if (mesh % element(e) % frozen) cycle
      mark(e) = mesh % element(e) % adaptation % mark
    end do

    ! transfer marks to ghosts
    if (mesh % n_ghost > 0) then
      mark_val(1:1,1:1,1:1,1:mesh%n_elem+mesh%n_ghost) => mark
      mark_buf = ElementTransferBuffer_3D(mesh, mark_val)
      call mark_buf % Transfer(mesh, mark_val, tag = 1731)
      call mark_buf % Merge(mark_val)
    end if

    ! set refinement marks .....................................................

    select case(mesh % child_type)

    case('s') ! global of local refinement by subividing

      do e = 1, mesh % n_elem
        associate(element => mesh % element(e))

          if (element % frozen) then

            ! no refinement
            element % adaptation % mark = 0

          else if (mark(e) > 0) then

            ! regular refinement
            element % adaptation % mark = 8000

          else

            comp_refined = .false.

            ! vertices with neighbors marked for refinement
            do k = 1, 8
              nn = element % vertex(k) % n_neighbor
              if (nn > 0) then
                l  = k + 18
                i1 = element % vertex(k) % i_neighbor
                i2 = i1 + nn - 1
                do i = i1, i2
                  if (mark(element % neighbor(i) % id) > 0) then
                    comp_refined(l) = .true.
                    exit
                  end if
                end do
              end if
            end do

          ! edges ...
          do k = 1, 12
            nn = element % edge(k) % n_neighbor
            if (nn > 0) then
              l  = k + 6
              i1 = element % edge(k) % i_neighbor
              i2 = i1 + nn - 1
              do i = i1, i2
                if (mark(element % neighbor(i) % id) > 0) then
                  comp_refined(l) = .true.
                  exit
                end if
              end do
              ! unmark edge vertices
              if (comp_refined(l)) then
                comp_refined(V_EDGE(1,k) + 18) = .false.
                comp_refined(V_EDGE(2,k) + 18) = .false.
              end if
            end if
          end do

          ! faces ...
          do k = 1, 6
            i = element % face(k) % i_neighbor
            if (i > 0) then
              comp_refined(k) = mark(element % neighbor(i) % id) > 0
            end if
            if (comp_refined(k)) then
              ! unmark face edges and vertices
              do i = 1, 4
                comp_refined(E_FACE(i,k) + 6 ) = .false.
                comp_refined(V_FACE(i,k) + 18) = .false.
              end do
            end if
          end do

          ! refinement mark
          select case(count(comp_refined))
          case(0)
            element % adaptation % mark = 0
          case(1)
            do l = 1, 26
              if (comp_refined(l)) exit
            end do
            select case(ElementComponentType(l))
            case(IS_VERTEX)
              element % adaptation % mark = 100 + ElementVertexID(l)
            case(IS_EDGE)
              element % adaptation % mark = 200 + ElementEdgeID(l)
            case(IS_FACE)
              element % adaptation % mark = 400 + ElementFaceID(l)
            end select
          case default
            ! selection of several component results in complete refinement
            element % adaptation % mark = 800
          end select

        end if

      end associate
    end do

    case('c') ! global or zonal cloning

      do e = 1, mesh % n_elem
        associate(element => mesh % element(e))
          if (element % frozen) then
            element % adaptation % mark = 0
          else if (mark(e) > 0) then
            element % adaptation % mark = 1000
          else if (any(mark(element % neighbor % id) > 0)) then
            element % adaptation % mark = 100
          else
            element % adaptation % mark = 0
          end if
        end associate
      end do

    end select

    ! set number of sublevels ..................................................

    do e = 1, mesh % n_elem
       mesh % element(e) % adaptation % sublevels = max(mark(e), 0) / 1000
    end do

  end subroutine ProcessAdaptationPattern_3D

  !=============================================================================

end module Process_Adaptation_Pattern__3D
