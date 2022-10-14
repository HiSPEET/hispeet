!> summary:  Export Schwarz domain data into VTK XML file
!> author:   Joerg Stiller
!> date:     2022/03/10
!> license:  Institute of Fluid Mechanics TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Export_VTK_Schwarz_Domain__3D
  use Kind_Parameters
  use Export_VTK_Volume_Data__3D
  implicit none
  private

  public :: ExportVTK_SchwarzDomain

contains

  !-----------------------------------------------------------------------------
  !> Export cuboidal Schwarz domain with scalar data to VTK

  subroutine ExportVTK_SchwarzDomain(x_cube, xi, v_sd, file)
    real(RNP),        intent(in) :: x_cube(0:,:) !< cuboidal core domain
    real(RNP),        intent(in) :: xi(0:)       !< standard core points
    class(*),         intent(in) :: v_sd(:,:,:)  !< subdomain variable (RSP,RDP)
    character(len=*), intent(in) :: file         !< VTK output file

    real(RNP), allocatable :: xi_sd(:), x(:,:,:,:,:), v(:,:,:,:,:)
    integer :: ne, no, np, ns, po
    integer :: d, e, i, j, k, l, m, n

    !$omp master

    ! dimensions ...............................................................

    po = ubound(xi,1)    ! polynomial degree
    np = po + 1          ! number of core points per direction
    ns = size(v_sd,1)    ! number of subdomain points per direction
    no = (ns - np) / 2   ! overlap
    ne = (ns - 1)**3     ! number of exported linear elements

    ! workspace ................................................................

    allocate(xi_sd(ns), x(0:1,0:1,0:1,ne,3), v(0:1,0:1,0:1,ne,1))

    ! subdomain points .........................................................

    ! standard coordinates of subdomain points
    xi_sd = [ xi(np-no:po) - 2, xi(0:po), xi(0:no-1) + 2 ]

    ! subdomain element points
    e = 0
    do k = 1, ns - 1
    do j = 1, ns - 1
    do i = 1, ns - 1
      e = e + 1
      do d = 1, 3
        do n = 0, 1
        do m = 0, 1
        do l = 0, 1
          x(l,m,n,e,d) =  x_cube(0,d)                &
                       +  x_cube(1,d) * xi_sd(i + l) &
                       +  x_cube(2,d) * xi_sd(j + m) &
                       +  x_cube(3,d) * xi_sd(k + n)
        end do
        end do
        end do
      end do
    end do
    end do
    end do

    ! subdomain variable .......................................................

    e = 0

    select type(v_sd)

    type is (real(RSP))
      do k = 1, ns - 1
      do j = 1, ns - 1
      do i = 1, ns - 1
        e = e + 1
        do n = 0, 1
        do m = 0, 1
        do l = 0, 1
          v(l,m,n,e,1) = v_sd(i+l,j+m,k+n)
        end do
        end do
        end do
      end do
      end do
      end do

    type is (real(RDP))
      do k = 1, ns - 1
      do j = 1, ns - 1
      do i = 1, ns - 1
        e = e + 1
        do n = 0, 1
        do m = 0, 1
        do l = 0, 1
          v(l,m,n,e,1) = v_sd(i+l,j+m,k+n)
        end do
        end do
        end do
      end do
      end do
      end do

    end select

    call ExportVTK_VolumeData(x, s=v, sname=['v'], file=file, subdiv = .false.)

    !$omp end master

  end subroutine ExportVTK_SchwarzDomain

  !=============================================================================

end module Export_VTK_Schwarz_Domain__3D
