!> summary:  Creation of a cuboid domain with regular Cartesian mesh
!> author:   Joerg Stiller
!> date:     2025/09/24
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================
!> summary:  Creation of a spherical shell segment with structured mesh
!> author:   Joerg Stiller
!> date:     2025/09/23
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> The shell segment is specified via mean spherical coordinates r₀, φ₀, θ₀
!> and extensions Δr, Δφ, Δθ such that
!>
!>       r₀ - Δr/2  ≤  r  ≤  r₀ + δr/2
!>       φ₀ - Δφ/2  ≤  φ  ≤  φ₀ + δφ/2
!>       θ₀ - Δθ/2  ≤  θ  ≤  θ₀ + δθ/2
!>
!> The following conditions must be satisfied:
!>
!>       0  <  δr  <  2 r₀
!>       0  <  δφ  <  2 π
!>       0  <  δθ  <    π
!===============================================================================

module Create_Spherical_Shell_Segment
  use Kind_Parameters
  use Constants
  use XMPI
  use Mesh__3D
  use Generate_Regular_Mesh__3D
  implicit none
  private

  public :: CreateSphericalShellSegment

contains

  !-----------------------------------------------------------------------------
  !> Creation of a cuboid domain with regular Cartesian mesh

  subroutine CreateSphericalShellSegment(comm, input, mesh)

    ! arguments ................................................................

    type(MPI_Comm),   intent(in)  :: comm  !< MPI communicator
    character(len=*), intent(in)  :: input !< input file
    type(Mesh_3D),    intent(out) :: mesh  !< spectral-element mesh

    ! input parameters .........................................................

    real(RNP) :: r0 = ONE  ! r₀
    real(RNP) :: p0 = ZERO ! φ₀
    real(RNP) :: t0 = PI/2 ! θ₀
    real(RNP) :: dr = ONE  ! δr
    real(RNP) :: dp = ONE  ! δφ
    real(RNP) :: dt = ONE  ! δθ

    namelist/spherical_shell_segment_prm/ r0, p0, t0, dr, dp, dt

    integer :: np(3) = 1           ! number of partitions per direction
    integer :: ep(3) = 2           ! elements per partition and direction
    integer :: pg    = 8           ! polynomial order of geometry
    logical :: structured = .true. ! define mesh to be (un)structured

    namelist/spherical_shell_segment_prm/ np, ep, pg, structured

    ! auxiliary variables ......................................................

    real(RNP) :: xo(3), dx(3)
    real(RNP) :: r, p, t
    logical   :: periodic(3) = .false.
    integer   :: rank, io
    integer   :: e, i, j, k

    ! initialization ...........................................................

    call MPI_Comm_rank(comm, rank)

    if (rank == 0) then
      open(newunit = io, file = input)
      read(io, nml = spherical_shell_segment_prm)
      close(io)
    end if

    call XMPI_Bcast(r0        , 0, comm)
    call XMPI_Bcast(p0        , 0, comm)
    call XMPI_Bcast(t0        , 0, comm)
    call XMPI_Bcast(dr        , 0, comm)
    call XMPI_Bcast(dp        , 0, comm)
    call XMPI_Bcast(dt        , 0, comm)
    call XMPI_Bcast(np        , 0, comm)
    call XMPI_Bcast(ep        , 0, comm)
    call XMPI_Bcast(pg        , 0, comm)
    call XMPI_Bcast(structured, 0, comm)

    ! Cartesian mesh in [-1,1]³ ................................................

    xo = -ONE
    dx =  TWO / (np * ep)

    call GenerateRegularMesh(mesh, np, ep, xo, dx, periodic, comm, pg)

    mesh % regular    = .false.
    mesh % structured = structured
    mesh % dx         = 0

    ! transform to spherical shell .............................................

    do e = 1, mesh % n_elem
      associate(x_e => mesh % element(e) % geometry % x_e)
        do k = 0, pg
        do j = 0, pg
        do i = 0, pg
          r = r0 + HALF * dr * x_e(i,j,k,1)
          p = p0 + HALF * dp * x_e(i,j,k,2)
          t = t0 - HALF * dt * x_e(i,j,k,3)
          x_e(i,j,k,1) = r * sin(t) * cos(p)
          x_e(i,j,k,2) = r * sin(t) * sin(p)
          x_e(i,j,k,3) = r * cos(t)
        end do
        end do
        end do
      end associate
    end do

    call mesh % BuildCuboids()

  end subroutine CreateSphericalShellSegment

  !=============================================================================

end module Create_Spherical_Shell_Segment
