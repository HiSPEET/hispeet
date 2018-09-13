!> summary:  Cartesian diffusion operator, constant isotropic, structured
!> author:   Joerg Stiller
!> date:     2016/12/08, revised 2017/05/04-
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Cartesian diffusion operator, constant isotropic, structured
!===============================================================================

module CART__DG_Diffusion_CIS_Operator

!### CHECK
!  use OpenMP_Binding
!  use XMPI
!### END CHECK

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE, HALF
  use Array_Assignments

  use CART__TPO_Diffusion
  use CART__DG_Element_Operators
  use CART__Mesh_Partition
  use CART__Face_Transfer_Buffer

  implicit none
  private

  public :: DiffusionOperator

contains

!-------------------------------------------------------------------------------
!> Diffusion operator, v = lambda M u + nu L u

subroutine DiffusionOperator(mesh, eop, lambda, nu, bc, u, v)

  ! arguments ..................................................................

  ! *** ACC: Assuming present(eop,bc,u,v)
  ! *** ACC: Check that all data required by GPU is on device (derived types?)

  class(MeshPartition),       intent(in)  :: mesh !< mesh partition
  class(DG_ElementOperators), intent(in)  :: eop  !< DG element operators

  real(RNP), intent(in)  :: lambda  !< Helmholtz parameter
  real(RNP), intent(in)  :: nu      !< diffusivity
  character, intent(in)  :: bc(6)   !< boundary conditions {P,D,N}
  real(RNP), intent(in)  :: u       !< approximate solution
  real(RNP), intent(out) :: v       !< result

  dimension :: u(0:eop%po, 0:eop%po, 0:eop%po, mesh%ne1, mesh%ne2, mesh%ne3)
  dimension :: v(0:eop%po, 0:eop%po, 0:eop%po, mesh%ne1, mesh%ne2, mesh%ne3)

  ! local variables ............................................................

  procedure(TPO_Diffusion_Proc), pointer, save :: StiffnessOperator

  ! subdomain boundary conditions
  character :: sdbc(6)

  ! work arrays for [u] and {du/dn} on mesh faces
  real(RNP), allocatable, save :: J_u(:,:,:), D_u(:,:,:)

  ! MPI transfer buffers -- declared save to become shared with OpenMP
  type(FaceTransferBuffer), allocatable, save :: buf__D_u  ! D_u transfer buffer
  type(FaceTransferBuffer), allocatable, save :: buf__J_u  ! J_u transfer buffer

  integer :: po, ne(3), np = -1
  integer :: i1, i2, i3
  integer :: j1, j2, j3
  integer :: i

!### CHECK
!  integer :: mpi_id, omp_id
!  mpi_id = mesh%part
!  omp_id = OMP_Get_Thread_Num()
!### END CHECK

  ! initialization .............................................................
!### CHECK
! block
! if (any(isNaN(u))) then
!   print '(9G0)', '@DO1[',mesh%comm,']: u has NaN'
! end if
! call MPI_Barrier(mesh%comm)
! end block
!### END CHECK

  ! mesh dimensions
  po    = ubound(u,3)  ! polynomial order
  ne(1) = ubound(u,4)  ! number of elements in direction 1
  ne(2) = ubound(u,5)  ! number of elements in direction 2
  ne(3) = ubound(u,6)  ! number of elements in direction 3

  ! procedure for evaluating the element operators
  if (np /= po + 1) then
    np  = po + 1
    call TPO_Diffusion_Assign(np, StiffnessOperator)
  end if

  ! subdomain BC
  sdbc = 'I'
  do i = 1, 6
    if (size(mesh%boundary(i)%face) > 0) sdbc(i) = bc(i)
  end do

  ! arrays and buffers for [u] and {du/dn}
  !$omp single
  allocate(J_u(0:po, 0:po, mesh%nf))
  allocate(D_u(0:po, 0:po, mesh%nf))
  allocate(buf__D_u)
  allocate(buf__J_u)
  !$omp end single

  call AssignScalar(J_u, ZERO)
  call AssignScalar(D_u, ZERO)

  ! index range of x1-faces
  i1 = 1
  j1 = i1 - 1 + mesh%nf1

  ! index range of x2-faces
  i2 = j1 + 1
  j2 = j1 + mesh%nf2

  ! index range of x3-faces
  i3 = j2 + 1
  j3 = j2 + mesh%nf3

  ! *** ACC: possibly include required eop components into following region

  !$acc data copyin(ne, sdbc) create(J_u, D_u)

  ! local traces of u and du/dx ................................................

  call TraceOperators( eop, po, ne, u,                                &
                       J_u(:,:,i1:j1), J_u(:,:,i2:j2), J_u(:,:,i3:j3), &
                       D_u(:,:,i1:j1), D_u(:,:,i2:j2), D_u(:,:,i3:j3)  )

  ! initialize exchange with neighbors .........................................

  ! start transfer operations
  ! *** ACC: Beware that all MPI operations will be executed only on host (CPU)!
  ! *** ACC: Take care that host is updated
  call buf__D_u % Transfer(mesh, D_u, tag=1000)
  call buf__J_u % Transfer(mesh, J_u, tag=2000)

  ! apply element operators ....................................................

  associate( Ms => eop % w, &
             Ls => eop % L, &
             dx => eop % dx )

    call StiffnessOperator(np, product(ne), Ms, Ls, lambda, nu, dx, u, v)

  end associate
!### CHECK
! block
! if (any(isNaN(v))) then
!   print '(9G0)', '@DO1[',mesh%comm,']: v has NaN'
! end if
! call MPI_Barrier(mesh%comm)
! end block
!### END CHECK

  ! apply boundary conditions ..................................................

  call ApplyBC1(po, ne, sdbc(1:2), J_u(:,:,i1:j1), D_u(:,:,i1:j1))
  call ApplyBC2(po, ne, sdbc(3:4), J_u(:,:,i2:j2), D_u(:,:,i2:j2))
  call ApplyBC3(po, ne, sdbc(5:6), J_u(:,:,i3:j3), D_u(:,:,i3:j3))

  ! finalize recv and add remote values ........................................

  call buf__D_u % Merge(mesh, D_u)
  call buf__J_u % Merge(mesh, J_u)

  ! add fluxes .................................................................

  ! wait for the stiffness operator to complete
  !$acc wait

  call AddFluxes( eop, po, ne, nu,                               &
                  J_u(:,:,i1:j1), J_u(:,:,i2:j2), J_u(:,:,i3:j3), &
                  D_u(:,:,i1:j1), D_u(:,:,i2:j2), D_u(:,:,i3:j3), &
                  v                                               )

  !$acc end data
!### CHECK
! block
! if (any(isNaN(v))) then
!   print '(9G0)', '@DO2[',mesh%comm,']: v has NaN'
! end if
! call MPI_Barrier(mesh%comm)
! end block
!### END CHECK

  ! wait for send operations to complete .......................................

  call buf__D_u % Finish()
  call buf__J_u % Finish()

  ! clean-up ...................................................................

  !$omp single
  deallocate(J_u, D_u)
  deallocate(buf__J_u)
  deallocate(buf__D_u)
  !$omp end single

end subroutine DiffusionOperator

!-------------------------------------------------------------------------------
!> Traces of u and du/dn on local element faces

subroutine TraceOperators(eop, po, ne, u, J1_u, J2_u, J3_u, D1_u, D2_u, D3_u)

  ! arguments ..................................................................

  class(DG_ElementOperators), intent(in) :: eop !< ! element operators

  integer, intent(in) :: po     !< polynomial order (for convenience)
  integer, intent(in) :: ne(3)  !< number of elements per direction

  real(RNP), intent(in)  :: u(0:po, 0:po, 0:po, ne(1), ne(2), ne(3))

  real(RNP), intent(out) :: J1_u(0:po, 0:po, 0:ne(1), 1:ne(2), 1:ne(3))
  real(RNP), intent(out) :: J2_u(0:po, 0:po, 1:ne(1), 0:ne(2), 1:ne(3))
  real(RNP), intent(out) :: J3_u(0:po, 0:po, 1:ne(1), 1:ne(2), 0:ne(3))

  real(RNP), intent(out) :: D1_u(0:po, 0:po, 0:ne(1), 1:ne(2), 1:ne(3))
  real(RNP), intent(out) :: D2_u(0:po, 0:po, 1:ne(1), 0:ne(2), 1:ne(3))
  real(RNP), intent(out) :: D3_u(0:po, 0:po, 1:ne(1), 1:ne(2), 0:ne(3))

  ! local variables ............................................................

  real(RNP), allocatable, save :: T1_u (:,:,:,:,:,:) ! u      trace on x1-faces
  real(RNP), allocatable, save :: T1_du(:,:,:,:,:,:) ! du/dx1 trace on x1-faces
  real(RNP), allocatable, save :: T2_u (:,:,:,:,:,:) ! u      trace on x2-faces
  real(RNP), allocatable, save :: T2_du(:,:,:,:,:,:) ! du/dx2 trace on x2-faces
  real(RNP), allocatable, save :: T3_u (:,:,:,:,:,:) ! u      trace on x3-faces
  real(RNP), allocatable, save :: T3_du(:,:,:,:,:,:) ! du/dx3 trace on x3-faces

  real(RNP) :: g(3), us(0:po,0:po), tmp1, tmp2
  integer   :: i, j, k, l, m, n

  ! intialization ..............................................................

  !$omp single

  allocate( T1_u (2, 0:po, 0:po, 0:ne(1), 1:ne(2), 1:ne(3)))
  allocate( T2_u (2, 0:po, 0:po, 1:ne(1), 0:ne(2), 1:ne(3)))
  allocate( T3_u (2, 0:po, 0:po, 1:ne(1), 1:ne(2), 0:ne(3)))

  allocate( T1_du(2, 0:po, 0:po, 0:ne(1), 1:ne(2), 1:ne(3)))
  allocate( T2_du(2, 0:po, 0:po, 1:ne(1), 0:ne(2), 1:ne(3)))
  allocate( T3_du(2, 0:po, 0:po, 1:ne(1), 1:ne(2), 0:ne(3)))

  !$omp end single

  call AssignScalar(T1_u , ZERO)
  call AssignScalar(T2_u , ZERO)
  call AssignScalar(T3_u , ZERO)
  call AssignScalar(T1_du, ZERO)
  call AssignScalar(T2_du, ZERO)
  call AssignScalar(T3_du, ZERO)

  ! generate traces ............................................................

  ! *** ACC: Take care that T1_u, ... are created on device
  ! *** ACC/OMP: Set T1_u, ... to ZERO using AssignScalar from module
  ! ***          Array_Assignments; remove source = ZERO from allocation

  associate( Ds => eop%D, dx => eop%dx)

    g = 2 / dx

    ! *** OMP/ACC: parallel

    !$omp do collapse(3)
    do n = 1, ne(3)
    do m = 1, ne(2)
    do l = 1, ne(1)

      ! x3 face: u bottom

      do j = 0, po
      do i = 0, po
        T3_u(2, i,j, l,m,n-1) = u(i,j,0,l,m,n)
      end do
      end do

      k_planes: do k = 0, po

        ! extract k-slice
        do j = 0, po
        do i = 0, po
          us(i,j) = u(i,j,k,l,m,n)    ! reuse cached data when k = 0
        end do
        end do

        ! u west/east
        do j = 0, po
          T1_u(2, j,k, l-1,m,n) = us( 0,j)
          T1_u(1, j,k, l  ,m,n) = us(po,j)
        end do

        ! du/dx1 west/east
        do j = 0, po
          tmp1 = ZERO
          tmp2 = ZERO
          do i = 0, po
            tmp1 = tmp1 + Ds( 0,i) * us(i,j)
            tmp2 = tmp2 + Ds(po,i) * us(i,j)
          end do
          T1_du(2, j,k, l-1,m,n) = g(1) * tmp1
          T1_du(1, j,k, l  ,m,n) = g(1) * tmp2
        end do

        ! u south/north
        do i = 0, po
          T2_u(2, i,k, l,m-1,n) = us(i, 0)
          T2_u(1, i,k, l,m  ,n) = us(i,po)
        end do

        ! du/dx2 south/north
        do i = 0, po
          tmp1 = ZERO
          tmp2 = ZERO
          do j = 0, po
            tmp1 = tmp1 + Ds( 0,j) * us(i,j)
            tmp2 = tmp2 + Ds(po,j) * us(i,j)
          end do
          T2_du(2, i,k, l,m-1,n) = g(2) * tmp1
          T2_du(1, i,k, l,m  ,n) = g(2) * tmp2
        end do

      end do k_planes

      ! x3 face: u top
      do j = 0, po
      do i = 0, po
        T3_u(1,i,j,l,m,n) = us(i,j)    ! reuse cached data
      end do
      end do

      ! du/dx3 bottom/top
      do j = 0, po
      do i = 0, po
        tmp1 = ZERO
        tmp2 = ZERO
        do k = 0, po
          tmp1 = tmp1 + Ds( 0,k) * u(i,j,k,l,m,n)
          tmp2 = tmp2 + Ds(po,k) * u(i,j,k,l,m,n)
        end do
        T3_du(2, i,j, l,m,n-1) = g(3) * tmp1
        T3_du(1, i,j, l,m,n  ) = g(3) * tmp2
      end do
      end do

    end do
    end do
    end do
    !$omp end do

    ! apply trace operators ....................................................

    ! x1
    !$omp do collapse(2)
    do n = 1, ne(3)
    do m = 1, ne(2)
    do l = 0, ne(1)
      do k = 0, po
      do j = 0, po
        J1_u(j,k,l,m,n) =  T1_u (1,j,k,l,m,n) - T1_u (2,j,k,l,m,n)
        D1_u(j,k,l,m,n) = (T1_du(1,j,k,l,m,n) + T1_du(2,j,k,l,m,n)) * HALF
      end do
      end do
    end do
    end do
    end do
    !$omp end do nowait

    ! x2
    !$omp do collapse(2)
    do n = 1, ne(3)
    do m = 0, ne(2)
    do l = 1, ne(1)
      do k = 0, po
      do i = 0, po
        J2_u(i,k,l,m,n) =  T2_u (1,i,k,l,m,n) - T2_u (2,i,k,l,m,n)
        D2_u(i,k,l,m,n) = (T2_du(1,i,k,l,m,n) + T2_du(2,i,k,l,m,n)) * HALF
      end do
      end do
    end do
    end do
    end do
    !$omp end do nowait

    ! x3
    !$omp do collapse(2)
    do n = 0, ne(3)
    do m = 1, ne(2)
    do l = 1, ne(1)
      do j = 0, po
      do i = 0, po
        J3_u(i,j,l,m,n) =  T3_u (1,i,j,l,m,n) - T3_u (2,i,j,l,m,n)
        D3_u(i,j,l,m,n) = (T3_du(1,i,j,l,m,n) + T3_du(2,i,j,l,m,n)) * HALF
      end do
      end do
    end do
    end do
    end do
    !$omp end do

  end associate

  !$omp single
  deallocate(T1_u, T1_du, T2_u, T2_du, T3_u, T3_du)
  !$omp end single


end subroutine TraceOperators

!-------------------------------------------------------------------------------
!> Apply boundary conditions in x1 direction

subroutine ApplyBC1(po, ne, bc, J1_u, D1_u)

  ! arguments ..................................................................

  integer,   intent(in)    :: po    !< polynomial order
  integer,   intent(in)    :: ne(3) !< number of elements per direction
  character, intent(in)    :: bc(2) !< BC on west and east boundaries
  real(RNP), intent(inout) :: J1_u(0:po,0:po,0:ne(1),ne(2),ne(3)) !< [u] @ F1
  real(RNP), intent(inout) :: D1_u(0:po,0:po,0:ne(1),ne(2),ne(3)) !< {∇u} @ F1

  ! local variables ............................................................

  real(RNP) :: cb
  integer   :: j, k, m, n

  ! periodic BC ................................................................

  if (all(bc == 'P')) then

    ! *** OMP/ACC: parallel

    !$omp do collapse(2)
    do n = 1, ne(3)
    do m = 1, ne(2)
      do k = 0, po
      do j = 0, po

        J1_u(j,k,    0,m,n)  =  J1_u(j,k,0,m,n) + J1_u(j,k,ne(1),m,n)
        J1_u(j,k,ne(1),m,n)  =  J1_u(j,k,0,m,n)

        D1_u(j,k,    0,m,n)  =  D1_u(j,k,0,m,n) + D1_u(j,k,ne(1),m,n)
        D1_u(j,k,ne(1),m,n)  =  D1_u(j,k,0,m,n)

      end do
      end do
    end do
    end do
    !$omp end do

  end if

  ! western boundary ...........................................................

  if (bc(1) == 'D' .or. bc(1) == 'N') then

    if (bc(1) == 'D') then
      cb = 2
    else
      cb = 0
    end if

    !$omp do collapse(2)
    do n = 1, ne(3)
    do m = 1, ne(2)
      do k = 0, po
      do j = 0, po

        J1_u(j,k,0,m,n)  =  cb * J1_u(j,k,0,m,n)
        D1_u(j,k,0,m,n)  =  cb * D1_u(j,k,0,m,n)

      end do
      end do
    end do
    end do
    !$omp end do

  end if

  ! eastern boundary ...........................................................

  if (bc(2) == 'D' .or. bc(2) == 'N') then

    if (bc(2) == 'D') then
      cb = 2
    else
      cb = 0
    end if

    !$omp do collapse(2)
    do n = 1, ne(3)
    do m = 1, ne(2)
      do k = 0, po
      do j = 0, po

        J1_u(j,k,ne(1),m,n)  =  cb * J1_u(j,k,ne(1),m,n)
        D1_u(j,k,ne(1),m,n)  =  cb * D1_u(j,k,ne(1),m,n)

      end do
      end do
    end do
    end do
    !$omp end do

  end if

end subroutine ApplyBC1

!-------------------------------------------------------------------------------
!> Apply boundary conditions in x2 direction

subroutine ApplyBC2(po, ne, bc, J2_u, D2_u)

  ! arguments ..................................................................

  integer,   intent(in)    :: po    !< polynomial order
  integer,   intent(in)    :: ne(3) !< number of elements per direction
  character, intent(in)    :: bc(2) !< BC on west and east boundaries
  real(RNP), intent(inout) :: J2_u(0:po,0:po,ne(1),0:ne(2),ne(3)) !< [u]  @ F2
  real(RNP), intent(inout) :: D2_u(0:po,0:po,ne(1),0:ne(2),ne(3)) !< {∇u} @ F2

  ! local variables ............................................................

  real(RNP) :: cb
  integer   :: i, k, l, n

  ! periodic BC ................................................................

  if (all(bc == 'P')) then

    !$omp do collapse(2)
    do n = 1, ne(3)
    do l = 1, ne(1)
      do k = 0, po
      do i = 0, po

        J2_u(i,k,l,    0,n)  =  J2_u(i,k,l,0,n) + J2_u(i,k,l,ne(2),n)
        J2_u(i,k,l,ne(2),n)  =  J2_u(i,k,l,0,n)

        D2_u(i,k,l,    0,n)  =  D2_u(i,k,l,0,n) + D2_u(i,k,l,ne(2),n)
        D2_u(i,k,l,ne(2),n)  =  D2_u(i,k,l,0,n)

      end do
      end do
    end do
    end do
    !$omp end do

  end if

  ! southern boundary ..........................................................

  if (bc(1) == 'D' .or. bc(1) == 'N') then

    if (bc(1) == 'D') then
      cb = 2
    else
      cb = 0
    end if

    !$omp do collapse(2)
    do n = 1, ne(3)
    do l = 1, ne(1)
      do k = 0, po
      do i = 0, po

        J2_u(i,k,l,0,n)  =  cb * J2_u(i,k,l,0,n)
        D2_u(i,k,l,0,n)  =  cb * D2_u(i,k,l,0,n)

      end do
      end do
    end do
    end do
    !$omp end do

  end if

  ! northern boundary ..........................................................

  if (bc(2) == 'D' .or. bc(2) == 'N') then

    if (bc(2) == 'D') then
      cb = 2
    else
      cb = 0
    end if

    !$omp do collapse(2)
    do n = 1, ne(3)
    do l = 1, ne(1)
      do k = 0, po
      do i = 0, po

        J2_u(i,k,l,ne(2),n)  =  cb * J2_u(i,k,l,ne(2),n)
        D2_u(i,k,l,ne(2),n)  =  cb * D2_u(i,k,l,ne(2),n)

      end do
      end do
    end do
    end do
    !$omp end do

  end if

end subroutine ApplyBC2

!-------------------------------------------------------------------------------
!> Apply boundary conditions in x3 direction

subroutine ApplyBC3(po, ne, bc, J3_u, D3_u)

  ! arguments ..................................................................

  integer,   intent(in)    :: po    !< polynomial order
  integer,   intent(in)    :: ne(3) !< number of elements per direction
  character, intent(in)    :: bc(2) !< BC on west and east boundaries
  real(RNP), intent(inout) :: J3_u(0:po,0:po,ne(1),ne(2),0:ne(3)) !< [u]  @ F3
  real(RNP), intent(inout) :: D3_u(0:po,0:po,ne(1),ne(2),0:ne(3)) !< {∇u} @ F3

  ! local variables ............................................................

  real(RNP) :: cb
  integer   :: i, j, l, m

  ! periodic BC ................................................................

  if (all(bc == 'P')) then

    !$omp do collapse(2)
    do m = 1, ne(2)
    do l = 1, ne(1)
      do j = 0, po
      do i = 0, po

        J3_u(i,j,l,m,    0)  =  J3_u(i,j,l,m,0) + J3_u(i,j,l,m,ne(3))
        J3_u(i,j,l,m,ne(3))  =  J3_u(i,j,l,m,0)

        D3_u(i,j,l,m,    0)  =  D3_u(i,j,l,m,0) + D3_u(i,j,l,m,ne(3))
        D3_u(i,j,l,m,ne(3))  =  D3_u(i,j,l,m,0)

      end do
      end do
    end do
    end do
    !$omp end do

  end if

  ! bottom boundary ............................................................

  if (bc(1) == 'D' .or. bc(1) == 'N') then

    if (bc(1) == 'D') then
      cb = 2
    else
      cb = 0
    end if

    !$omp do collapse(2)
    do m = 1, ne(2)
    do l = 1, ne(1)
      do j = 0, po
      do i = 0, po

        J3_u(i,j,l,m,0)  =  cb * J3_u(i,j,l,m,0)
        D3_u(i,j,l,m,0)  =  cb * D3_u(i,j,l,m,0)

      end do
      end do
    end do
    end do
    !$omp end do

  end if

  ! top boundary ...............................................................

  if (bc(2) == 'D' .or. bc(2) == 'N') then

    if (bc(2) == 'D') then
      cb = 2
    else
      cb = 0
    end if

    !$omp do collapse(2)
    do m = 1, ne(2)
    do l = 1, ne(1)
      do j = 0, po
      do i = 0, po

        J3_u(i,j,l,m,ne(3))  =  cb * J3_u(i,j,l,m,ne(3))
        D3_u(i,j,l,m,ne(3))  =  cb * D3_u(i,j,l,m,ne(3))

      end do
      end do
    end do
    end do
    !$omp end do

  end if

end subroutine ApplyBC3

!-------------------------------------------------------------------------------
!> Compute and add fluxes through element boundaries

subroutine AddFluxes(eop, po, ne, nu, J1_u, J2_u, J3_u, D1_u, D2_u, D3_u, v)

  ! arguments ..................................................................

  class(DG_ElementOperators), intent(in) :: eop !< ! element operators

  integer,   intent(in) :: po     !< polynomial order
  integer,   intent(in) :: ne(3)  !< number of elements per direction
  real(RNP), intent(in) :: nu     !< diffusivity

  real(RNP), intent(in) :: J1_u(0:po, 0:po, 0:ne(1), 1:ne(2), 1:ne(3))
  real(RNP), intent(in) :: J2_u(0:po, 0:po, 1:ne(1), 0:ne(2), 1:ne(3))
  real(RNP), intent(in) :: J3_u(0:po, 0:po, 1:ne(1), 1:ne(2), 0:ne(3))

  real(RNP), intent(in) :: D1_u(0:po, 0:po, 0:ne(1), 1:ne(2), 1:ne(3))
  real(RNP), intent(in) :: D2_u(0:po, 0:po, 1:ne(1), 0:ne(2), 1:ne(3))
  real(RNP), intent(in) :: D3_u(0:po, 0:po, 1:ne(1), 1:ne(2), 0:ne(3))

  real(RNP), intent(inout) :: v(0:po, 0:po, 0:po, ne(1), ne(2), ne(3))

  ! local variables ............................................................

  real(RNP), dimension(:,:), allocatable :: M1, M2, M3
  real(RNP), dimension(:),   allocatable :: delta_0, delta_P
  real(RNP) :: c(3)
  integer :: i, j, k, l, m, n

  !.............................................................................

  associate( Ms => eop % w,  &
             Ds => eop % D,  &
             dx => eop % dx, &
             mu => eop % mu  )


    ! auxiliaries .............................................................

    ! face mass matrices  . . . . . . . . . . . . . . . . . . . . . . . . . . .

    allocate(M1(0:po,0:po), M2(0:po,0:po), M3(0:po,0:po))

    c(1) = dx(2) * dx(3) / 4
    c(2) = dx(3) * dx(1) / 4
    c(3) = dx(1) * dx(2) / 4

    do j = 0, po
    do i = 0, po
      M1(i,j) = c(1) * Ms(i) * Ms(j)
      M2(i,j) = c(2) * Ms(i) * Ms(j)
      M3(i,j) = c(3) * Ms(i) * Ms(j)
    end do
    end do

    ! delta function . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    allocate(delta_0(0:po), source = ZERO)
    delta_0(0) = ONE

    allocate(delta_P(0:po), source = ZERO)
    delta_P(po) = ONE

    ! add fluxes ...............................................................

    c = ONE / dx

    ! *** OMP/ACC: parallel

    !$omp do collapse(3)
    do n = 1, ne(3)
    do m = 1, ne(2)
    do l = 1, ne(1)

      do k = 0, po
      do j = 0, po
      do i = 0, po

        v(i,j,k,l,m,n)                                                         &

        = v(i,j,k,l,m,n)                                                       &

        + nu *                                                                 &
          (                                                                    &
            M1(j,k) * (                                                        &
                      - (c(1)*Ds( 0,i) + mu(1)*delta_0(i)) * J1_u(j,k,l-1,m,n) &
                      +                        delta_0(i)  * D1_u(j,k,l-1,m,n) &
                      - (c(1)*Ds(po,i) - mu(1)*delta_P(i)) * J1_u(j,k,l  ,m,n) &
                      -                        delta_P(i)  * D1_u(j,k,l  ,m,n) &
                      )                                                        &

          + M2(i,k) * (                                                        &
                      - (c(2)*Ds( 0,j) + mu(2)*delta_0(j)) * J2_u(i,k,l,m-1,n) &
                      +                        delta_0(j)  * D2_u(i,k,l,m-1,n) &
                      - (c(2)*Ds(po,j) - mu(2)*delta_P(j)) * J2_u(i,k,l,m  ,n) &
                      -                        delta_P(j)  * D2_u(i,k,l,m  ,n) &
                      )                                                        &

          + M3(i,j) * (                                                        &
                      - (c(3)*Ds( 0,k) + mu(3)*delta_0(k)) * J3_u(i,j,l,m,n-1) &
                      +                        delta_0(k)  * D3_u(i,j,l,m,n-1) &
                      - (c(3)*Ds(po,k) - mu(3)*delta_P(k)) * J3_u(i,j,l,m,n  ) &
                      -                        delta_P(k)  * D3_u(i,j,l,m,n  ) &
                      )                                                        &
          )
      end do
      end do
      end do

    end do
    end do
    end do
    !$omp end do

  end associate

end subroutine AddFluxes

!===============================================================================

end module CART__DG_Diffusion_CIS_Operator
