program OperatorSpectrum
  use Kind_Parameters, only: RDP
  use CG_Element_Operators_1D
  use Eigenproblems
  implicit none

  type(CG_ElementOperators1D), allocatable :: eop
  real(RDP),    allocatable :: S(:,:), Lambda_L(:)
  real(RDP),    allocatable :: A(:,:)
  complex(RDP), allocatable :: Lambda_A(:)
  real(RDP) :: dx = 2
  integer :: i, nc, nd, po = 15

  eop = CG_ElementOperators1D(po)

  ! diffusion
  nd = po - 1
  allocate(S(nd,nd), Lambda_L(nd))
  call eop % GetEllipticEigensystem(dx, S, Lambda_L)

  print '(/,A,I0,A)', 'po = ', po, ': Λ(L)'
  do i = 1, nd
    print '(2X,ES12.5)', Lambda_L(i)
  end do

  ! convection
  nc = po
  allocate(A(nc,nc), Lambda_A(nc))
  do i = 1, nc
    A(i,:) = eop % D(i,1:)
  end do
  call SolveNonsymmetricEigenproblem(A, Lambda_A)

  print '(/,A,I0,A)', 'po = ', po, ': Λ(A), |Λ(A)|'
  do i = 1, nc
    print '(2X,ES12.5,1X,SP,1X,ES12.5,"i",2X,S,ES12.5)', &
        Lambda_A(i), abs(Lambda_A(i))
  end do

  deallocate(eop, A, S, Lambda_L, Lambda_A)

  ! max |Λ|

  print '(/,A)', 'po, max|Λ(L)|, max|Λ(A)|'
  do po = 2, 32
    nd = po - 1
    nc = po
    allocate(eop, S(nd,nd), Lambda_L(nd), A(nc,nc), Lambda_A(nc))
    eop = CG_ElementOperators1D(po)
    call eop % GetEllipticEigensystem(dx, S, Lambda_L)
    do i = 1, nc
      A(i,:) = eop % D(i,1:)
    end do
    call SolveNonsymmetricEigenproblem(A, Lambda_A)
    print '(I5,2(3X,ES12.5))', &
        po,                    &
        maxval(abs(Lambda_L)), &
        maxval(abs(Lambda_A))
    deallocate(eop, A, S, Lambda_L, Lambda_A)
  end do


end program OperatorSpectrum
