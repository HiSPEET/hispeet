!> summary:  Evaluation of  v = [(S₃xS₂xS₁) Λ (S₃xS₂xS₁)ᵀ] u
!> author:   Joerg Stiller
!> date:     2020/05/29
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> The procedures evaluate a spectral tensor-product operator of the form
!>
!>     (S₃ x S₂ x S₁) Λ (S₃ x S₂ x S₁)ᵀ
!>
!> where S₁,S₂,S₃ are column matrices of the eigenvectors to 1D suboperators
!> for directions 1,2,3 and Lambda is the diagonal matrix of corresponding
!> eigenvalues. Given the operand u(1:n1,1:n2,1:n3,1:ne), the result is
!> computed as
!>
!>     v(i,j,k,e) = sum[l=1:n1, m=1:n2, n=1:n3]
!>                    S₃(k,n) * S₂(j,m) * S₁(i,l) * z(l,m,n,e)
!>
!> where
!>
!>     z(i,j,k,e) = sum[l=1:n1, m=1:n2, n=1:n3]
!>                    Λ(i,j,k) * S₃(n,k) * S₂(m,j) * S₁(l,i) * u(l,m,n,e)
!>
!> for each point (i,j,k) in every element e.
!>
!> When compiled with OpenACC and run on a GPU, any call to this procedure
!> must be followed by an "acc wait" before accessing u or v from CPU or
!> other than the default accelerator queue.
!===============================================================================

module TPO__Spectral__3D
  use TPO__Spectral__3D_CI
  use TPO__Spectral__3D_CA
end module TPO__Spectral__3D
