program main
  
  use numeric_kinds, only: dp
  use davidson, only: generalized_eigensolver
  use test_utils, only: apply_mtx_to_vect, apply_stx_to_vect, compute_matrix_on_the_fly, compute_stx_on_the_fly, write_matrix
  use array_utils, only: diagonal, norm, generate_diagonal_dominant

  implicit none

  integer, parameter :: dim = 50
  real(dp), dimension(3) :: eigenvalues_DPR
  real(dp), dimension(dim, 3) :: eigenvectors_DPR
  real(dp), dimension(dim) :: xs, zs
  real(dp), dimension(dim, dim) :: mtx, stx
  integer :: iter_i, j

  do j=1, dim
     mtx(:, j) = compute_matrix_on_the_fly(j, dim)
     stx(:, j) = compute_stx_on_the_fly(j, dim)
  end do

  
  call generalized_eigensolver(apply_mtx_to_vect, eigenvalues_DPR, eigenvectors_DPR, 3, "DPR", 1000, &
       1d-8, iter_i, 20, apply_stx_to_vect)

  print "(a, 3f8.4)", "eigenvalues: ", eigenvalues_DPR
  print *, "Test 1"
  print *, "Check that eigenvalue equation:  H V = l B V holds"
  print *, "DPR method:"
  do j=1,3
     xs = matmul(mtx, eigenvectors_DPR(:, j)) - (eigenvalues_DPR(j) * matmul(stx, eigenvectors_DPR(:, j)))
     print "(a, e10.3)", "error: ", norm(xs)
     print "(a, i2, a, e12.5, a, l)", "eigenvalue ", j, ": ", eigenvalues_DPR(j), " succeeded: ", norm(xs) < 1d-8
  end do
  
  print *, "Test 2"
  print *, "If V are the eigenvector then V * V^T = I"
  zs = diagonal(matmul(eigenvectors_DPR, transpose(eigenvectors_DPR)))
  ! There are only 3 eigenvectors
  print "(a, l)", "DPR method: ", norm(zs(:3)) < sqrt(3.d0)


end program main
