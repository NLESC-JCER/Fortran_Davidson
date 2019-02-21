program main
  
  use numeric_kinds, only: dp
  use davidson, only: generalized_eigensolver, norm, lapack_generalized_eigensolver, generate_diagonal_dominant
  use test_utils, only: compute_matrix_on_the_fly, compute_stx_on_the_fly, diagonal, write_matrix

  implicit none

  real(dp), dimension(3) :: eigenvalues_DPR
  real(dp), dimension(50, 3) :: eigenvectors_DPR
  real(dp), dimension(50, 50) :: mtx, stx
  real(dp), dimension(50) :: xs, zs
  integer :: iter_i, j

  ! Matrix to check the algorithm
  mtx = generate_diagonal_dominant(50, 1d-3)
  stx = generate_diagonal_dominant(50, 1d-3, 1d0)
  call write_matrix("matrix_free.txt", mtx)
  call write_matrix("stx_free.txt", stx)

  ! NOTE:
  ! compute_matrix_on_the_fly and compute_stx_on_the_fly call some global variables hardcoded just for testing
  
  call generalized_eigensolver(compute_matrix_on_the_fly
  , eigenvalues_DPR, eigenvectors_DPR, 3, "DPR", 1000, &
       1d-8, iter_i, 20, compute_stx_on_the_fly)

  print *, "eigenvalues: ", eigenvalues_DPR
  print *, "Test 1"
  print *, "Check that eigenvalue equation:  H V = l B V holds"
  print *, "DPR method:"
  do j=1,3
     xs = matmul(mtx, eigenvectors_DPR(:, j)) - (eigenvalues_DPR(j) * matmul(stx, eigenvectors_DPR(:, j)))
     print *, "error: ", norm(xs)
     print *, "eigenvalue ", j, ": ", eigenvalues_DPR(j), " succeeded: ", norm(xs) < 1d-8
  end do
  
  print *, "Test 2"
  print *, "If V are the eigenvector then V * V^T = I"
  zs = diagonal(matmul(eigenvectors_DPR, transpose(eigenvectors_DPR)))
  ! There are only 3 eigenvectors
  print *, "DPR method: ", norm(zs(:3)) < sqrt(3.d0)


end program main
