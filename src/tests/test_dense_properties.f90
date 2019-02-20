program main
  
  use numeric_kinds, only: dp
  use davidson, only: generalized_eigensolver, norm, lapack_generalized_eigensolver, generate_diagonal_dominant
  use test_utils, only: diagonal

  implicit none

  real(dp), dimension(3) :: eigenvalues_DPR, eigenvalues_GJD
  real(dp), dimension(50, 3) :: eigenvectors_DPR, eigenvectors_GJD
  real(dp), dimension(50, 50) :: mtx
  real(dp), dimension(50) :: xs, ys, zs
  real(dp) :: test_norm_eigenvalues
  integer :: iter_i, j

  ! mtx = read_matrix("tests/matrix.txt", 100)
  mtx = generate_diagonal_dominant(50, 1d-3)

  call generalized_eigensolver(mtx, eigenvalues_GJD, eigenvectors_GJD, 3, "GJD", 1000, 1d-8, iter_i)
  call generalized_eigensolver(mtx, eigenvalues_DPR, eigenvectors_DPR, 3, "DPR", 1000, 1d-8, iter_i)

  print *, "Test 1"
  test_norm_eigenvalues = norm(eigenvalues_GJD - eigenvalues_DPR)
  print *, "Check that eigenvalues norm computed by different methods are the same: ", test_norm_eigenvalues < 1e-8
  
  print *, "Test 2"
  print *, "Check that eigenvalue equation:  H V = l V holds"
  print *, "DPR method:"
  do j=1,3
     xs = matmul(mtx, eigenvectors_DPR(:, j)) - (eigenvalues_DPR(j) * eigenvectors_DPR(:, j))
     print *, "eigenvalue ", j, ": ", norm(xs) < 1d-8
  end do
  print *, "GJD method:"
  do j=1,3
     xs = matmul(mtx, eigenvectors_GJD(:, j)) - (eigenvalues_GJD(j) * eigenvectors_GJD(:, j))
     print *, "eigenvalue ", j, ": ", norm(xs) < 1d-8
  end do

  print *, "Test 3"
  print *, "If V are the eigenvector then V * V^T = I"
  ys = diagonal(matmul(eigenvectors_GJD, transpose(eigenvectors_GJD)))
  zs = diagonal(matmul(eigenvectors_DPR, transpose(eigenvectors_DPR)))
  ! There are only 3 eigenvectors
  print *, "GJD method: ", norm(xs(:3)) < sqrt(3.d0)
  print *, "DPR method: ", norm(ys(:3)) < sqrt(3.d0)
  

end program main
