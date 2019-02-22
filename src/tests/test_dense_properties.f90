program main
  
  use numeric_kinds, only: dp
  use davidson, only: generalized_eigensolver, norm, generate_diagonal_dominant
  use test_utils, only: diagonal

  implicit none

  integer, parameter :: dim = 50
  integer, parameter :: lowest = 3
  real(dp), dimension(3) :: eigenvalues_DPR, eigenvalues_GJD
  real(dp), dimension(dim, 3) :: eigenvectors_DPR, eigenvectors_GJD
  real(dp), dimension(dim, dim) :: mtx
  real(dp), dimension(dim) :: xs, ys, zs
  real(dp) :: test_norm_eigenvalues
  integer :: iter_i, j

  ! mtx = read_matrix("tests/matrix.txt", 100)
  mtx = generate_diagonal_dominant(dim, 1d-3)

  call generalized_eigensolver(mtx, eigenvalues_GJD, eigenvectors_GJD, lowest, "GJD", 1000, 1d-8, iter_i)
  call generalized_eigensolver(mtx, eigenvalues_DPR, eigenvectors_DPR, lowest, "DPR", 1000, 1d-8, iter_i)

  print *, "Test 1"
  test_norm_eigenvalues = norm(eigenvalues_GJD - eigenvalues_DPR)
  print *, "Check that eigenvalues norm computed by different methods are the same: ", test_norm_eigenvalues < 1e-8
  
  print *, "Test 2"
  print *, "Check that eigenvalue equation:  H V = l V holds"
  print *, "DPR method:"
  do j=1,lowest
     xs = matmul(mtx, eigenvectors_DPR(:, j)) - (eigenvalues_DPR(j) * eigenvectors_DPR(:, j))
     print *, "eigenvalue ", j, ": ", norm(xs) < 1d-8
  end do
  print *, "GJD method:"
  do j=1,lowest
     xs = matmul(mtx, eigenvectors_GJD(:, j)) - (eigenvalues_GJD(j) * eigenvectors_GJD(:, j))
     print *, "eigenvalue ", j, ": ", norm(xs) < 1d-8
  end do

  print *, "Test 3"
  print *, "If V are the eigenvector then V * V^T = I"
  ys = diagonal(matmul(eigenvectors_GJD, transpose(eigenvectors_GJD)))
  zs = diagonal(matmul(eigenvectors_DPR, transpose(eigenvectors_DPR)))
  ! There are only 3 eigenvectors
  print *, "GJD method: ", norm(xs(:3)) < sqrt(real(lowest))
  print *, "DPR method: ", norm(ys(:3)) < sqrt(real(lowest))
  

end program main
