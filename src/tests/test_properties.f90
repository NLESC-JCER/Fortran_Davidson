program main
  
  use numeric_kinds, only: dp
  use davidson, only: eigensolver, norm, lapack_eigensolver, generate_diagonal_dominant
  use test_utils, only: diagonal

  implicit none

  real(dp), dimension(3) :: eigenvalues_DPR, eigenvalues_GJD
  real(dp), dimension(50, 3) :: eigenvectors_DPR, eigenvectors_GJD
  real(dp), dimension(50, 50) :: mtx
  real(dp) :: test_norm_eigenvalues
  real(dp), dimension(50) :: test_DPR, test_GJD
  integer :: iter_i
  character(len=20) :: arg1

  ! mtx = read_matrix("tests/matrix.txt", 100)
  mtx = generate_diagonal_dominant(50, 1d-3)

  call eigensolver(mtx, eigenvalues_GJD, eigenvectors_GJD, 3, "GJD", 100, 1d-8, iter_i)
  call eigensolver(mtx, eigenvalues_DPR, eigenvectors_DPR, 3, "DPR", 100, 1d-8, iter_i)

  print *, "Test 1"
  test_norm_eigenvalues = norm(eigenvalues_GJD - eigenvalues_DPR)
  print *, test_norm_eigenvalues
  print *, "Eigenvalues norm computed by different methods are the same: ", test_norm_eigenvalues < 1e-6

  print *, "Test 2"
  test_DPR = diagonal(matmul(eigenvectors_DPR, transpose(eigenvectors_DPR)))
  test_GJD = diagonal(matmul(eigenvectors_GJD, transpose(eigenvectors_GJD)))
  print *, "If V are the eigenvalues, then V V^T = I"
  print *, "DPR: ", norm(test_DPR(:3) - 1) < 1e-6
  print *, "GJD: ", norm(test_GJD(:3) - 1) < 1e-6

  print *, eigenvalues_DPR(:3)
  print *, eigenvalues_GJD(:3)
  
  
end program main
