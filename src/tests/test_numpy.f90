program main
  use numeric_kinds, only: dp
  use davidson, only: eigensolver, generate_diagonal_dominant
  use test_utils, only:  write_matrix, write_vector

  implicit none

  real(dp), dimension(3) :: eigenvalues_DPR, eigenvalues_GJD
  real(dp), dimension(50, 3) :: eigenvectors_DPR, eigenvectors_GJD
  real(dp), dimension(50, 50) :: mtx
  integer :: iter_i

  ! mtx = read_matrix("tests/matrix.txt", 100)
  mtx = generate_diagonal_dominant(50, 1d-3)

  call eigensolver(mtx, eigenvalues_GJD, eigenvectors_GJD, 3, "GJD", 1000, 1d-8, iter_i)
  call eigensolver(mtx, eigenvalues_DPR, eigenvectors_DPR, 3, "DPR", 1000, 1d-8, iter_i)

  call write_matrix("matrix.txt", mtx)
  call write_vector("eigenvalues_GJD.txt",eigenvalues_GJD)
  call write_vector("eigenvalues_DPR.txt",eigenvalues_DPR)

  call write_matrix("eigenvectors_GJD.txt", eigenvectors_GJD)
  call write_matrix("eigenvectors_DPR.txt", eigenvectors_DPR)
end program main
