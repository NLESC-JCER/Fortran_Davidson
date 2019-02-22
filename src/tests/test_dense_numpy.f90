program main
  use numeric_kinds, only: dp
  use davidson, only: generalized_eigensolver, generate_diagonal_dominant
  use test_utils, only:  write_matrix, write_vector

  implicit none

  integer, parameter :: dim = 50
  integer, parameter :: lowest = 3
  real(dp), dimension(3) :: eigenvalues_DPR, eigenvalues_GJD, eigenvalues_DPR_gen, eigenvalues_GJD_gen
  real(dp), dimension(dim, 3) :: eigenvectors_DPR, eigenvectors_GJD, eigenvectors_DPR_gen, eigenvectors_GJD_gen
  real(dp), dimension(dim, dim) ::  mtx, stx
  integer :: iter_i

  ! mtx = read_matrix("tests/matrix.txt", 100)
  mtx = generate_diagonal_dominant(dim, 1d-3)
  
  ! call eigenvalue solver
  call generalized_eigensolver(mtx, eigenvalues_GJD, eigenvectors_GJD, lowest, "GJD", 1000, 1d-8, iter_i)
  call generalized_eigensolver(mtx, eigenvalues_DPR, eigenvectors_DPR, lowest, "DPR", 1000, 1d-8, iter_i)

  call write_matrix("test_dense_spec_matrix.txt", mtx)
  call write_vector("test_dense_spec_eigenvalues_GJD.txt",eigenvalues_GJD)
  call write_vector("test_dense_spec_eigenvalues_DPR.txt",eigenvalues_DPR)
  call write_matrix("test_dense_spec_eigenvectors_GJD.txt", eigenvectors_GJD)
  call write_matrix("test_dense_spec_eigenvectors_DPR.txt", eigenvectors_DPR)

  ! call generalized eigenvalue solver
  stx = generate_diagonal_dominant(dim, 1d-3, 1d0)
  call generalized_eigensolver(mtx, eigenvalues_GJD_gen, eigenvectors_GJD_gen, lowest, "GJD", 1000, 1d-8, iter_i, 10, stx)
  call generalized_eigensolver(mtx, eigenvalues_DPR_gen, eigenvectors_DPR_gen, lowest, "DPR", 1000, 1d-8, iter_i, 10, stx)

  call write_matrix("test_dense_gen_matrix.txt", mtx)
  call write_matrix("test_dense_gen_stx.txt", stx)
  call write_vector("test_dense_gen_eigenvalues_GJD.txt",eigenvalues_GJD_gen)
  call write_vector("test_dense_gen_eigenvalues_DPR.txt",eigenvalues_DPR_gen)
  call write_matrix("test_dense_gen_eigenvectors_GJD.txt", eigenvectors_GJD_gen)
  call write_matrix("test_dense_gen_eigenvectors_DPR.txt", eigenvectors_DPR_gen)

  
end program main
