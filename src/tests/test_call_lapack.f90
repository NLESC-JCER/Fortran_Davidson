program main
  use numeric_kinds, only: dp
  use test_utils, only:  write_matrix, write_vector
  use davidson, only: generate_diagonal_dominant, lapack_generalized_eigensolver
  
  implicit none
  
  real(dp), dimension(50) :: eigenvalues
  real(dp), dimension(50, 50) :: eigenvectors
  real(dp), dimension(50, 50) :: copy, mtx, stx
  integer :: iter_i

  ! mtx = read_matrix("tests/matrix.txt", 100)
  mtx = generate_diagonal_dominant(50, 1d-3)
  copy = mtx
  stx = generate_diagonal_dominant(50, 1d-3)
  call write_matrix("test_lapack_matrix.txt", mtx)  
  call write_matrix("test_lapack_stx.txt", stx)  

  ! standard eigenvalue problem
  call lapack_generalized_eigensolver(copy, eigenvalues, eigenvectors)

  call write_vector("test_lapack_eigenvalues.txt",eigenvalues)
  call write_matrix("test_lapack_eigenvectors.txt",eigenvectors)

  ! General eigenvalue problem
  call lapack_generalized_eigensolver(copy, eigenvalues, eigenvectors, stx)
  call write_vector("test_lapack_eigenvalues_gen.txt",eigenvalues)
  call write_matrix("test_lapack_eigenvectors_gen.txt",eigenvectors)

  
  
end program main
