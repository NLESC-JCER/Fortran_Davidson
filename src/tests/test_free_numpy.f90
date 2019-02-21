program main
  use numeric_kinds, only: dp
  use davidson, only: generalized_eigensolver, generate_diagonal_dominant
  use test_utils, only: compute_matrix_on_the_fly, compute_stx_on_the_fly, write_matrix, write_vector
  use global_variables, only: global_matrix, global_stx

  implicit none

  real(dp), dimension(3) :: eigenvalues_DPR
  real(dp), dimension(50, 3) :: eigenvectors_DPR
  real(dp), dimension(50, 50) :: mtx, stx
  integer :: iter_i
  
  call write_matrix("matrix_free.txt", global_matrix)
  call write_matrix("stx_free.txt", global_stx)

  ! NOTE:
  ! compute_matrix_on_the_fly and compute_stx_on_the_fly call some global variables hardcoded just for testing
  
  ! call eigenvalue solver
  call generalized_eigensolver(compute_matrix_on_the_fly, eigenvalues_DPR, eigenvectors_DPR, 3, "DPR", 1000, &
       1d-8, iter_i, 20, compute_stx_on_the_fly)

  call write_vector("eigenvalues_DPR_free.txt",eigenvalues_DPR)
  call write_matrix("eigenvectors_DPR_free.txt", eigenvectors_DPR)

end program main
