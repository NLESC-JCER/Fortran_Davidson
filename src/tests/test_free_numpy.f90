program main
  use numeric_kinds, only: dp
  use davidson, only: generalized_eigensolver, generate_diagonal_dominant
  use test_utils, only: apply_mtx_to_vect, apply_stx_to_vect, compute_matrix_on_the_fly, compute_stx_on_the_fly, write_matrix, write_vector

  implicit none

  integer, parameter :: dim = 50
  integer, parameter :: lowest = 3
  real(dp), dimension(3) :: eigenvalues_DPR
  real(dp), dimension(dim, 3) :: eigenvectors_DPR
  real(dp), dimension(dim, dim) :: mtx, stx
  integer :: iter_i, j

  do j=1, dim
     mtx(:, j) = compute_matrix_on_the_fly(j, dim)
     stx(:, j) = compute_stx_on_the_fly(j, dim)
  end do

  ! Write matrices down to test the eigenvalues against numpy
  call write_matrix("matrix_free.txt", mtx)
  call write_matrix("stx_free.txt", stx)

  call generalized_eigensolver(apply_mtx_to_vect, eigenvalues_DPR, eigenvectors_DPR, lowest, &
       "DPR", 1000, 1d-8, iter_i, 20, apply_stx_to_vect)

  call write_vector("eigenvalues_DPR_free.txt",eigenvalues_DPR)
  call write_matrix("eigenvectors_DPR_free.txt", eigenvectors_DPR)

end program main
