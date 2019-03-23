program main
  
  use array_utils, only: sort_symmetric_matrix, diagonal, norm
  use davidson, only: generalized_eigensolver
  use numeric_kinds, only: dp
  use test_utils, only: read_matrix


  implicit none
 
  integer, parameter :: dim = 864
  real(dp), dimension(3) :: eigenvalues_DPR
  real(dp), dimension(dim, 3) :: eigenvectors_DPR
  real(dp), dimension(dim, dim) :: mtx
  integer :: iter_i
  
  mtx = read_matrix("/home/felipe/Primer/Fortran_Davidson/src/tests/data/sorted_bse.txt", dim)

  ! call sort_symmetric_matrix(mtx)

  print *, sum((mtx - transpose(mtx)) ** 2)
  call generalized_eigensolver(mtx, eigenvalues_DPR, eigenvectors_DPR, 3, "DPR", 100, 1d-5, iter_i, 20)

  print *, eigenvalues_DPR
  
end program main
