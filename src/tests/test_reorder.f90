program main
  
  use array_utils, only: sort_symmetric_matrix, diagonal, norm
  use davidson, only: generalized_eigensolver
  use numeric_kinds, only: dp
  use test_utils, only: read_matrix


  implicit none
 
  integer, parameter :: dim = 864
  integer, parameter :: lowest = 6
  real(dp), dimension(lowest) :: eigenvalues_DPR
  real(dp), dimension(dim, lowest) :: eigenvectors_DPR
  real(dp), dimension(dim, dim) :: mtx
  integer :: iter_i

  ! mtx = read_matrix("/home/felipe/Primer/Fortran_Davidson/src/tests/data/bse_singlet.dat", dim)
  mtx = read_matrix("/home/felipe/Primer/Fortran_Davidson/src/tests/data/sorted_bse.dat", dim)

  ! call sort_symmetric_matrix(mtx)

  call generalized_eigensolver(mtx, eigenvalues_DPR, eigenvectors_DPR, lowest, "DPR", 15, 1d-4, iter_i, lowest*3)

  print *, eigenvalues_DPR
  
end program main
