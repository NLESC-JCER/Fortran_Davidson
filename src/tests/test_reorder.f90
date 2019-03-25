program main

  use array_utils, only: sort_symmetric_matrix, diagonal, norm
  use davidson, only: generalized_eigensolver
  use numeric_kinds, only: dp
  use test_utils, only: read_matrix

  implicit none

  integer, parameter :: dim = 864
  integer, parameter :: lowest = 6
  real(dp), dimension(lowest) :: eigenvalues_DPR, expected
  real(dp), dimension(dim, lowest) :: eigenvectors_DPR
  real(dp), dimension(dim, dim) :: mtx
  integer :: iter_i

  mtx = read_matrix("../src/tests/data/sorted_bse.dat", dim)

  call generalized_eigensolver(mtx, eigenvalues_DPR, eigenvectors_DPR, lowest, "DPR", 15, 1d-4, iter_i, lowest*3)

  expected = [0.30445426, 0.31341032, 0.31360998, 0.33246853, 0.34212415, 0.35761287]

  print *, "Error expected eigenvalues: ", norm(eigenvalues_DPR - expected)

end program main
