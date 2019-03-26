program main

  use array_utils, only: diagonal, norm
  use davidson, only: generalized_eigensolver
  use numeric_kinds, only: dp
  use test_utils, only: read_matrix

  implicit none

  integer, parameter :: dim = 864
  integer, parameter :: lowest = 6
  real(dp), dimension(lowest) :: eigenvalues_DPR, eigenvalues_GJD, expected
  real(dp), dimension(dim, lowest) :: eigenvectors_DPR, eigenvectors_GJD
  real(dp), dimension(dim, dim) :: mtx
  integer :: iter_i

  mtx = read_matrix("../src/tests/data/bse_singlet.dat", dim)

  print *, "DPR Method"
  call generalized_eigensolver(mtx, eigenvalues_DPR, eigenvectors_DPR, lowest, "DPR", 50, 1d-4, iter_i, lowest*3)

  expected = [0.30445426, 0.31341032, 0.31360998, 0.33246853, 0.34212415, 0.35761287]

  print "(a, i2)", "steps to convergence: ", iter_i
  print "(a, 6f8.4)", "Eigenvalues:", eigenvalues_DPR
  print "(a, e10.3)", "Relative error with expected eigenvalues: ", norm(eigenvalues_DPR - expected)

  print *, "GJD Method"

  call generalized_eigensolver(mtx, eigenvalues_GJD, eigenvectors_GJD, lowest, "GJD", 10, 1d-4, iter_i, lowest*2)

  print "(a, i2)", "steps to convergence: ", iter_i
  print "(a, 6f8.4)", "Eigenvalues:", eigenvalues_GJD
  print "(a, e10.3)", "Relative error with expected eigenvalues: ", norm(eigenvalues_GJD - expected)
  
end program main
