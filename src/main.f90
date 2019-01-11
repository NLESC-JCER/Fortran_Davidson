program main
  use numeric_kinds, only: dp
  use davidson, only: lapack_eigensolver, lapack_qr, eigensolver, eye, generate_triangular

  implicit none
  integer :: i
  real(dp) :: arr(2, 5) = reshape([(i, i=0,9)], [2, 5])
  real(dp) :: matrix(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])
  real(dp), dimension(2) :: eigenvalues
  real(dp), dimension(3, 2) :: eigenvectors
  real(dp), dimension(3, 2) :: copy
  ! real(dp), dimension(10, 10) :: mtx, rs
  ! rs = generate_triangular(10)
  ! mtx = eye(10, 10) + rs * 0.01

  ! call eigensolver(mtx, eigenvalues, eigenvectors, 3, "DPR", 1000, 1e-8)

  ! copy = matrix
  call lapack_eigensolver(matrix, eigenvalues, eigenvectors)
  print *, "eigenvectors: ", eigenvectors

  print *, eigenvalues
end program main
