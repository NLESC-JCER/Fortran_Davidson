program main
  use numeric_kinds, only: dp
  use davidson, only: lapack_eigensolver
  
  implicit none
  integer :: i, j, info
  real(dp) :: arr(5, 2) = reshape([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], [5, 2])
  real(dp) :: matrix(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])
  real(dp), dimension(3) :: eigenvalues
  real(dp), dimension(3, 3) :: eigenvectors
  
  call lapack_eigensolver(matrix, eigenvalues, eigenvectors)

  print *, "info: ", info
  print *, "eigenvalues: ",  eigenvalues
  print *, "eigenvectors: ", eigenvectors(:, 1)

end program main
