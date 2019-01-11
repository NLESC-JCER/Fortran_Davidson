program main
  use numeric_kinds, only: dp
  use davidson, only: lapack_eigensolver, lapack_qr
  
  implicit none
  integer :: i, j, info
  real(dp) :: arr(5, 2) = reshape([(i, i=0,9)], [5, 2])
  real(dp) :: matrix(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])
  real(dp), dimension(3) :: eigenvalues
  real(dp), dimension(3, 3) :: eigenvectors, rs, copy

  copy = matrix
  call lapack_eigensolver(matrix, eigenvalues, eigenvectors)

  print *, "info: ", info
  print *, "eigenvalues: ",  eigenvalues
  print *, "eigenvectors: ", eigenvectors

  rs =  matmul(transpose(eigenvectors), matmul(copy, eigenvectors))
  print *, "check: ", rs(1,1) -eigenvalues(1) < 1e-8
  
  call lapack_qr(arr)

  print *, "QR: ", arr(:, 1)
  
end program main
