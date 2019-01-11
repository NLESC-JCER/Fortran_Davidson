program main
  use numeric_kinds, only: dp
  use davidson, only: lapack_eigensolver, lapack_qr, eye
  
  implicit none
  integer :: i, j, info
  real(dp) :: arr(2, 5) = reshape([(i, i=0,9)], [2, 5])
  real(dp) :: matrix(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])
  real(dp), dimension(3) :: eigenvalues
  real(dp), dimension(3, 3) :: eigenvectors, rs, copy
  real(dp), dimension(2, 3) :: xs
  
  copy = matrix
  call lapack_eigensolver(matrix, eigenvalues, eigenvectors)

  print *, "info: ", info
  print *, "eigenvalues: ",  eigenvalues
  print *, "eigenvectors: ", eigenvectors

  rs =  matmul(transpose(eigenvectors), matmul(copy, eigenvectors))
  print *, "check: ", rs(1,1) -eigenvalues(1) < 1e-8
  
  call lapack_qr(arr)

  print *, "QR: ", arr(:, 1)

  xs = eye(3,2)

  print *, "eye: ", xs(:, 1)
  
  
end program main
