program main
  use numeric_kinds, only: dp
  
  implicit none
  integer :: i, j, info
  real(dp) :: arr(5, 2) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6, 0], [5, 2])
  real(dp) :: matrix(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])
  real(dp), dimension(3) :: eigenvalues_re, eigenvalues_im
  real(dp), dimension(3, 3) :: eigenvectors, vl
  real(dp), dimension(3 * 5) :: work
  real(dp), dimension(2, 3) :: xs

  print *, arr(:, 1)
  call DGEEV("N", "V", 3, matrix, 3, eigenvalues_re, eigenvalues_im, vl, 1, &
  eigenvectors, 3, work, 5 * 3, info)

  print *, "info: ", info
  print *, "eigenvalues: ",  eigenvalues_re([3,1,2])
  print *, "eigenvectors: ", eigenvectors(1, :)
  
end program main
