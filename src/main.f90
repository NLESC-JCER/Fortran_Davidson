program main
  use numeric_kinds, only: dp
  
  implicit none
  integer :: i, j, info
  real(dp) :: matrix(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])
  real(dp), dimension(3) :: eigenvalues_re, eigenvalues_im
  real(dp), dimension(3, 3) :: eigenvectors, vl
  real(dp), dimension(3 * 5) :: work
  
  call DGEEV("N", "V", 3, matrix, 3, eigenvalues_re, eigenvalues_im, vl, 1, &
  eigenvectors, 3, work, 5 * 3, info)

  print *, "info: ", info
  print *, "eigenvalues: ",  eigenvalues_re
  print *, "eigenvectors: ", eigenvectors(1, :)

  print *, "square: ", sqrt(sum(eigenvalues_re ** 2))
  
end program main
