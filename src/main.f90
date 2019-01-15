program main
  use numeric_kinds, only: dp
  use davidson, only: eigensolver, generate_diagonal_dominant, lapack_dgesv, lapack_eigensolver

  implicit none
  integer :: i, j
  real(dp) :: arr(2, 5) = reshape([(i, i=0,9)], [2, 5])
  real(dp) :: brr(2, 2) = reshape([(i, i=10,13)], [2, 2])
  real(dp) :: matrix(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])
  real(dp), dimension(3) :: eigenvalues
  real(dp), dimension(20, 3) :: eigenvectors
  real(dp), dimension(3, 3) :: crr, copy
  real(dp), dimension(20, 20) :: mtx, rs
  real(dp), dimension(3, 1) :: xs

  mtx = generate_diagonal_dominant(20, 1d-4)

  open(unit=3541, file="matrix.txt")
  do i=1,20
     do j=1,20
        if (i<=j) then
           write(3541, *) mtx(i, j)
        end if
     end do
  end do
  close(3541)

  call eigensolver(mtx, eigenvalues, eigenvectors, 3, "DPR", 100, 1d-8)
  print *, eigenvalues
  print *, eigenvectors

end program main
