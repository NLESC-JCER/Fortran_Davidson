module test_utils
  use numeric_kinds, only: dp
  implicit none

contains

  subroutine write_vector(path_file, vector)
    !> Write vector to path_file
    character(len=*), intent(in) :: path_file
    real(dp), dimension(:) :: vector

    open(unit=314, file=path_file, status="REPLACE")
    write(314, *) vector
    close(314)
    
  end subroutine write_vector

  function cast_to_double(vector) result(rs)
    !> convert an array of int into double
    integer, dimension(:), intent(in) :: vector
    real(dp), dimension(size(vector))  :: rs
    integer :: i
    do i=1,size(vector)
       rs(i) = real(vector(i))
    end do

  end function cast_to_double
  
end module test_utils

program main
  use numeric_kinds, only: dp
  use davidson, only: generalized_eigensolver
  use lapack_wrapper, only: lapack_generalized_eigensolver
  use array_utils, only: norm, generate_diagonal_dominant
  use test_utils, only: write_vector, cast_to_double
  use benchmark, only: compute_benchmark

  implicit none

  integer, parameter :: dim = 50
  real(dp), dimension(3) :: eigenvalues_DPR, eigenvalues_GJD
  real(dp), dimension(dim, 3) :: eigenvectors_DPR, eigenvectors_GJD
  real(dp), dimension(dim, dim) :: mtx
  real(dp), dimension(dim, dim) :: stx
  real(dp) :: test_norm_eigenvalues
  real(dp), dimension(dim) :: xs
  integer :: iter_i, j

  mtx = generate_diagonal_dominant(dim, 1d-4)
  stx = generate_diagonal_dominant(dim, 1d-4, 1d0)

  call generalized_eigensolver(mtx, eigenvalues_GJD, eigenvectors_GJD, 3, "GJD", 1000, 1d-8, iter_i, 10, stx)
  call generalized_eigensolver(mtx, eigenvalues_DPR, eigenvectors_DPR, 3, "DPR", 1000, 1d-8, iter_i, 10, stx)

  print *, "Test 1"
  test_norm_eigenvalues = norm(eigenvalues_GJD - eigenvalues_DPR)
  print *, "Check that eigenvalues norm computed by different methods are the same: ", test_norm_eigenvalues < 1e-6
  
  print *, "Test 2"
  print *, "Check that eigenvalue equation:  H V = l S V "
  print *, "DPR method:"
  do j=1,3
     xs = matmul(mtx, eigenvectors_DPR(:, j)) - (eigenvalues_DPR(j) * matmul(stx, eigenvectors_DPR(:, j)))
     print *, "eigenvalue ", j, ": ", eigenvalues_DPR(j), " : ", norm(xs) , " : ", norm(xs)< 1.d-7
  end do
  print *, "GJD method:"
  do j=1,3
     xs = matmul(mtx, eigenvectors_GJD(:, j)) - (eigenvalues_GJD(j) * matmul( stx, eigenvectors_GJD(:, j)))
     print *, "eigenvalue ", j, ": ",eigenvalues_GJD(j), " : ",  norm(xs), " : ", norm(xs)< 1.d-7
  end do  
  
end program main
