module test_utils
  use numeric_kinds, only: dp
  implicit none

contains
  function diagonal(matrix)
    !> return the diagonal of a matrix
    real(dp), dimension(:, :), intent(in) :: matrix
    real(dp), dimension(size(matrix, 1)) :: diagonal

    ! local variables
    integer :: i, j, m

    ! dimension of the matrix
    m = size(matrix, 1)
    
    do i=1,m
       do j=1,m
          if  (i == j) then
             diagonal(i) = matrix(i, j)
          end if
       end do
    end do

  end function diagonal

  function read_matrix(path_file, dim) result(mtx)
    !> read a row-major square matrix from a file
    !> \param path_file: path to the file
    !> \param dim: dimension of the square matrix
    !> \return matrix 

    character(len=*), intent(in) :: path_file
    integer, intent(in) :: dim
    real(dp), dimension(dim, dim) :: mtx
    integer :: i
    
    open(unit=3541, file=path_file, status="OLD")
    do i=1,dim
       read(3541, *) mtx(i, :)
    end do
    close(3541)
    
  end function read_matrix

  subroutine write_vector(path_file, vector)
    !> Write vector to path_file
    character(len=*), intent(in) :: path_file
    real(dp), dimension(:) :: vector

    open(unit=314, file=path_file, status="NEW")
    write(314, *) vector
    close(314)
    
  end subroutine write_vector
  
end module test_utils

program main
  use numeric_kinds, only: dp
  use davidson, only: eigensolver, norm, lapack_eigensolver
  use test_utils, only: diagonal , read_matrix, write_vector
  use benchmark, only: compute_benchmark

  implicit none
  
  integer :: i, j
  real(dp), dimension(3) :: eigenvalues_DPR, eigenvalues_GJD
  real(dp), dimension(100, 3) :: eigenvectors_DPR, eigenvectors_GJD
  real(dp), dimension(100, 100) :: mtx
  real(dp) :: test_norm_eigenvalues, test_norm_eigenvectors
  real(dp), dimension(100) :: test_DPR, test_GJD, arr
  real(dp), dimension(7, 2) :: times
  integer, dimension(7) :: dims
  character(len=20) :: arg

  mtx = read_matrix("tests/matrix.txt", 100)

  call eigensolver(mtx, eigenvalues_GJD, eigenvectors_GJD, 3, "GJD", 100, 1d-8)
  call eigensolver(mtx, eigenvalues_DPR, eigenvectors_DPR, 3, "DPR", 100, 1d-8)

  print *, "Test 1"
  test_norm_eigenvalues = norm(eigenvalues_GJD - eigenvalues_DPR)
  print *, "Eigenvalues norm computed by different methods are the same: ", test_norm_eigenvalues < 1e-8

  print *, "Test 2"
  test_DPR = diagonal(matmul(eigenvectors_DPR, transpose(eigenvectors_DPR)))
  test_GJD = diagonal(matmul(eigenvectors_GJD, transpose(eigenvectors_GJD)))
  print *, "If V are the eigenvalues, then V V^T = I"
  print *, "DPR: ", norm(test_DPR(:3) - 1) < 1e-8
  print *, "GJD: ", norm(test_GJD(:3) - 1) < 1e-8

  ! RUn benchmark
  call get_command_argument(1, arg)
  if (arg == "benchmark") then
     print *, "Running Benchmark! "
     dims = [10, 50, 100, 500, 1000, 5000, 10 ** 4] !, 10 ** 5, 5 * 10 ** 5, 10 ** 6, 5 * 10 ** 6, 10 ** 7]
     call compute_benchmark(dims, 3, times)
     print *, times(:, 1)
  end if
  
  call write_vector("times_DPR.txt", times(:, 1))
  
end program main
