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
  use davidson, only: eigensolver, norm, lapack_eigensolver, generate_diagonal_dominant
  use test_utils, only: diagonal , read_matrix, write_vector, cast_to_double
  use benchmark, only: compute_benchmark

  implicit none

  real(dp), dimension(3) :: eigenvalues_DPR, eigenvalues_GJD
  real(dp), dimension(50, 3) :: eigenvectors_DPR, eigenvectors_GJD
  real(dp), dimension(50, 50) :: mtx
  real(dp) :: test_norm_eigenvalues, sparsity
  real(dp), dimension(50) :: test_DPR, test_GJD
  real(dp), dimension(5, 2) :: times
  integer, dimension(5, 2) :: iters
  integer, dimension(5) :: dims
  integer :: iter_i
  character(len=20) :: arg1

  ! mtx = read_matrix("tests/matrix.txt", 100)
  mtx = generate_diagonal_dominant(50, 1d-3)

  call eigensolver(mtx, eigenvalues_GJD, eigenvectors_GJD, 3, "GJD", 100, 1d-8, iter_i)
  call eigensolver(mtx, eigenvalues_DPR, eigenvectors_DPR, 3, "DPR", 100, 1d-8, iter_i)
  
  print *, "Test 1"
  test_norm_eigenvalues = norm(eigenvalues_GJD - eigenvalues_DPR)
  print *, test_norm_eigenvalues
  print *, "Eigenvalues norm computed by different methods are the same: ", test_norm_eigenvalues < 1e-6

  print *, "Test 2"
  test_DPR = diagonal(matmul(eigenvectors_DPR, transpose(eigenvectors_DPR)))
  test_GJD = diagonal(matmul(eigenvectors_GJD, transpose(eigenvectors_GJD)))
  print *, "If V are the eigenvalues, then V V^T = I"
  print *, "DPR: ", norm(test_DPR(:3) - 1) < 1e-6
  print *, "GJD: ", norm(test_GJD(:3) - 1) < 1e-6

  print *, eigenvalues_DPR(:3)
  print *, eigenvalues_GJD(:3)
  
  ! RUn benchmark
  call get_command_argument(1, arg1)
  if (arg1 == "benchmark") then
     print *, "Running Benchmark! "
     dims = [10, 50, 100, 500, 1000] !, 5000, 7500, 10 ** 4]
     sparsity = 1d-4
     call compute_benchmark(dims, 3, sparsity, times, iters)
  end if
  
  call write_vector("times_DPR.txt", times(:, 1))
  call write_vector("times_GJD.txt", times(:, 2))

  call write_vector("cycles_DPR.txt", cast_to_double(iters(:, 1)))
  call write_vector("cycles_GJD.txt", cast_to_double(iters(:, 2)))
  
end program main
