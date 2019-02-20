module test_utils
  use numeric_kinds, only: dp
  implicit none

contains

  function compute_vector_on_fly(i, dim) result (vector)
    !> \param[in] i index of the i-th column
    !> \param[in] dim dimension of the resulting column
    !> \return the i-th column of a square matrix of dimension dim

    integer, intent(in) :: i, dim
    real(dp), dimension(dim) :: vector

    ! local variables
    real(dp), dimension(dim, dim) :: mtx
    mtx = read_matrix("matrix_free.txt", dim)
    
    vector = mtx(:, i)
    
  end function compute_vector_on_fly


  function compute_vector_generalized_eigenvalue(i, dim) result (vector)
    !> \param[in] i index of the i-th column
    !> \param[in] dim dimension of the resulting column
    !> \return the i-th column of a square matrix of dimension dim

    integer, intent(in) :: i, dim
    real(dp), dimension(dim) :: vector

    ! Generate random vector
    call random_number(vector)
    ! scale the vector
    vector = 0.d0
    ! generate diagonal dominant component
    vector(i) = 1d0
    
  end function compute_vector_generalized_eigenvalue

  
  
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
    real(dp), dimension(:), intent(in) :: vector
    integer :: i

    open(unit=314, file=path_file, status="REPLACE")
    do i=1,size(vector)
       write(314, *) vector(i)
    end do
    close(314)
    
  end subroutine write_vector

  subroutine write_matrix(path_file, mtx)
    !> Write matrix to path_file
    character(len=*), intent(in) :: path_file
    real(dp), dimension(:, :), intent(in) :: mtx
    integer :: i, j

    open(unit=314, file=path_file, status="REPLACE")
    do i=1,size(mtx, 1)
       do j=1,size(mtx, 2)
          write(314, *) mtx(i, j)
       end do
    end do
    close(314)
    
  end subroutine write_matrix

  
end module test_utils

