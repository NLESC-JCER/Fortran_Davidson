module array_utils

  use numeric_kinds, only: dp
  use lapack_wrapper, only: lapack_generalized_eigensolver, lapack_matmul, lapack_matrix_vector, &
       lapack_qr, lapack_solver, lapack_sort
  implicit none

  !> \private
  private
  !> \public
  public :: concatenate, diagonal,eye, generate_diagonal_dominant, norm, &
       sort_symmetric_matrix

contains

    pure function eye(m, n, alpha)
    !> Create a matrix with ones in the diagonal and zero everywhere else
    !> \param m: number of rows
    !> \param n: number of colums
    !> \param alpha: optional diagonal value
    !> \return matrix of size n x m
    integer, intent(in) :: n, m
    real(dp), dimension(m, n) :: eye
    real(dp), intent(in), optional :: alpha
    
    !local variable
    integer :: i, j
    real(dp) :: x

    ! check optional values
    x = 1.d0
    if (present(alpha)) x = alpha
    
    do i=1, m
       do j=1, n
          if (i /= j) then
             eye(i, j) = 0
          else
             eye(i, i) = x
          end if
       end do
    end do
    
  end function eye

  pure function norm(vector)
    !> compute the norm-2 of a vector
    real(dp), dimension(:), intent(in) :: vector
    real(dp) :: norm

    norm = sqrt(sum(vector ** 2.d0))

  end function norm
  
  subroutine concatenate(arr, brr)

    !> Concatenate two matrices
    !> \param arr: first array
    !> \param brr: second array
    !> \return arr concatenate brr (overwrites arr)

    real(dp), dimension(:, :), intent(inout), allocatable :: arr
    real(dp), dimension(:, :), intent(in) :: brr
    real(dp), dimension(:, :), allocatable :: tmp_array
    integer :: new_dim, dim_cols, dim_rows
    
    ! dimension
    dim_rows = size(arr, 1)
    dim_cols = size(arr, 2)

    ! Number of columns of the new matrix
    new_dim = dim_cols + size(brr, 2)

    ! move to temporal array
    allocate(tmp_array(dim_rows, new_dim))
    tmp_array(:, :dim_cols) = arr
   
    ! Move to new expanded matrix
    deallocate(arr)
    call move_alloc(tmp_array, arr)

    arr(:, dim_cols + 1:) = brr

  end subroutine concatenate
    
  function generate_diagonal_dominant(m, sparsity, diag_val) result(arr)
    !> Generate a diagonal dominant square matrix of dimension m
    !> \param m dimension of the matrix
    !> \param sparsity magnitude order of the off-diagonal values
      
    integer, intent(in) :: m ! size of the square matrix
    real(dp), optional :: diag_val
    integer :: i, j
    real(dp) :: sparsity 
    real(dp), dimension(m, m) :: arr
    call random_number(arr)

    arr = arr * sparsity
    do j=1, m
       do i=1, m
          if (i > j) then
             arr(i, j) = arr(j, i)
          else if(i == j) then
            if (present(diag_val))then
              arr(i,i) = diag_val
            else
             arr(i, i) = i
            end if
          end if
       end do
    end do

  end function generate_diagonal_dominant

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

  subroutine sort_symmetric_matrix(matrix)
    !> Sort symmetric matrix base on the diagonal values

    real(dp), dimension(:, :), intent(inout) :: matrix

    ! local variables
    integer :: i
    real(dp), dimension(size(matrix, 1), size(matrix, 2)) ::  copy
    real(dp), dimension(size(matrix, 1)) ::  d
    integer, dimension(size(matrix, 1)) :: keys

    ! sort diagonal
    d = diagonal(matrix)
    keys = lapack_sort('I', d)

    ! reorder matrix
    copy = matrix

    do i=1,size(d)
       matrix(i, :) = copy(keys(i), :)
    end do

    copy = matrix
    do i=1,size(d)
       matrix(:, i) = copy(:, keys(i))
    end do

  end subroutine sort_symmetric_matrix

end module array_utils
