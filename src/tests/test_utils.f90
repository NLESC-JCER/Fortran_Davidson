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

