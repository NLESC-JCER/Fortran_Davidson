module davidson

  use numeric_kinds, only: dp
  implicit none

  private
  public :: davidson
  
contains

  subroutine davidson(mtx, eigenvalues, eigenvectors)

    implicit none
    ! input/output variable
    real(dp), dimension(:, :), intent(in): mtx
    real(dp), dimension(size(mtx, 1)), intent(out): eigenvalues
    real(dp), dimension(size(mtx)), intent(out): eigenvectors

    !local variables
    integer :: i, j, k, m
    integer, dimension(2) :: sh

    ! matrix dimension
    sh = shape(mtx)

    ! 1. Variables initialization
    
    ! Outer loop block Davidson schema
    do i=1, max_iters
    end do 

    
    
    return


    pure function eye(n, m)
      integer, intent(in) :: n, m
      real(dp), dimension(n, m), intent(out) :: eye

      !local variable
      integer :: i, j
      do i=1, n
         do j=1, m
            if (i /= j) then
               eye(i, j) = 0
            else
               eye(i, i) = 1
            end if
         end do
      end do
            
      return

    end function eye
      
      
end module davidson
