
module sort

  use numeric_kinds, only: dp
  implicit none

  contains
  
  pure function argsort(a) result(b)
    !> Taken from: https://github.com/certik/fortran-utils
    !> Returns the indices that would sort an array.
    !> \param a: Array to compute sorted indices
    
    ! Arguments
    real(dp), intent(in):: a(:)   ! array of numbers
    integer :: b(size(a))         ! indices into the array 'a' that sort it
    !
    ! Example
    ! -------
    !
    ! rargsort([4.1_dp, 2.1_dp, 2.05_dp, -1.5_dp, 4.2_dp]) ! Returns [4, 3, 2, 1, 5]
    
    integer :: N                           ! number of numbers/vectors
    integer :: i,imin                      ! indices: i, i of smallest
    integer :: temp1                       ! temporary
    real(dp) :: temp2
    real(dp) :: a2(size(a))
    a2 = a
    N=size(a)
    do i = 1, N
       b(i) = i
    end do
    do i = 1, N-1
       ! find ith smallest in 'a'
       imin = minloc(a2(i:),1) + i - 1
    ! swap to position i in 'a' and 'b', if not already there
       if (imin /= i) then
          temp2 = a2(i); a2(i) = a2(imin); a2(imin) = temp2
          temp1 = b(i); b(i) = b(imin); b(imin) = temp1
       end if
    end do
  end function argsort

  
end module sort
