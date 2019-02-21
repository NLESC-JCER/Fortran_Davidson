module matrix_free

  use numeric_kinds, only: dp
  implicit none

  public compute_stx_on_the_fly
  
contains

    function compute_matrix_on_the_fly(i, dim) result (vector)
    !> \param[in] i index of the i-th column
    !> \param[in] dim dimension of the resulting column
    !> \return the i-th column of a square matrix of dimension dim
    
    integer, intent(in) :: i, dim
    real(dp), dimension(dim) :: vector
    
    ! local variable
    integer :: j
    real(dp) :: x, y
    
    x = exp(real(i)/real(dim))
    
    do j=1,dim
       y = exp(real(j)/real(dim))
       if (j >= i) then
          vector(j) = cos(log(sqrt(atan2(x, y)))) * 1e-4
       else
          vector(j) = cos(log(sqrt(atan2(y, x)))) * 1e-4
       endif
    end do
    
    vector(i) = vector(i) + real(i)
    
  end function compute_matrix_on_the_fly

  function compute_stx_on_the_fly(i, dim) result (vector)
    !> \param[in] i index of the i-th column
    !> \param[in] dim dimension of the resulting column
    !> \return the i-th column of a square matrix of dimension dim

    integer, intent(in) :: i, dim
    real(dp), dimension(dim) :: vector

    vector = 0d0
    vector(i) = 1d0
    
  end function compute_stx_on_the_fly

end module matrix_free
  
program main
  use numeric_kinds, only: dp
  use davidson, only: generalized_eigensolver,  generate_diagonal_dominant, norm
  use matrix_free, only: compute_matrix_on_the_fly, compute_stx_on_the_fly

  implicit none

  real(dp), dimension(3) :: eigenvalues_DPR
  real(dp), dimension(100, 3) :: eigenvectors_DPR
  real(dp), dimension(100) :: xs, zs
  real(dp), dimension(10, 10) ::  arr
  integer :: iter_i, j

  do j=1, 10
     arr(:, j) = compute_matrix_on_the_fly(j, 10)
  end do

  ! call generalized_eigensolver(compute_matrix_on_the_fly, eigenvalues_DPR, eigenvectors_DPR, 3, "DPR", 1000, &
  !      1d-8, iter_i, 20, compute_stx_on_the_fly)

  ! print *, "eigenvalues: ", eigenvalues_DPR
  
end program main
 
