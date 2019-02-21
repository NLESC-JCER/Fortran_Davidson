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

  integer, parameter :: dim = 1000
  integer, parameter :: lowest = 3
  real(dp), dimension(lowest) :: eigenvalues_DPR
  real(dp), dimension(dim, lowest) :: eigenvectors_DPR
  real(dp), dimension(dim) :: xs
  real(dp), dimension(dim, dim) :: mtx, stx
  integer :: iter_i, j

  do j=1, dim
     mtx(:, j) = compute_matrix_on_the_fly(j, dim)
     stx(:, j) = compute_stx_on_the_fly(j, dim)
  end do

  call generalized_eigensolver(compute_matrix_on_the_fly, eigenvalues_DPR, eigenvectors_DPR, &
       lowest, "DPR", 1000, 1d-8, iter_i, 20, compute_stx_on_the_fly)

  do j=1,lowest
     xs = matmul(mtx, eigenvectors_DPR(:, j)) - (eigenvalues_DPR(j) * matmul(stx, eigenvectors_DPR(:, j)))
     print *, "error: ", norm(xs)
     print *, "eigenvalue ", j, ": ", eigenvalues_DPR(j), " succeeded: ", norm(xs) < 1d-8
  end do

  
end program main
 
