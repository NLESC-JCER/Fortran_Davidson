!> \namespace Davidson eigensolver
!> \author Felipe Zapata
module davidson


  
  use numeric_kinds, only: dp
  ! use F90_LAPACK, only: dgeev
  implicit none

  !> \private
  private :: eye
  public :: eigensolver
  
contains

  subroutine eigensolver(mtx, eigenvalues, eigenvectors, lowest, method, max_iters)
    !> \param mtx: Matrix to diagonalize
    !> \param lowest: Number of lowest eigenvalues/eigenvectors to compute
    !> \param method: Method to compute the correction vector. Available
    !> methods are,
    !> DPR: Diagonal-Preconditioned-Residue
    !> \param max_iters: Maximum number of iterations
    !> \return eigenvalues and eigenvectors of the matrix `mtx`


    ! input/output variable
    integer, intent(in) :: lowest
    integer, optional :: max_iters
    real(dp), dimension(:, :), intent(in) :: mtx
    real(dp), dimension(size(mtx, 1)), intent(out) :: eigenvalues
    real(dp), dimension(size(mtx, 1), lowest), intent(out) :: eigenvectors
    character(len=10), optional :: method
    
    !local variables
    integer :: k, m
    integer, dimension(2) :: sh

    ! Basis of subspace of approximants
    real(dp), dimension(lowest):: guess_eigenvalues
    real(dp), dimension(:, :), allocatable :: projected, V

    ! Check optional arguments
    if (.not. present(max_iters)) max_iters=1000
    if (.not. present(method)) method="DPR"
        
    ! Initial subpsace
    k = lowest * 2
    
    ! matrix dimension
    sh = shape(mtx)

    ! 1. Variables initialization
    guess_eigenvalues = 0
    allocate(V(sh(1), k))
    V = eye(sh(1), k) ! Initial orthonormal basis
    
    ! ! Outer loop block Davidson schema
    ! outer_loop: do m=k, max_iters, k

    ! 2. Generate subpace matrix problem
    projected = matmul(transpose(V), matmul(mtx, V))

    ! 3. compute the eigenvalues and eigenvectors of the matrix
    ! call lapack_eigensolver(eigenvalues, eigenvectors)
    
    ! 4. Check for convergence

    ! 5. Calculate the correction vector
    
    ! end do outer_loop

    deallocate(V)
    
    return

  end subroutine eigensolver

  pure function eye(n, m)
    !> Create a diagonal matrix
    !> \param n: number of rows
    !> \param m: number of colums
    !> \return matrix of size n x m
    integer, intent(in) :: n, m
    real(dp), dimension(n, m) :: eye
    
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
