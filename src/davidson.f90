!> \namespace Davidson eigensolver
!> \author Felipe Zapata
module davidson

  use numeric_kinds, only: dp
  use sort, only: argsort
  ! use F90_LAPACK, only: dgeev
  implicit none

  !> \private
  private :: eye, lapack_eigensolver
  !> \public
  public :: eigensolver
  
  interface
     module function compute_correction(mtx, V, eigenvalues, eigenvectors, method) &
          result(correction)
       real(dp), dimension(:), intent(in) :: eigenvalues
       real(dp), dimension(:, :), intent(in) :: mtx, V, eigenvectors
       character(len=10), optional :: method
       real(dp), dimension(size(mtx, 1)) :: correction
       
     end function compute_correction
     
  end interface
  
contains
  
  subroutine eigensolver(mtx, eigenvalues, eigenvectors, lowest, method, max_iters, &
    tolerance)
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
    real(dp), dimension(:), intent(out) :: eigenvalues
    real(dp), dimension(:, :), intent(out) :: eigenvectors
    real(dp), optional :: tolerance
    character(len=10), optional :: method
    
    !local variables
    integer :: k, m
    integer, dimension(2) :: sh
    real(dp) :: residue

    ! Basis of subspace of approximants
    real(dp), dimension(lowest):: guess_eigenvalues
    real(dp), dimension(size(mtx, 1)) :: correction
    real(dp), dimension(:, :), allocatable :: projected, V

    ! Check optional arguments
    if (.not. present(max_iters)) max_iters=1000
    if (.not. present(method)) method="DPR"
    if (.not. present(tolerance)) tolerance=1e-8
        
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

    ! 3. compute the `lowest` eigenvalues and their corresponding eigenvectors
    ! of the project matrix using lapack
    call lapack_eigensolver(projected, eigenvalues, eigenvectors)
    
    ! 4. Check for convergence
    residue = sqrt(sum(guess_eigenvalues - eigenvalues(:lowest) ** 2))
    
    if (residue < tolerance) then
       print *, "done!"
    end if
    
    ! 5. Calculate the correction vector
    correction = compute_correction(mtx, V, eigenvalues, eigenvectors, method)

    ! 6
    
    ! Update guess
    guess_eigenvalues = eigenvalues(:lowest)
    ! end do outer_loop

    ! Free memory
    deallocate(V)
    
    ! Select the lowest eigenvalues and their corresponding eigenvectos
    eigenvalues = eigenvalues(:lowest)
    eigenvectors = eigenvectors(:lowest, :)
    
    return

  end subroutine eigensolver

  subroutine lapack_eigensolver(mtx, eigenvalues, eigenvectors)
    !> Call the DGEEV subroutine lapack to compute ALL the eigenvalues
    !> and corresponding eigenvectors of mtx
    !> \param mtx: Matrix to diaogonalize
    !> \param eigenvalues
    !> \param eigenvectors
    ! input/output
    real(dp), dimension(:, :), intent(in) :: mtx
    real(dp), dimension(:), intent(inout) :: eigenvalues
    real(dp), dimension(:, :), intent(inout) :: eigenvectors
    real(dp), dimension(5 * size(eigenvalues)) :: work ! check dgeev documentation

    ! Local variables
    integer :: dim, info
    integer, dimension(size(mtx, 1)) :: indices ! index of the sort eigenvalues
    real(dp), dimension(size(mtx, 1)) :: eigenvalues_im ! imaginary part
    real(dp), dimension(size(mtx, 1), size(mtx, 1)) :: vl ! check dgeev documentation

    ! dimension of the guess space
    dim = size(mtx, 1)
    
    call DGEEV("N", "V", dim, mtx, dim, eigenvalues, eigenvalues_im, vl, 1, &
         eigenvectors, dim, work, size(work), info)

    ! Return the indices of the lowest eigenvalues
    indices = argsort(eigenvalues)

    ! Sort the eigenvalues and eigenvectors of the basis
    eigenvalues = eigenvalues(indices)
    eigenvectors = eigenvectors(indices, :)
    
    return

  end subroutine lapack_eigensolver
  
  pure function eye(m, n)
    !> Create a matrix with ones in the diagonal and zero everywhere else
    !> \param m: number of rows
    !> \param n: number of colums
    !> \return matrix of size n x m
    integer, intent(in) :: n, m
    real(dp), dimension(n, m) :: eye
    
    !local variable
    integer :: i, j
    do i=1, m
       do j=1, n
          if (i /= j) then
             eye(i, j) = 0
          else
             eye(i, i) = 1
          end if
       end do
    end do
    
  end function eye
  
  
end module davidson


submodule (davidson) correction_methods
  !> submodule containing the implementations of different kind
  !> of correction vector for the Davidson's diagonalization

  implicit none
  
contains

  module function compute_correction(mtx, V, eigenvalues, eigenvectors, method) &
       result(correction)
    real(dp), dimension(:), intent(in) :: eigenvalues
    real(dp), dimension(:, :), intent(in) :: mtx, V, eigenvectors
    character(len=10), optional :: method
    real(dp), dimension(size(mtx, 1)) :: correction
    
    select case (method)
    case ("DPR")
       correction = compute_DPR(mtx, V, eigenvalues, eigenvectors)
    case default
       correction = compute_DPR(mtx, V, eigenvalues, eigenvectors)
    end select
    
  end function compute_correction

  function compute_DPR(mtx, V, eigenvalues, eigenvectors) result(correction)
    real(dp), dimension(:), intent(in) :: eigenvalues
    real(dp), dimension(:, :), intent(in) :: mtx, V, eigenvectors
    real(dp), dimension(size(mtx, 1)) :: correction

  end function compute_DPR
  
     
end submodule correction_methods
