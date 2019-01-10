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
     module function compute_correction(mtx, V, eigenvalues, eigenvectors, lowest, method) &
          result(correction)
       !> compute the correction vector using `method`
       
       integer, intent(in) :: lowest
       real(dp), dimension(:), intent(in) :: eigenvalues
       real(dp), dimension(:, :), intent(in) :: mtx, V, eigenvectors
       character(len=10), optional :: method
       real(dp), dimension(size(mtx, 1), lowest) :: correction
       
     end function compute_correction
     
  end interface
  
contains
  
  subroutine eigensolver(mtx, eigenvalues, eigenvectors, lowest, method, max_iters, &
       tolerance)
    !> The current implementation uses a general  davidson algorithm, meaning
    !> that it compute all the eigenvalues simultaneusly. The family
    !> of Davidson algorithm only differ in the way that the correction
    !> vector is computed.
    
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
    integer :: dim_sub
    integer, dimension(2) :: sh
    real(dp) :: residue

    ! Basis of subspace of approximants
    real(dp), dimension(lowest):: guess_eigenvalues
    real(dp), dimension(size(mtx, 1), lowest * 2) :: correction
    real(dp), dimension(:, :), allocatable :: projected, V

    ! Check optional arguments
    if (.not. present(max_iters)) max_iters=1000
    if (.not. present(method)) method="DPR"
    if (.not. present(tolerance)) tolerance=1e-8
        
    ! subpsace dimension
    dim_sub = lowest * 2
    
    ! matrix dimension
    sh = shape(mtx)

    ! 1. Variables initialization
    guess_eigenvalues = 0
    V = eye(sh(1), dim_sub) ! Initial orthonormal basis
    
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
    correction = compute_correction(mtx, V, eigenvalues, eigenvectors, lowest, method)

    ! 6. Add the correction vectors to the current basis
    call concatenate(V, correction)
    
    ! Update guess
    guess_eigenvalues = eigenvalues(:lowest)
    ! end do outer_loop

    ! Free memory
    deallocate(V, projected)
    
    ! Select the lowest eigenvalues and their corresponding eigenvectos
    eigenvalues = eigenvalues(:lowest)
    eigenvectors = eigenvectors(:lowest, :)
    
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

  subroutine concatenate(arr, brr)
    !> Concatenate two matrices
    !> \param arr: first array
    !> \param brr: second array
    !> \return arr concatenate brr

    real(dp), dimension(:, :), intent(inout), allocatable :: arr
    real(dp), dimension(:, :), intent(in) :: brr
    real(dp), dimension(:, :), allocatable :: tmp_array
    integer :: new_dim, dim_col

    ! dimension
    dim_col = size(arr, 2)
    new_dim = dim_col + size(brr, 2)

    ! move to temporal array
    allocate(tmp_array(size(arr, 1), new_dim))
    tmp_array(: dim_col, :) = arr

    ! Move to new expanded matrix
    deallocate(arr)
    call move_alloc(tmp_array, arr)

    arr(dim_col + 1: , :) = brr

  end subroutine concatenate
  

end module davidson


submodule (davidson) correction_methods
  !> submodule containing the implementations of different kind
  !> of correction vector for the Davidson's diagonalization

  implicit none
  
contains

  module function compute_correction(mtx, V, eigenvalues, eigenvectors, lowest, method) &
       result(correction)
    integer, intent(in) :: lowest
    real(dp), dimension(:), intent(in) :: eigenvalues
    real(dp), dimension(:, :), intent(in) :: mtx, V, eigenvectors
    character(len=10), optional :: method
    real(dp), dimension(size(mtx, 1), lowest) :: correction
    
    select case (method)
    case ("DPR")
       correction = compute_DPR(mtx, V, eigenvalues, eigenvectors, lowest)
    case default
       correction = compute_DPR(mtx, V, eigenvalues, eigenvectors, lowest)
    end select
    
  end function compute_correction

  function compute_DPR(mtx, V, eigenvalues, eigenvectors, lowest) result(correction)
    integer, intent(in) :: lowest
    real(dp), dimension(:), intent(in) :: eigenvalues
    real(dp), dimension(:, :), intent(in) :: mtx, V, eigenvectors
    real(dp), dimension(size(mtx, 1), lowest) :: correction
    
  end function compute_DPR
  
     
end submodule correction_methods
