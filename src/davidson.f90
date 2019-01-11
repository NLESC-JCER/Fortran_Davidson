!> \namespace Davidson eigensolver
!> \author Felipe Zapata
module davidson

  use numeric_kinds, only: dp
  use sort, only: argsort
  implicit none

  !> \private
  private :: concatenate
  !> \public
  public :: eigensolver, eye, lapack_eigensolver, lapack_qr, generate_triangular
  
  interface
     module function compute_correction(mtx, V, eigenvalues, ritz_vectors, lowest, method) &
          result(correction)
       !> compute the correction vector using `method`
       
       integer, intent(in) :: lowest
       real(dp), dimension(:), intent(in) :: eigenvalues
       real(dp), dimension(:, :), intent(in) :: mtx, V, ritz_vectors
       character(len=10), optional :: method
       real(dp), dimension(lowest * 2, size(mtx, 2)) :: correction
       
     end function compute_correction
     
  end interface
  
contains
  
  subroutine eigensolver(mtx, eigenvalues, ritz_vectors, lowest, method, max_iters, &
       tolerance)
    !> The current implementation uses a general  davidson algorithm, meaning
    !> that it compute all the eigenvalues simultaneusly. The family
    !> of Davidson algorithm only differ in the way that the correction
    !> vector is computed.
    
    !> \param mtx: Matrix to diagonalize
    !> \param eigenvalues: Computed eigenvalues
    !> \param ritz_vectors: approximation to the eigenvectors
    !> \param lowest: Number of lowest eigenvalues/ritz_vectors to compute
    !> \param method: Method to compute the correction vector. Available
    !> methods are,
    !> DPR: Diagonal-Preconditioned-Residue
    !> \param max_iters: Maximum number of iterations
    !> \param tolerance: norm-2 error of the eigenvalues
    !> \return eigenvalues and ritz_vectors of the matrix `mtx`


    ! input/output variable
    integer, intent(in) :: lowest
    real(dp), dimension(:, :), intent(in) :: mtx
    real(dp), dimension(:), intent(out) :: eigenvalues
    real(dp), dimension(:, :), intent(out) :: ritz_vectors
    integer, intent(in) :: max_iters
    real(dp), intent(in) :: tolerance
    character(len=10), intent(in) :: method

    
    !local variables
    integer :: dim_sub, max_dim
    real(dp) :: residue

    ! Basis of subspace of approximants
    real(dp), dimension(lowest):: guess_eigenvalues
    real(dp), dimension(:, :), allocatable :: correction, eigenvectors_sub, projected, V

    ! Initial subpsace dimension
    dim_sub = lowest * 2

    ! maximum dimension of the basis for the subspace
    max_dim = size(mtx, 2) / 2
    
    ! 1. Variables initialization
    guess_eigenvalues = 0
    V = eye(size(mtx, 1), dim_sub) ! Initial orthonormal basis

    ! ! Outer loop block Davidson schema
    ! outer_loop: do m=k, max_iters, k

    
    ! 2. Generate subpace matrix problem by projecting into V
    projected = matmul(transpose(V), matmul(mtx, V))

    ! 3. compute the eigenvalues and their corresponding ritz_vectors
    ! for the projected matrix using lapack
    allocate(eigenvectors_sub(dim_sub, size(V,1)))
    call lapack_eigensolver(projected, eigenvalues, eigenvectors_sub)

    ! 4. Check for convergence
    residue = sqrt(sum((guess_eigenvalues - eigenvalues) ** 2))
    print *, "residues: ", residue
    
    if (residue < tolerance) then
       print *, "done!"
    end if
    
    ! 5. Calculate the correction vector
    correction = compute_correction(mtx, V, eigenvalues, ritz_vectors, lowest, method)

    ! 6. Add the correction vectors to the current basis
    if (size(V, 1) <= max_dim) then
       ! append correction to the current basis
       call concatenate(V, correction) 
    else
       ! Reduce the basis of the subspace to the current correction
       deallocate(V)
       call move_alloc(correction, V)
    end if

    ! 7. Orthogonalize basis
    call lapack_qr(V)
    
    ! 8. Update guess
    guess_eigenvalues = eigenvalues
    ! end do outer_loop
    print *, guess_eigenvalues
    
    ! Free memory
    if (allocated(correction)) then
       deallocate(correction)
    end if

    deallocate(projected, V)
    
    ! Select the lowest eigenvalues and their corresponding ritz_vectors
    ! They are sort in increasing order
    eigenvalues = eigenvalues(:lowest)
    ritz_vectors = ritz_vectors(:lowest, :)

    print *, "done!"
    
  end subroutine eigensolver

  subroutine lapack_eigensolver(mtx, eigenvalues, eigenvectors)
    !> Call the DGEEV subroutine lapack to compute ALL the eigenvalues
    !> and corresponding eigenvectors of mtx
    !> \param mtx: Matrix to diaogonalize
    !> \param eigenvalues: lowest eigenvalues
    !> \param eigenvectors: corresponding eigenvectors
    ! input/output
    real(dp), dimension(:, :), intent(in) :: mtx
    real(dp), dimension(:), intent(inout) :: eigenvalues
    real(dp), dimension(:, :), intent(inout) :: eigenvectors
    real(dp), dimension(:), allocatable :: work ! workspace, see lapack documentation

    ! Local variables
    integer :: i, dim, info, lwork, lowest
    integer, dimension(size(mtx, 2)) :: indices ! index of the sort eigenvalues
     ! ALL the eigenvalues of the subpace (re, im)
    real(dp), dimension(size(mtx, 2)) :: eigenvalues_work_re, eigenvalues_work_im
    real(dp), dimension(size(mtx, 2), size(mtx, 2)) :: eigenvectors_work ! All the eigenvectors
    real(dp), dimension(size(mtx, 2), size(mtx, 2)) :: vl ! check dgeev documentation
    
    ! ! dimension of the guess space
    dim = size(mtx, 2)
    lowest = size(eigenvalues) ! Number of lowest eigenvalues to take
    
    ! Query size of the optimal workspace
    allocate(work(1))
    call DGEEV("N", "V", dim, mtx, dim, eigenvalues_work_re, eigenvalues_work_im, vl, 1, &
         eigenvectors_work, dim, work, -1, info)

    ! Allocate memory for the workspace
    lwork = max(1, int(work(1)))
    deallocate(work)
    allocate(work(lwork))

    ! Compute Eigenvalues
    call DGEEV("N", "V", dim, mtx, dim, eigenvalues_work_re, eigenvalues_work_im, vl, 1, &
         eigenvectors_work, dim, work, lwork, info)

    ! Return the indices of the  eigenvalues sorted
    indices = argsort(eigenvalues_work_re)

    ! Sort the eigenvalues and eigenvectors of the basis
    eigenvalues_work_re = eigenvalues_work_re(indices)
    eigenvectors_work = eigenvectors_work(:, indices)
    
    ! Take the lowest eigenvalues/eigenvectors pairs
    eigenvalues = eigenvalues_work_re(1:lowest)
    eigenvectors = eigenvectors_work(:, [(i, i=1,lowest)])
    
    ! release memory
    deallocate(work)
    
  end subroutine lapack_eigensolver


  subroutine lapack_qr(basis)
    !> Orthoghonalize the basis using the QR factorization
    !> \param basis
    !> \return orthogonal basis    
    real(dp), dimension(:, :), intent(inout) :: basis
    real(dp), dimension(:), allocatable :: work ! workspace, see lapack documentation
    real(dp), dimension(:), allocatable :: tau ! see DGEQRF documentation
    integer :: info, lwork, m, n

    ! Matrix shape
    n = size(basis, 1)
    m = size(basis, 2)

    ! 1. Call the QR decomposition
    ! 1.1 Query size of the workspace (Check lapack documentation)
    allocate(work(1))
    call DGEQRF(m, n, basis, max(1, m), tau, work, -1, info)

    ! 1.2 Allocate memory for the workspace
    lwork = max(1, int(work(1)))
    deallocate(work)
    allocate(work(lwork))

    ! 1.3 Call QR factorization
    allocate(tau(min(m, n)))
    call DGEQRF(m, n, basis, max(1, m), tau, work, lwork, info)
    deallocate(work)
    
    ! 2. Generates an orthonormal matrix
    ! 2.1 Query size of the workspace (Check lapack documentation)
    allocate(work(1))
    call DORGQR(m, n, min(m, n), basis, max(1, m), tau, work, -1, info)

    ! 2.2 Allocate memory fo the workspace
    lwork = max(1, int(work(1)))
    deallocate(work)
    allocate(work(lwork))

    ! 2.3 compute the matrix Q
    call DORGQR(m, n, min(m, n), basis, max(1, m), tau, work, lwork, info)
    
    ! release memory
    deallocate(work, tau)
    
  end subroutine lapack_qr

  pure function eye(n, m)
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
             eye(j, i) = 0
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
    integer :: new_dim, dim_cols, dim_rows
    
    ! dimension
    dim_cols = size(arr, 1)
    dim_rows = size(arr, 2)
    ! Number of columns of the new matrix
    new_dim = dim_cols + size(brr, 1)

    ! move to temporal array
    allocate(tmp_array(new_dim, dim_rows))
    tmp_array(:dim_cols, :) = arr

    ! Move to new expanded matrix
    deallocate(arr)
    call move_alloc(tmp_array, arr)

    arr(dim_cols + 1: , :) = brr

  end subroutine concatenate

  function generate_triangular(m) result(trn)
      
    integer, intent(in) :: m ! size of the square matrix
    integer :: i, j
    real(dp), dimension(m, m) :: trn
    call random_number(trn)
    
    do i=1, m
       do j=1, m
          if (i < j) then
             trn(j, i) = trn(i, j)
          else if (i == j) then
             trn(i, i) = i
          end if
       end do
    end do
    
  end function generate_triangular


end module davidson


submodule (davidson) correction_methods
  !> submodule containing the implementations of different kind
  !> of correction vectors for the Davidson's diagonalization algorithm

  implicit none
  
contains

  module function compute_correction(mtx, V, eigenvalues, ritz_vectors, lowest, method) &
       result(correction)
    integer, intent(in) :: lowest
    real(dp), dimension(:), intent(in) :: eigenvalues
    real(dp), dimension(:, :), intent(in) :: mtx, V, ritz_vectors
    character(len=10), optional :: method
    real(dp), dimension(lowest, size(mtx, 2)) :: correction
    
    select case (method)
    case ("DPR")
       correction = compute_DPR(mtx, V, eigenvalues, ritz_vectors, lowest)
    case default
       correction = compute_DPR(mtx, V, eigenvalues, ritz_vectors, lowest)
    end select
    
  end function compute_correction

  function compute_DPR(mtx, V, eigenvalues, ritz_vectors, lowest) result(correction)
    integer, intent(in) :: lowest
    real(dp), dimension(:), intent(in) :: eigenvalues
    real(dp), dimension(:, :), intent(in) :: mtx, V, ritz_vectors
    real(dp), dimension(lowest, size(mtx, 2)) :: correction
    
  end function compute_DPR
  
     
end submodule correction_methods
