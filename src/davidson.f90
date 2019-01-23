!> \namespace Davidson eigensolver
!> \author Felipe Zapata
module davidson

  use numeric_kinds, only: dp
  implicit none

  !> \private
  private :: eye, lapack_solver, lapack_qr, concatenate
  !> \public
  public :: eigensolver, generate_diagonal_dominant, norm, lapack_eigensolver
  
  interface
     module function compute_correction(mtx, V, eigenvalues, eigenvectors, method) &
          result(correction)
       !> compute the correction vector using a given `method` for the Davidson algorithm
       !> See correction_methods submodule for the implementations
       !> \param mtx: Original matrix
       !> \param V: Basis of the iteration subspace
       !> \param eigenvalues: of the reduce problem
       !> \param eigenvectors: of the reduce problem
       !> \param method: name of the method to compute the correction
       
       real(dp), dimension(:), intent(in) :: eigenvalues
       real(dp), dimension(:, :), intent(in) :: mtx, V, eigenvectors
       character(len=*), optional, intent(in) :: method
       real(dp), dimension(size(mtx, 1), size(V, 2)) :: correction
       
       
     end function compute_correction
     
  end interface
  
contains
  
  subroutine eigensolver(mtx, eigenvalues, ritz_vectors, lowest, method, max_iters, &
       tolerance, iters)
    !> The current implementation uses a general  davidson algorithm, meaning
    !> that it compute all the eigenvalues simultaneusly using a block approach.
    !> The family of Davidson algorithm only differ in the way that the correction
    !> vector is computed.
    
    !> \param mtx: Matrix to diagonalize
    !> \param eigenvalues: Computed eigenvalues
    !> \param ritz_vectors: approximation to the eigenvectors
    !> \param lowest: Number of lowest eigenvalues/ritz_vectors to compute
    !> \param method: Method to compute the correction vector. Available
    !> methods are,
    !>    DPR: Diagonal-Preconditioned-Residue
    !>    GJD: Generalized Jacobi Davidson
    !> \param max_iters: Maximum number of iterations
    !> \param tolerance: norm-2 error of the eigenvalues
    !> \param method: Method to compute the correction vectors
    !> \param iters: Number of iterations until convergence
    !> \return eigenvalues and ritz_vectors of the matrix `mtx`
    implicit none
    ! input/output variable
    integer, intent(in) :: lowest
    real(dp), dimension(:, :), intent(in) :: mtx
    real(dp), dimension(lowest), intent(out) :: eigenvalues
    real(dp), dimension(:, :), intent(out) :: ritz_vectors
    integer, intent(in) :: max_iters
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: method
    integer, intent(out) :: iters
    
    !local variables
    integer :: i, j, dim_sub, max_dim
    
    ! Basis of subspace of approximants
    real(dp), dimension(size(mtx, 1)) :: guess, rs
    real(dp), dimension(lowest):: errors

    ! Working arrays
    real(dp), dimension(:), allocatable :: eigenvalues_sub
    real(dp), dimension(:, :), allocatable :: correction, eigenvectors_sub, projected, V

    ! Iteration subpsace dimension
    dim_sub = lowest + (lowest / 2)

    ! maximum dimension of the basis for the subspace
    max_dim = size(mtx, 2) / 2
    
    ! 1. Variables initialization
    V = eye(size(mtx, 1), dim_sub) ! Initial orthonormal basis

    ! ! Outer loop block Davidson schema
    outer_loop: do i=1, max_iters

       ! 2. Generate subpace matrix problem by projecting into V
       projected = lapack_matmul('T', 'N', V, lapack_matmul('N', 'N', mtx, V))
       
       ! 3. compute the eigenvalues and their corresponding ritz_vectors
       ! for the projected matrix using lapack
       call check_deallocate_matrix(eigenvectors_sub)

       if (allocated(eigenvalues_sub)) then
          deallocate(eigenvalues_sub)
       end if

       ! allocate(eigenvectors_sub(size(projected, 1), size(projected,2)))
       allocate(eigenvalues_sub(size(projected, 1)))
       allocate(eigenvectors_sub(size(projected, 1), size(projected, 2)))

       call lapack_eigensolver(projected, eigenvalues_sub, eigenvectors_sub)

       ! 4. Check for convergence
       ritz_vectors =  matmul(V, eigenvectors_sub(:, 1:lowest))
       ! ritz_vectors = lapack_matmul('N', 'N', V, eigenvectors_sub(:, :lowest))
       do j=1,lowest
          guess = eigenvalues_sub(j) * ritz_vectors(:, j)
          rs  = matmul(mtx, ritz_vectors(:, j)) - guess
          ! rs = lapack_matrix_vector('N', mtx, ritz_vectors(:, j)) - guess
          errors(j) = norm(rs)
       end do

       if (all(errors < tolerance)) then
          iters = i
          exit
       end if
       
       ! 5. Add the correction vectors to the current basis
       if (size(V, 2) <= max_dim) then
          ! append correction to the current basis
          call check_deallocate_matrix(correction)
          allocate(correction(size(mtx, 1), size(V, 2)))
          correction = compute_correction(mtx, V, eigenvalues_sub, eigenvectors_sub, method)
          ! 6. Increase Basis size
          call concatenate(V, correction)
          
       else
          ! 6. Otherwise reduce the basis of the subspace to the current correction
          V = matmul(V, eigenvectors_sub(:, :dim_sub))
          ! V = lapack_matmul('N', 'N', V, eigenvectors_sub(:, :dim_sub))
       end if

       ! 7. Orthogonalize basis
       call lapack_qr(V)

    end do outer_loop

    !  8. Check convergence
    if (i > max_iters / dim_sub) then
       print *, "Warning: Algorithm did not converge!!"
    end if

    
    ! Select the lowest eigenvalues and their corresponding ritz_vectors
    ! They are sort in increasing order
    eigenvalues = eigenvalues_sub(:lowest)

    
    ! Free memory
    call check_deallocate_matrix(correction)

    deallocate(eigenvalues_sub, eigenvectors_sub, V)
    
  end subroutine eigensolver

  subroutine lapack_eigensolver(mtx, eigenvalues, eigenvectors)
    !> Call the DSYEV subroutine lapack to compute ALL the eigenvalues
    !> and corresponding eigenvectors of mtx
    !> \param mtx: Matrix to diaogonalize
    !> \param eigenvalues: lowest eigenvalues
    !> \param eigenvectors: corresponding eigenvectors
    !> \return eigenvalues/eigenvectors
    
    ! input/output
    implicit none
    real(dp), dimension(:, :), intent(in) :: mtx
    real(dp), dimension(size(mtx, 1)), intent(inout) :: eigenvalues
    real(dp), dimension(size(mtx, 1), size(mtx, 2)), intent(inout) :: eigenvectors

    ! Local variables
    integer :: dim, info, lwork
     ! ALL the eigenvalues of the subpace (re, im)
    real(dp), dimension(size(mtx, 1)) :: eigenvalues_work
    real(dp), dimension(:), allocatable :: work ! workspace, see lapack documentation
    
    ! ! dimension of the guess space
    dim = size(mtx, 1)
    
    ! Query size of the optimal workspace
    allocate(work(1))
    
    call DSYEV("V", "U", dim, mtx, dim, eigenvalues_work, work, -1, info)
    call check_lapack_call(info, "DSYEV")

    ! Allocate memory for the workspace
    lwork = max(1, int(work(1)))
    deallocate(work)
    allocate(work(lwork))

    ! Compute Eigenvalues
    call DSYEV("V", "U", dim, mtx, dim, eigenvalues_work, work, lwork, info)
    call check_lapack_call(info, "DSYEV")

    ! Sort the eigenvalues and eigenvectors of the basis
    eigenvalues = eigenvalues_work
    eigenvectors = mtx
    
    ! release memory
    deallocate(work)
    
  end subroutine lapack_eigensolver


  subroutine lapack_qr(basis)
    !> Orthoghonalize the basis using the QR factorization.
    !> QR factorization of the M-by-N (M>N) matrx A=Q*R in the form where
    !> Q is square M-by-M matrix and R is an upper triangular M-by-N matrix.
    !> The equality A=Q*R can be re-written also as a product Q1*R1 where Q1
    !> is a rectangular M-by-N submatrix of the matrix Q and R1 is M-by-M
    !> submatrix of the R. Let us note that columns of Q1 are orthonormal
    !> (they are orthogonal to each other and have norms equal to 1).
    !> The equality A=Q1*R1 can be treated as every column of A is a linear
    !> combination of Q1 columns, i.e. they span the same linear space.
    !> In other words, columns of Q1 is the result of ortogonalization of columns A.
    !> DGEQRF does not not compute Q directly, DORGQR must be call subsequently.
    
    !> \param basis
    !> \return orthogonal basis    

    implicit none
    real(dp), dimension(:, :), intent(inout) :: basis
    real(dp), dimension(:), allocatable :: work ! workspace, see lapack documentation
    real(dp), dimension(size(basis, 2)) :: tau ! see DGEQRF documentation
    integer :: info, lwork, m, n

    ! Matrix shape
    m = size(basis, 1)
    n = size(basis, 2)

    ! 1. Call the QR decomposition
    ! 1.1 Query size of the workspace (Check lapack documentation)
    allocate(work(1))
    call DGEQRF(m, n, basis, m, tau, work, -1, info)
    call check_lapack_call(info, "DGEQRF")

    ! 1.2 Allocate memory for the workspace
    lwork = max(1, int(work(1)))
    deallocate(work)
    allocate(work(lwork))

    ! 1.3 Call QR factorization
    call DGEQRF(m, n, basis, m, tau, work, lwork, info)
    call check_lapack_call(info, "DGEQRF")
    deallocate(work)

    ! 2. Generates an orthonormal matrix
    ! 2.1 Query size of the workspace (Check lapack documentation)
    allocate(work(1))
    call DORGQR(m, n, min(m, n), basis, m, tau, work, -1, info)
    call check_lapack_call(info, "DORGQR")

    ! 2.2 Allocate memory fo the workspace
    lwork = max(1, int(work(1)))
    deallocate(work)
    allocate(work(lwork))

    ! 2.3 compute the matrix Q
    call DORGQR(m, n, min(m, n), basis, m, tau, work, lwork, info)
    call check_lapack_call(info, "DORGQR")
    
    ! release memory
    deallocate(work)
    
  end subroutine lapack_qr
  
  subroutine lapack_solver(arr, brr)
    !> Call lapack DSYSV subroutine to solve a AX=B Linear system
    !> \param arr: matrix with the coefficients of the linear system
    !> \param brr: Vector with the constant terms
    !> \returns: Solution vector X (overwriten brr)

    implicit none
    
    real(dp), dimension(:, :), intent(inout) :: arr, brr
    
    ! local variables
    real(dp), dimension(:), allocatable :: work
    integer :: n, info, lwork
    integer, dimension(size(arr, 1)) :: ipiv
    
    n = size(arr, 1)

    ! query spacework size
    allocate(work(1))
    call DSYSV("U", n, 1, arr, n, ipiv, brr, n, work, -1, info)
    call check_lapack_call(info, "DSYSV")

    ! Allocate memory fo the workspace
    lwork = max(1, int(work(1)))
    deallocate(work)
    allocate(work(lwork))
    
    ! run linear solver
    call DSYSV("U", n, 1, arr, n, ipiv, brr, n, work, lwork, info)
    ! If the diagonalization fails due to a singular value try to recover
    ! by replacing the 0 value with a tiny one
    if (info > 0) then
       arr(info, info) = tiny(arr(1, 1))
       call DSYSV("U", n, 1, arr, n, ipiv, brr, n, work, lwork, info)
       call check_lapack_call(info, "DSYSV")
    end if

    deallocate(work)
    
  end subroutine lapack_solver

    function lapack_matmul(transA, transB, arr, brr, alpha) result (mtx)
    !> perform the matrix multiplication alpha * arr ^ (transA) * brr ^ (transB)
    !> see Lapack DGEMM for further details
    !> \param transA: 'T' transpose A, 'N' do not tranpose
    !> \param transB: 'T' transpose B, 'N' do not tranpose
    !> \param arr: first matrix to multiply
    !> \param brr: second matrix
    !> \param alpha: optional scalar number   
    !> \return matrix multiplication

    implicit none
    
    character(len=1), intent(in) :: transA, transB
    real(dp), dimension(:, :), intent(in) :: arr, brr
    real(dp), optional, intent(in) :: alpha 
    real(dp), dimension(:, :), allocatable :: mtx

    ! local variables
    real(dp) :: x
    integer :: m, n, k, lda, ldb
    x = 1.d0

    ! check optional variable
    if (present(alpha)) x=alpha
    
    if (transA == 'T') then
       k = size(arr, 1)
       m = size(arr, 2)
       lda = k
    else
       k = size(arr, 2)
       m = size(arr, 1)
       lda = m
    end if
    
    if (transB == 'T') then
       n = size(brr, 1)
       ldb = n
    else
       n = size(brr, 2)
       ldb = k
    end if

    ! resulting array
    allocate(mtx(m, n))
    mtx = 0.d0

    call DGEMM(transA, transB, m, n, k, x, arr, lda, brr, ldb, 0.d0, mtx, m)
    

  end function lapack_matmul


  function lapack_matrix_vector(transA, mtx, vector, alpha) result(rs)
    !> perform the Matrix vector multiplication alpha * mtx ^ (transA) * vector
    !> see DGEMV for details
    !> \param transA: 'T' transpose A; 'N' do not transpose
    !> \param mtx: matrix to multiply
    !> \param vector: vector to multiply
    !> \param alpha: optional scalar value
    !> \return resulting vector

    implicit none
    
    character(len=1), intent(in) :: transA
    real(dp), dimension(:, :), intent(in) :: mtx
    real(dp), dimension(:), intent(in) :: vector
    real(dp), optional, intent(in) :: alpha 
    real(dp), dimension(:), allocatable :: rs

    ! local variable
    integer :: m, n
    real(dp) :: scalar
    scalar = 1
    
    ! check optional variable
    if (present(alpha)) scalar=alpha
    
    ! number of row of mtx
    m = size(mtx, 1)
    n = size(mtx, 2)

    allocate(rs(m))
    
    call DGEMV(transA, m, n, scalar, mtx, m, vector, 1, 0.d0, rs, 1)
    

  end function lapack_matrix_vector

  subroutine check_deallocate_matrix(mtx)
    !> deallocate a matrix if allocated
    real(dp), dimension(:, :), allocatable, intent(inout) ::  mtx
    
    if (allocated(mtx)) then
       deallocate(mtx)
    end if

  end subroutine check_deallocate_matrix

  subroutine check_lapack_call(info, name)
    !> Check if a subroutine finishes sucessfully
    !> \param info: Termination signal
    !> \param name: Name of the subroutine
    integer :: info
    character(len=*), intent(in) :: name
    
    if (info /= 0) then
       print *, "call to subroutine: ", name, " has failed!"
       print *, "info: ", info
       error stop
    end if
    
  end subroutine check_lapack_call
  
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
  
  function generate_diagonal_dominant(m, sparsity) result(arr)
    !> Generate a diagonal dominant square matrix of dimension m
    !> \param m: dimension of the matrix
    !> \param sparsity: magnitude order of the off-diagonal values
      
    integer, intent(in) :: m ! size of the square matrix
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
             arr(i, i) = i
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

  
end module davidson


submodule (davidson) correction_methods
  !> submodule containing the implementations of different kind
  !> algorithms to compute the correction vectors for the Davidson's diagonalization

  implicit none
  
contains

  module function compute_correction(mtx, V, eigenvalues, eigenvectors, method) &
       result(correction)
    !> see interface in davidson module
    real(dp), dimension(:), intent(in) :: eigenvalues
    real(dp), dimension(:, :), intent(in) :: mtx, V, eigenvectors
    character(len=*), optional,intent(in) :: method

    ! local variables
    character(len=10) :: opt 
    real(dp), dimension(size(mtx, 1), size(V, 2)) :: correction

    !check optional arguments
    opt="DPR"
    if (present(method)) opt=trim(method)
    
    select case (method)
    case ("DPR")
       correction = compute_DPR(mtx, V, eigenvalues, eigenvectors)
    case ("GJD")
       correction = compute_GJD(mtx, V, eigenvalues, eigenvectors)
    end select
    
  end function compute_correction

  function compute_DPR(mtx, V, eigenvalues, eigenvectors) result(correction)
    !> compute Diagonal-Preconditioned-Residue (DPR) correction
    real(dp), dimension(:), intent(in) :: eigenvalues
    real(dp), dimension(:, :), intent(in) :: mtx, V, eigenvectors
    real(dp), dimension(size(mtx, 1), size(V, 2)) :: correction
    ! local variables
    integer :: j, m
    real(dp), dimension(size(mtx, 1), size(mtx, 2)) :: diag, arr
    real(dp), dimension(size(mtx, 1)) :: brr

    ! shape of matrix
    m = size(mtx, 1)
    
    do j=1, size(V, 2)
       diag = eye(m , m, eigenvalues(j))
       arr = mtx - diag
       brr = matmul(V, eigenvectors(:, j))
       correction(:, j) = matmul(arr, brr) / (eigenvalues(j) - mtx(j, j))       
    end do

  end function compute_DPR

  function compute_GJD(mtx, V, eigenvalues, eigenvectors) result(correction)
    !> Compute the Generalized Jacobi Davidson (GJD) correction
    use ieee_arithmetic, only: ieee_is_normal
    
    real(dp), dimension(:), intent(in) :: eigenvalues
    real(dp), dimension(:, :), intent(in) :: mtx, V, eigenvectors
    real(dp), dimension(size(mtx, 1), size(V, 2)) :: correction

    ! local variables
    integer :: k, m
    real(dp), dimension(size(mtx, 1), 1) :: rs
    real(dp), dimension(size(mtx, 1), size(V, 2)) :: ritz_vectors
    real(dp), dimension(size(mtx, 1), size(mtx, 2)) :: arr, xs, ys
    real(dp), dimension(size(mtx, 1), 1) :: brr

    ! Diagonal matrix
    m = size(mtx, 1)
    ritz_vectors = matmul(V, eigenvectors) !lapack_matmul('N', 'N', V, eigenvectors)
    do k=1, size(V, 2)
       rs(:, 1) = ritz_vectors(:, k)
       xs = eye(m, m) - matmul(rs, transpose(rs))!lapack_matmul('N', 'T', rs, rs)
       ys = substract_from_diagonal(mtx, eigenvalues(k))
       arr = matmul(xs, matmul(ys, xs))
       brr = -rs
       
       call lapack_solver(arr, brr)
       correction(:, k) = brr(:, 1)
    end do


    
  end function compute_GJD

  function substract_from_diagonal(mtx, alpha) result(arr)
    !> susbstract an scalar from the diagonal of a matrix
    !> \param mtx: square matrix
    !> \param alpha: scalar to substract
    real(dp), dimension(:, :), intent(in) :: mtx
    real(dp), dimension(size(mtx, 1), size(mtx, 2)) :: arr
    real(dp), intent(in) :: alpha
    integer :: i

    arr = mtx
    do i=1,size(mtx, 1)
       arr(i, i) = arr(i, i) - alpha
    end do
    
  end function substract_from_diagonal
  
end submodule correction_methods
