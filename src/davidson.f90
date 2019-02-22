!> \namespace Davidson eigensolver
!> \author Felipe Zapata
module davidson

  use numeric_kinds, only: dp
  use lapack_wrapper, only: lapack_generalized_eigensolver, lapack_matmul, lapack_matrix_vector, &
       lapack_qr, lapack_solver
  implicit none

  !> \private
  private
  !> \public
  public :: eye, generalized_eigensolver, generate_diagonal_dominant, norm, free_matmul
  
  interface generalized_eigensolver
  !> \brief Solve a (general) eigenvalue problem using different types of Davidson algorithms.

  !> \param[in] mtx: Matrix to diagonalize
  !> \param[in, opt] stx: Optional matrix for the general eigenvalue problem:
  !> \f$ mtx \lambda = V stx \lambda \f$
     
  !> \param[out] eigenvalues Computed eigenvalues
  !> \param[out] ritz_vectors approximation to the eigenvectors
  !> \param[in] lowest Number of lowest eigenvalues/ritz_vectors to compute
  !> \param[in] method Method to compute the correction vector. Available
  !> methods are,
  !>    DPR: Diagonal-Preconditioned-Residue
  !>    GJD: Generalized Jacobi Davidson
  !> \param[in] max_iters: Maximum number of iterations
  !> \param[in] tolerance norm-2 error of the eigenvalues
  !> \param[in] method: Method to compute the correction vectors
  !> \param[in, opt] max_dim_sub: maximum dimension of the subspace search   
  !> \param[out] iters: Number of iterations until convergence
  !> \return eigenvalues and ritz_vectors of the matrix `mtx`

     module procedure generalized_eigensolver_dense
     module procedure generalized_eigensolver_free
     
  end interface generalized_eigensolver
  
  interface
     
     module subroutine generalized_eigensolver_dense(mtx, eigenvalues, ritz_vectors, lowest, method, max_iters, &
          tolerance, iters, max_dim_sub, stx)
       !> The current implementation uses a general  davidson algorithm, meaning
       !> that it compute all the eigenvalues simultaneusly using a block approach.
       !> The family of Davidson algorithm only differ in the way that the correction
       !> vector is computed.
       
       !> \param[in] mtx: Matrix to diagonalize
       !> \param[in, opt] Optional matrix to solve the general eigenvalue problem:
       !> \f$ mtx \lambda = V stx \lambda \f$
       !> \param[out] eigenvalues Computed eigenvalues
       !> \param[out] ritz_vectors approximation to the eigenvectors
       !> \param[in] lowest Number of lowest eigenvalues/ritz_vectors to compute
       !> \param[in] method Method to compute the correction vector. Available
       !> methods are,
       !>    DPR: Diagonal-Preconditioned-Residue
       !>    GJD: Generalized Jacobi Davidson
       !> \param[in] max_iters: Maximum number of iterations
       !> \param[in] tolerance norm-2 error of the eigenvalues
       !> \param[in] method: Method to compute the correction vectors
       !> \param[in, opt] max_dim_sub: maximum dimension of the subspace search   
       !> \param[out] iters: Number of iterations until convergence
       !> \return eigenvalues and ritz_vectors of the matrix `mtx`
       use numeric_kinds, only: dp
       implicit none
       ! input/output variable
       integer, intent(in) :: lowest
       real(dp), dimension(:, :), intent(in) :: mtx
       real(dp), dimension(:, :), intent(in), optional :: stx
       real(dp), dimension(lowest), intent(out) :: eigenvalues
       real(dp), dimension(:, :), intent(out) :: ritz_vectors
       integer, intent(in) :: max_iters
       integer, intent(in), optional :: max_dim_sub
       real(dp), intent(in) :: tolerance
       character(len=*), intent(in) :: method
       integer, intent(out) :: iters
       
     end subroutine generalized_eigensolver_dense

     module function compute_correction_generalized_dense(mtx, V, eigenvalues, eigenvectors, method, stx) &
          result(correction)
       !> compute the correction vector using a given `method` for the Davidson algorithm
       !> See correction_methods submodule for the implementations
       !> \param[in] mtx: Original matrix
       !> \param[in] stx: Matrix to compute the general eigenvalue problem
       !> \param[in] V: Basis of the iteration subspace
       !> \param[in] eigenvalues: of the reduce problem
       !> \param[in] eigenvectors: of the reduce problem
       !> \param[in] method: name of the method to compute the correction
 
       real(dp), dimension(:), intent(in) :: eigenvalues
       real(dp), dimension(:, :), intent(in) :: mtx, V, eigenvectors
       real(dp), dimension(:, :), intent(in), optional :: stx
       character(len=*), optional, intent(in) :: method
       real(dp), dimension(size(mtx, 1), size(V, 2)) :: correction
       
     end function compute_correction_generalized_dense
     
  end interface
  
contains    

  subroutine generalized_eigensolver_free(fun_mtx_gemv, eigenvalues, ritz_vectors, lowest, method, max_iters, &
       tolerance, iters, max_dim_sub, fun_stx_gemv)
    !> \brief use a pair of functions fun_mtx and fun_stx to compute on the fly the matrices to solve
    !>  the general eigenvalue problem
    !> The current implementation uses a general  davidson algorithm, meaning
    !> that it compute all the eigenvalues simultaneusly using a block approach.
    !> The family of Davidson algorithm only differ in the way that the correction
    !> vector is computed.
    
    !> \param[in] fun_mtx_gemv: Function to apply the matrix to a buncof vectors
    !> \param[in, opt] fun_stx_gemv: (optional) function to apply the pencil to a bunch of vectors.
    !> \param[out] eigenvalues Computed eigenvalues
    !> \param[out] ritz_vectors approximation to the eigenvectors
    !> \param[in] lowest Number of lowest eigenvalues/ritz_vectors to compute
    !> \param[in] method Method to compute the correction vector. Available
    !> methods are,
    !>    DPR: Diagonal-Preconditioned-Residue
    !>    GJD: Generalized Jacobi Davidson
    !> \param[in] max_iters: Maximum number of iterations
    !> \param[in] tolerance norm-2 error of the eigenvalues
    !> \param[in] method: Method to compute the correction vectors
    !> \param[in, opt] max_dim_sub: maximum dimension of the subspace search   
    !> \param[out] iters: Number of iterations until convergence
    !> \return eigenvalues and ritz_vectors of the matrix `mtx`

    implicit none

    ! input/output variable
    integer, intent(in) :: lowest
    real(dp), dimension(lowest), intent(out) :: eigenvalues
    real(dp), dimension(:, :), intent(out) :: ritz_vectors
    integer, intent(in) :: max_iters
    integer, intent(in), optional :: max_dim_sub
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: method
    integer, intent(out) :: iters
    
    ! Function to compute the target matrix on the fly
    interface

       function fun_mtx_gemv(input_vect) result(output_vect)
         !> \brief Function to compute the optional mtx on the fly
         !> \param[in] i column/row to compute from mtx
         !> \param vec column/row from mtx
         use numeric_kinds, only: dp
         real (dp), dimension(:,:), intent(in) :: input_vect
         real (dp), dimension(size(input_vect,1),size(input_vect,2)) :: output_vect

       end function fun_mtx_gemv
       
       function fun_stx_gemv(input_vect) result(output_vect)
         !> \brief Fucntion to compute the optional stx matrix on the fly
         !> \param[in] i column/row to compute from stx
         !> \param vec column/row from stx
         use numeric_kinds, only: dp
         real(dp), dimension(:,:), intent(in) :: input_vect
         real (dp), dimension(size(input_vect,1),size(input_vect,2)) :: output_vect
         
       end function fun_stx_gemv

    end interface
    
    !local variables
    integer :: dim_mtx, dim_sub, max_dim, i, j
    
    ! ! Basis of subspace of approximants
    real(dp), dimension(size(ritz_vectors, 1),1) :: guess, rs
    real(dp), dimension(size(ritz_vectors, 1)  ) :: diag_mtx, diag_stx
    real(dp), dimension(lowest):: errors
    
    ! ! Working arrays
    real(dp), dimension(:), allocatable :: eigenvalues_sub
    real(dp), dimension(:, :), allocatable :: correction, eigenvectors_sub, mtx_proj, stx_proj, V
    
    ! Iteration subpsace dimension
    dim_sub = lowest * 2
    
    ! maximum dimension of the basis for the subspace
    if (present(max_dim_sub)) then
       max_dim  = max_dim_sub
    else
       max_dim = lowest * 10
    endif
    
    ! dimension of the matrix
    dim_mtx = size(ritz_vectors, 1)

    ! extract the diagonals of the matrices
    diag_mtx = extract_diagonal_free(fun_mtx_gemv,dim_mtx)
    diag_stx = extract_diagonal_free(fun_stx_gemv,dim_mtx)
    
    ! 1. Variables initialization
    V = eye(dim_mtx, dim_sub) ! Initial orthonormal basis
    
    ! 2. Generate subspace matrix problem by projecting into V
    mtx_proj = lapack_matmul('T', 'N', V, fun_mtx_gemv(V))
    stx_proj = lapack_matmul('T', 'N', V, fun_stx_gemv(V))

    
    ! Outer loop block Davidson schema
    outer_loop: do i=1, max_iters

       ! 3. compute the eigenvalues and their corresponding ritz_vectors
       ! for the projected matrix using lapack
       call check_deallocate_matrix(eigenvectors_sub)
       
       if (allocated(eigenvalues_sub)) then
          deallocate(eigenvalues_sub)
       end if
       
       allocate(eigenvalues_sub(size(mtx_proj, 1)))
       allocate(eigenvectors_sub(size(mtx_proj, 1), size(mtx_proj, 2)))
       
       call lapack_generalized_eigensolver(mtx_proj, eigenvalues_sub, eigenvectors_sub, stx_proj)
       
       ! 4. Check for convergence
       ritz_vectors = lapack_matmul('N', 'N', V, eigenvectors_sub(:, :lowest))
       do j=1,lowest
          guess = eigenvalues_sub(j) * fun_stx_gemv(reshape(ritz_vectors(:, j),(/dim_mtx,1/) ) )
          rs = fun_mtx_gemv(reshape(ritz_vectors(:, j), (/dim_mtx,1/))) - guess
          errors(j) = norm(reshape(rs,(/dim_mtx/)))
       end do
       
       if (all(errors < tolerance)) then
          iters = i
          exit
       end if
       
       ! 5. Add the correction vectors to the current basis
       if (size(V, 2) <= max_dim) then
          
          ! append correction to the current basis
          call check_deallocate_matrix(correction)
          allocate(correction(size(ritz_vectors, 1), size(V, 2)))
          
          correction = compute_DPR_free(fun_mtx_gemv, fun_stx_gemv, V, eigenvalues_sub, eigenvectors_sub, diag_mtx, diag_stx)
          
          ! 6. Increase Basis size
          call concatenate(V, correction)
          
          ! 7. Orthogonalize basis
          call lapack_qr(V)

          ! 8. Update the the projection 
          call update_projection_free(fun_mtx_gemv, V, mtx_proj)
          call update_projection_free(fun_stx_gemv, V, stx_proj)
          
       else
          
          ! 6. Otherwise reduce the basis of the subspace to the current correction
          V = lapack_matmul('N', 'N', V, eigenvectors_sub(:, :dim_sub))
          
          ! we refresh the projected matrices
          mtx_proj = lapack_matmul('T', 'N', V, fun_mtx_gemv(V))
          stx_proj = lapack_matmul('T', 'N', V, fun_stx_gemv(V))
          
       end if
       
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
    deallocate(eigenvalues_sub, eigenvectors_sub, V, mtx_proj)
    
    ! free optional matrix
    call check_deallocate_matrix(stx_proj)
    
  end subroutine generalized_eigensolver_free
    
    function compute_DPR_free(fun_mtx_gemv, fun_stx_gemv, V, eigenvalues, eigenvectors, diag_mtx, diag_stx) result(correction)

      !> compute the correction vector using the DPR method for a matrix free diagonalization
      !> See correction_methods submodule for the implementations
      !> \param[in] fun_mtx: function to compute matrix
      !> \param[in] fun_stx: function to compute the matrix for the generalized case
      !> \param[in] V: Basis of the iteration subspace
      !> \param[in] eigenvalues: of the reduce problem
      !> \param[in] eigenvectors: of the reduce problem
      !> \return correction matrix
      
      real(dp), dimension(:), intent(in) :: eigenvalues
      real(dp), dimension(:, :), intent(in) :: V, eigenvectors
      real(dp), dimension(:), intent(in) :: diag_mtx, diag_stx


      ! Function to compute the target matrix on the fly
      interface

         function fun_mtx_gemv(input_vect) result(output_vect)
           !> \brief Function to compute the optional mtx on the fly
           !> \param[in] i column/row to compute from mtx
           !> \param vec column/row from mtx
           use numeric_kinds, only: dp
           real (dp), dimension(:,:), intent(in) :: input_vect
           real (dp), dimension(size(input_vect,1),size(input_vect,2)) :: output_vect

         end function fun_mtx_gemv
         
         function fun_stx_gemv(input_vect) result(output_vect)
           !> \brief Fucntion to compute the optional stx matrix on the fly
           !> \param[in] i column/row to compute from stx
           !> \param vec column/row from stx
           use numeric_kinds, only: dp
           real(dp), dimension(:,:), intent(in) :: input_vect
           real (dp), dimension(size(input_vect,1),size(input_vect,2)) :: output_vect
           
         end function fun_stx_gemv

      end interface
      
      ! local variables
      !real(dp), dimension(size(V, 1),1) :: vector
      real(dp), dimension(size(V, 1), size(V, 2)) :: correction, vectors
      real(dp), dimension(size(V, 2),size(V, 2)) :: diag_eigenvalues
      integer :: ii, j
      integer :: dim
            
      dim = size(V,1)

      ! create a diagonal matrix of eigenvalues
      diag_eigenvalues = 0E0
      do ii =1, size(V,2)
        diag_eigenvalues(ii,ii) = eigenvalues(ii)
      end do

      ! create all the ritz vectors
      vectors = lapack_matmul('N','N', V, eigenvectors)

      ! initialize the correction vectors
      correction = fun_mtx_gemv(vectors) - lapack_matmul('N','N',fun_stx_gemv(vectors), diag_eigenvalues)

      do j=1, size(V, 2)
         do ii=1,size(correction,1)
            correction(ii, j) = correction(ii, j) / (eigenvalues(j) * diag_stx(ii)  - diag_mtx(ii))
         end do
      end do
      
    end function compute_DPR_free
  

  function extract_diagonal_free(fun_A_gemv,dim) result(out)
    !> \brief extract the diagonal of the matrix
    !> \param dim: dimension of the matrix


    implicit none
    integer, intent(in) :: dim
    real(dp), dimension(dim) :: out


    interface
         function fun_A_gemv(input_vect) result(output_vect)
           !> \brief Function to compute the optional mtx on the fly
           !> \param[in] i column/row to compute from mtx
           !> \param vec column/row from mtx
           use numeric_kinds, only: dp
           real (dp), dimension(:,:), intent(in) :: input_vect
           real (dp), dimension(size(input_vect,1),size(input_vect,2)) :: output_vect

         end function fun_A_gemv
    end interface

    ! local variable
    integer :: ii
    real(dp), dimension(dim,1) :: tmp_array
    
    do ii = 1,dim
      tmp_array = 0E0
      tmp_array(ii,1) = 1.0
      tmp_array = fun_A_gemv(tmp_array)
      out(ii) = tmp_array(ii,1)
    end do

  end function extract_diagonal_free


  subroutine update_projection_free(fun_A_gemv, V, A_proj)
    !> \brief update the projected matrices
    !> \param A: full matrix
    !> \param V: projector
    !> \param A_proj: projected matrix

    implicit none
    real(dp), dimension(:, :), intent(in) :: V
    real(dp), dimension(:, :), intent(inout), allocatable :: A_proj
    real(dp), dimension(:, :), allocatable :: tmp_array

    interface
         function fun_A_gemv(input_vect) result(output_vect)
           !> \brief Function to compute the optional mtx on the fly
           !> \param[in] i column/row to compute from mtx
           !> \param vec column/row from mtx
           use numeric_kinds, only: dp
           real (dp), dimension(:,:), intent(in) :: input_vect
           real (dp), dimension(size(input_vect,1),size(input_vect,2)) :: output_vect

         end function fun_A_gemv
    end interface
    
    ! local variables
    integer :: nvec, old_dim

    ! dimension of the matrices
    nvec = size(V,2)
    old_dim = size(A_proj,1)    

    ! move to temporal array
    allocate(tmp_array(nvec, nvec))
    tmp_array(:old_dim, :old_dim) = A_proj
    tmp_array(:,old_dim+1:) = lapack_matmul('T', 'N', V, fun_A_gemv(V(:, old_dim+1:)))
    tmp_array( old_dim+1:,:old_dim ) = transpose(tmp_array(:old_dim, old_dim+1:))
    
    ! Move to new expanded matrix
    deallocate(A_proj)
    call move_alloc(tmp_array, A_proj)

  end subroutine update_projection_free

  function free_matmul(fun, array) result (mtx)
    !> \brief perform a matrix-matrix multiplication by generating a matrix on the fly using `fun`
    !> \param[in] fun function to compute a matrix on the fly
    !> \param[in] array matrix to multiply with fun
    !> \return resulting matrix

    ! input/output
    implicit none
    real(dp), dimension(:, :), intent(in) :: array
    real(dp), dimension(size(array, 1), size(array, 2)) :: mtx

    interface
       function fun(i, dim) result(vec)
         !> \brief Fucntion to compute the matrix `mtx` on the fly
         !> \param[in] i column/row to compute from `mtx`
         !> \param vec column/row from mtx
         use numeric_kinds, only: dp
         integer, intent(in) :: i
         integer, intent(in) :: dim         
         real(dp), dimension(dim) :: vec

       end function fun

    end interface

    ! local variables
    real(dp), dimension(size(array, 1)) :: vec
    integer :: dim1, dim2, i, j

    ! dimension of the square matrix computed on the fly
    dim1 = size(array, 1)
    dim2 = size(array, 2)

    !$OMP PARALLEL DO &
    !$OMP PRIVATE(i, j, vec)
    do i = 1, dim1
       vec = fun(i, dim1)
       do j = 1, dim2
          mtx(i, j) = dot_product(vec, array(:, j))
       end do
    end do
    !$OMP END PARALLEL DO
    
  end function free_matmul

  subroutine check_deallocate_matrix(mtx)
    !> deallocate a matrix if allocated
    real(dp), dimension(:, :), allocatable, intent(inout) ::  mtx
    
    if (allocated(mtx)) then
       deallocate(mtx)
    end if

  end subroutine check_deallocate_matrix
  
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
    
  function generate_diagonal_dominant(m, sparsity, diag_val) result(arr)
    !> Generate a diagonal dominant square matrix of dimension m
    !> \param m dimension of the matrix
    !> \param sparsity magnitude order of the off-diagonal values
      
    integer, intent(in) :: m ! size of the square matrix
    real(dp), optional :: diag_val
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
            if (present(diag_val))then
              arr(i,i) = diag_val
            else
             arr(i, i) = i
            end if
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


submodule (davidson) davidson_dense
  !> Submodule containing the implementation of the Davidson diagonalization method
  !> for dense matrices
  use lapack_wrapper, only: lapack_generalized_eigensolver, lapack_matmul, lapack_matrix_vector, &
       lapack_qr, lapack_solver
  
  implicit none

  contains
  
    module subroutine generalized_eigensolver_dense(mtx, eigenvalues, ritz_vectors, lowest, method, max_iters, &
        tolerance, iters, max_dim_sub, stx)
    !> The current implementation uses a general  davidson algorithm, meaning
    !> that it compute all the eigenvalues simultaneusly using a block approach.
    !> The family of Davidson algorithm only differ in the way that the correction
    !> vector is computed.
    
    !> \param[in] mtx: Matrix to diagonalize
    !> \param[in, opt] Optional matrix to solve the general eigenvalue problem:
    !> \f$ mtx \lambda = V stx \lambda \f$
    !> \param[out] eigenvalues Computed eigenvalues
    !> \param[out] ritz_vectors approximation to the eigenvectors
    !> \param[in] lowest Number of lowest eigenvalues/ritz_vectors to compute
    !> \param[in] method Method to compute the correction vector. Available
    !> methods are,
    !>    DPR: Diagonal-Preconditioned-Residue
    !>    GJD: Generalized Jacobi Davidson
    !> \param[in] max_iters: Maximum number of iterations
    !> \param[in] tolerance norm-2 error of the eigenvalues
    !> \param[in] method: Method to compute the correction vectors
    !> \param[in, opt] max_dim_sub: maximum dimension of the subspace search   
    !> \param[out] iters: Number of iterations until convergence
    !> \return eigenvalues and ritz_vectors of the matrix `mtx`

    implicit none
    ! input/output variable
    integer, intent(in) :: lowest
    real(dp), dimension(:, :), intent(in) :: mtx
    real(dp), dimension(:, :), intent(in), optional :: stx
    real(dp), dimension(lowest), intent(out) :: eigenvalues
    real(dp), dimension(:, :), intent(out) :: ritz_vectors
    integer, intent(in) :: max_iters
    integer, intent(in), optional :: max_dim_sub
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
    real(dp), dimension(:, :), allocatable :: correction, eigenvectors_sub, mtx_proj, stx_proj, V

    ! generalize problem
    logical :: gev 

    ! Iteration subpsace dimension
    dim_sub = lowest * 2

    ! maximum dimension of the basis for the subspace
    if (present(max_dim_sub)) then
       max_dim  = max_dim_sub
    else
       max_dim = lowest * 10
    endif

    ! generalied problem
    gev = present(stx)

    ! 1. Variables initialization
    V = eye(size(ritz_vectors, 1), dim_sub) ! Initial orthonormal basis


   ! 2. Generate subpace matrix problem by projecting into V
   mtx_proj = lapack_matmul('T', 'N', V, lapack_matmul('N', 'N', mtx, V))

   if(gev) then
    stx_proj = lapack_matmul('T', 'N', V, lapack_matmul('N', 'N', stx, V))
   end if

    ! ! Outer loop block Davidson schema
    outer_loop: do i=1, max_iters

       ! 3. compute the eigenvalues and their corresponding ritz_vectors
       ! for the projected matrix using lapack
       call check_deallocate_matrix(eigenvectors_sub)

       if (allocated(eigenvalues_sub)) then
          deallocate(eigenvalues_sub)
       end if

       allocate(eigenvalues_sub(size(mtx_proj, 1)))
       allocate(eigenvectors_sub(size(mtx_proj, 1), size(mtx_proj, 2)))


       if (gev) then
        call lapack_generalized_eigensolver(mtx_proj, eigenvalues_sub, eigenvectors_sub, stx_proj)
       else
        call lapack_generalized_eigensolver(mtx_proj, eigenvalues_sub, eigenvectors_sub)
       end if

       ! 4. Check for convergence
       ritz_vectors = lapack_matmul('N', 'N', V, eigenvectors_sub(:, :lowest))
       do j=1,lowest
          if(gev) then
            guess = eigenvalues_sub(j) * lapack_matrix_vector('N',stx,ritz_vectors(:, j))
          else
            guess = eigenvalues_sub(j) * ritz_vectors(:, j)
          end if
          rs = lapack_matrix_vector('N', mtx, ritz_vectors(:, j)) - guess
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

          if(gev) then
            correction = compute_correction_generalized_dense(mtx, V, eigenvalues_sub, eigenvectors_sub, method, stx)
          else
            correction = compute_correction_generalized_dense(mtx, V, eigenvalues_sub, eigenvectors_sub, method)
          end if


          ! 6. Increase Basis size
          call concatenate(V, correction)
       
          ! 7. Orthogonalize basis
          call lapack_qr(V)
          

          ! 8. Update the the projection 
          call update_projection_dense(mtx, V, mtx_proj)
          if (gev) then
             call update_projection_dense(stx, V, stx_proj)
           end if

       else

          ! 6. Otherwise reduce the basis of the subspace to the current correction
          V = lapack_matmul('N', 'N', V, eigenvectors_sub(:, :dim_sub))

          ! we refresh the projected matrices
          mtx_proj = lapack_matmul('T', 'N', V, lapack_matmul('N', 'N', mtx, V))

          if(gev) then
           stx_proj = lapack_matmul('T', 'N', V, lapack_matmul('N', 'N', stx, V))
          end if

       end if

      
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
    deallocate(eigenvalues_sub, eigenvectors_sub, V, mtx_proj)

    ! free optional matrix
    if (gev) then
       call check_deallocate_matrix(stx_proj)
    endif
    
  end subroutine generalized_eigensolver_dense

  subroutine update_projection_dense(A, V, A_proj)
    !> \brief update the projected matrices
    !> \param A: full matrix
    !> \param V: projector
    !> \param A_proj: projected matrix

    implicit none
    real(dp), dimension(:, :), intent(in) :: A
    real(dp), dimension(:, :), intent(in) :: V
    real(dp), dimension(:, :), intent(inout), allocatable :: A_proj
    real(dp), dimension(:, :), allocatable :: tmp_array

    ! local variables
    integer :: nvec, old_dim

    ! dimension of the matrices
    nvec = size(V,2)
    old_dim = size(A_proj,1)    

    ! move to temporal array
    allocate(tmp_array(nvec, nvec))
    tmp_array(:old_dim, :old_dim) = A_proj
    tmp_array(:,old_dim+1:) = lapack_matmul('T', 'N', V, lapack_matmul('N', 'N', A, V(:, old_dim+1:)))
    tmp_array( old_dim+1:,:old_dim ) = transpose(tmp_array(:old_dim, old_dim+1:))
    
    ! Move to new expanded matrix
    deallocate(A_proj)
    call move_alloc(tmp_array, A_proj)
 
  end subroutine update_projection_dense
  
end submodule davidson_dense

submodule (davidson) correction_methods_generalized_dense
  !> submodule containing the implementations of different kind
  !> algorithms to compute the correction vectors for the Davidson's diagonalization

  implicit none
  
contains

  module function compute_correction_generalized_dense(mtx, V, eigenvalues, eigenvectors, method, stx) &
       result(correction)
    !> see interface in davidson module
    real(dp), dimension(:), intent(in) :: eigenvalues
    real(dp), dimension(:, :), intent(in) :: mtx, V, eigenvectors
    real(dp), dimension(:, :), intent(in), optional :: stx
    character(len=*), optional,intent(in) :: method
    logical :: gev 

    ! local variables
    character(len=10) :: opt 
    real(dp), dimension(size(mtx, 1), size(V, 2)) :: correction

    !check optional arguments
    gev = present(stx)
    opt="DPR"
    if (present(method)) opt=trim(method)
    
    select case (method)
    case ("DPR")
      if(gev) then
       correction = compute_DPR_generalized_dense(mtx, V, eigenvalues, eigenvectors, stx)
      else
        correction = compute_DPR_generalized_dense(mtx, V, eigenvalues, eigenvectors)
      end if
    case ("GJD")
      if(gev) then
       correction = compute_GJD_generalized_dense(mtx, V, eigenvalues, eigenvectors, stx)
      else
        correction = compute_GJD_generalized_dense(mtx, V, eigenvalues, eigenvectors)
      end if
    end select
    
  end function compute_correction_generalized_dense

  function compute_DPR_generalized_dense(mtx, V, eigenvalues, eigenvectors, stx) result(correction)
    !> compute Diagonal-Preconditioned-Residue (DPR) correction
    real(dp), dimension(:), intent(in) :: eigenvalues
    real(dp), dimension(:, :), intent(in) :: mtx, V, eigenvectors
    real(dp), dimension(:, :), intent(in), optional ::  stx 
    real(dp), dimension(size(mtx, 1), size(V, 2)) :: correction

    ! local variables
    integer :: ii,j, m
    real(dp), dimension(size(mtx, 1), size(mtx, 2)) :: diag, arr
    real(dp), dimension(size(mtx, 1)) :: vec
    logical :: gev

    ! shape of matrix
    m = size(mtx, 1)
    gev = (present(stx))

    do j=1, size(V, 2)
       if(gev) then
          diag = eigenvalues(j) * stx
       else
          diag = eye(m , m, eigenvalues(j))
       end if
       arr = mtx - diag
       vec = lapack_matrix_vector('N', V, eigenvectors(:, j))
      
       correction(:, j) = lapack_matrix_vector('N', arr, vec) 

       do ii=1,size(correction,1)
        if (gev) then
          correction(ii, j) = correction(ii, j) / (eigenvalues(j) * stx(ii,ii) - mtx(ii, ii))
        else 
          correction(ii, j) = correction(ii, j) / (eigenvalues(j)  - mtx(ii, ii))
        end if
       end do
    end do

  end function compute_DPR_generalized_dense

  function compute_GJD_generalized_dense(mtx, V, eigenvalues, eigenvectors, stx) result(correction)
    !> Compute the Generalized Jacobi Davidson (GJD) correction
    
    real(dp), dimension(:), intent(in) :: eigenvalues
    real(dp), dimension(:, :), intent(in) :: mtx, V, eigenvectors
    real(dp), dimension(:, :), intent(in), optional :: stx
    real(dp), dimension(size(mtx, 1), size(V, 2)) :: correction

    ! local variables
    integer :: k, m
    logical :: gev
    real(dp), dimension(size(mtx, 1), 1) :: rs
    real(dp), dimension(size(mtx, 1), size(V, 2)) :: ritz_vectors
    real(dp), dimension(size(mtx, 1), size(mtx, 2)) :: arr, xs, ys
    real(dp), dimension(size(mtx, 1), 1) :: brr

    ! Diagonal matrix
    m = size(mtx, 1)
    ritz_vectors = lapack_matmul('N', 'N', V, eigenvectors)

    gev = present(stx)

    do k=1, size(V, 2)
       rs(:, 1) = ritz_vectors(:, k)
       xs = eye(m, m) - lapack_matmul('N', 'T', rs, rs)
       if(gev) then
        ys = mtx - eigenvalues(k)*stx
       else
         ys = substract_from_diagonal(mtx, eigenvalues(k))
       end if
       arr = lapack_matmul('N', 'N', xs, lapack_matmul('N', 'N', ys, xs))
       brr = -rs
       
       call lapack_solver(arr, brr)
       correction(:, k) = brr(:, 1)
    end do
    
  end function compute_GJD_generalized_dense

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
  
end submodule correction_methods_generalized_dense
