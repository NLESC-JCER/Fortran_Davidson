!> \namespace Davidson eigensolver
!> \author Felipe Zapata
!> The current implementation uses a general  davidson algorithm, meaning
!> that it compute all the eigenvalues simultaneusly using a variable size block approach.
!> The family of Davidson algorithm only differ in the way that the correction
!> vector is computed.
!> Computed pairs of eigenvalues/eigenvectors are deflated using algorithm
!> described at: https://doi.org/10.1023/A:101919970


module davidson_dense
  !> Submodule containing the implementation of the Davidson diagonalization method
  !> for dense matrices
  use numeric_kinds, only: dp
  use lapack_wrapper, only: lapack_generalized_eigensolver, lapack_matmul, lapack_matrix_vector, &
       lapack_qr, lapack_solver, lapack_sort
  use array_utils, only: concatenate, diagonal, eye, generate_preconditioner, norm, modified_gram_schmidt

  implicit none
  
  !> \private
  private
  !> \public
  public :: generalized_eigensolver_dense

  interface
     module function compute_correction_generalized_dense(matrix, eigenvalues, ritz_vectors, residues, method, second_matrix) &
          result(correction)
       !> compute the correction vector using a given `method` for the Davidson algorithm
       !> See correction_methods submodule for the implementations
       !> \param[in] matrix: Original matrix
       !> \param[in] second_matrix: Matrix to compute the general eigenvalue problem
       !> \param[in] V: Basis of the iteration subspace
       !> \param[in] eigenvalues: of the reduce problem
       !> \param[in] ritz_vectors: guess eigenvectors
       !> \param[in] residues: matrix of residues
       !> \param[in] method: name of the method to compute the correction
       
       real(dp), dimension(:), intent(in) :: eigenvalues
       real(dp), dimension(:, :), intent(in) :: matrix, ritz_vectors, residues
       real(dp), dimension(:, :), intent(in), optional :: second_matrix
       character(len=*), optional, intent(in) :: method
       real(dp), dimension(size(matrix, 1), size(ritz_vectors, 2)) :: correction
       
     end function compute_correction_generalized_dense

  end interface
  
contains

    subroutine generalized_eigensolver_dense(matrix, eigenvalues, eigenvectors, lowest, method, max_iterations, &
        tolerance, iters, max_dim_sub, second_matrix)
     !> Implementation storing in memory the initial densed matrix.
      
    !> \param[in] matrix: Matrix to diagonalize
    !> \param[in, opt] Optional matrix to solve the general eigenvalue problem:
    !> \f$ matrix \lambda = V second_matrix \lambda \f$
    !> \param[out] eigenvalues Computed eigenvalues
    !> \param[out] eigenvectors approximation to the eigenvectors
    !> \param[in] lowest Number of lowest eigenvalues/eigenvectors to compute
    !> \param[in] method Method to compute the correction vector. Available
    !> methods are,
    !>    DPR: Diagonal-Preconditioned-Residue
    !>    GJD: Generalized Jacobi Davidson
    !> \param[in] max_iterations: Maximum number of iterations
    !> \param[in] tolerance norm-2 error of the eigenvalues
    !> \param[in] method: Method to compute the correction vectors
    !> \param[in, opt] max_dim_sub: maximum dimension of the subspace search   
    !> \param[out] iters: Number of iterations until convergence
    !> \return eigenvalues and eigenvectors of the matrix `matrix`

    implicit none
    ! input/output variable
    integer, intent(in) :: lowest
    real(dp), dimension(:, :), intent(in) :: matrix
    real(dp), dimension(:, :), intent(in), optional :: second_matrix
    real(dp), dimension(lowest), intent(out) :: eigenvalues
    real(dp), dimension(:, :), intent(out) :: eigenvectors
    integer, intent(in) :: max_iterations
    integer, intent(in), optional :: max_dim_sub
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: method
    integer, intent(out) :: iters
    
    !local variables
    integer :: i, j, initial_dimension, max_dim
    integer :: n_converged ! Number of converged eigenvalue/eigenvector pairs
    
    ! Basis of subspace of approximants
    real(dp), dimension(size(matrix, 1)) :: guess
    real(dp), dimension(lowest):: errors

    ! Working arrays
    real(dp), dimension(:), allocatable :: eigenvalues_sub
    real(dp), dimension(:, :), allocatable :: correction, eigenvectors_sub, matrix_proj, second_matrix_proj, V
    real(dp), dimension(:, :), allocatable ::residues, ritz_vectors

    ! Diagonal matrix
    real(dp), dimension(size(matrix, 1)) :: d
    
    ! generalize problem
    logical :: gev 

    ! indices of the eigenvalues/eigenvectors pair that have not converged
    logical, dimension(lowest) :: has_converged

    ! Iteration subpsace dimension
    initial_dimension = lowest * 2

    ! Initial number of converged eigenvalue/eigenvector pairs
    n_converged = 0
    has_converged = .False.
    
    ! maximum dimension of the basis for the subspace
    if (present(max_dim_sub)) then
       max_dim  = max_dim_sub
    else
       max_dim = lowest * 10
    endif

    ! generalied problem
    gev = present(second_matrix)

    ! 1. Variables initialization
    ! Select the initial ortogonal subspace based on lowest elements
    ! of the diagonal of the matrix
    d = diagonal(matrix)
    V = generate_preconditioner(d, initial_dimension)

   ! 2. Generate subpace matrix problem by projecting into V
   matrix_proj = lapack_matmul('T', 'N', V, lapack_matmul('N', 'N', matrix, V))

   if(gev) then
    second_matrix_proj = lapack_matmul('T', 'N', V, lapack_matmul('N', 'N', second_matrix, V))
   end if

    ! ! Outer loop block Davidson schema
    outer_loop: do i=1, max_iterations

       ! 3. compute the eigenvalues and their corresponding ritz_vectors
       ! for the projected matrix using lapack
       call check_deallocate_matrix(eigenvectors_sub)

       if (allocated(eigenvalues_sub)) then
          deallocate(eigenvalues_sub)
       end if

       allocate(eigenvalues_sub(size(matrix_proj, 1)))
       allocate(eigenvectors_sub(size(matrix_proj, 1), size(matrix_proj, 2)))


       if (gev) then
        call lapack_generalized_eigensolver(matrix_proj, eigenvalues_sub, eigenvectors_sub, second_matrix_proj)
       else
        call lapack_generalized_eigensolver(matrix_proj, eigenvalues_sub, eigenvectors_sub)
       end if

       ! 4.1 Compute residues
       ritz_vectors = lapack_matmul('N', 'N', V, eigenvectors_sub)

       call check_deallocate_matrix(residues)
       allocate(residues(size(matrix, 1), size(V, 2)))
       do j=1, size(V, 2)
          if(gev) then
             guess = eigenvalues_sub(j) * lapack_matrix_vector('N',second_matrix,ritz_vectors(:, j))
          else
             guess = eigenvalues_sub(j) * ritz_vectors(:, j)
          end if
          residues(:, j) = lapack_matrix_vector('N', matrix, ritz_vectors(:, j)) - guess
       end do

       ! 4.2 Check for convergence
       do j = 1,lowest
          errors(j) = norm(residues(:, j))
          if (errors(j) < tolerance) then
             has_converged(j) = .true.
          end if
       end do

       
       ! Count converged pairs of eigenvalues/eigenvectors
       n_converged = n_converged + count(errors < tolerance)

       ! Select the lowest eigenvalues and their corresponding ritz_vectors
       ! They are sort in increasing order
       eigenvalues = eigenvalues_sub(:lowest)
       eigenvectors = ritz_vectors(:,:lowest)
       
       if (all(has_converged)) then
          iters = i
          exit
       end if
       
       ! 5. Add the correction vectors to the current basis
       if (size(V, 2) <= max_dim) then

          ! append correction to the current basis
          call check_deallocate_matrix(correction)
          allocate(correction(size(matrix, 1), size(V, 2)))

          if(gev) then
             correction = compute_correction_generalized_dense(matrix, eigenvalues_sub, ritz_vectors, &
             residues, method, second_matrix)
          else
            correction = compute_correction_generalized_dense(matrix, eigenvalues_sub, ritz_vectors, residues, method)
          end if


          ! 6. Increase Basis size
          call concatenate(V, correction)
       
          ! 7. Orthogonalize basis
          ! call lapack_qr(V)
          call modified_gram_schmidt(V)

          ! update the projected matrices
          call update_projection_dense(matrix, V, matrix_proj)
          if(gev) then
            call update_projection_dense(second_matrix, V, second_matrix_proj)
          end if

       else

          ! 6. Otherwise reduce the basis of the subspace to the current correction
          V = lapack_matmul('N', 'N', V, eigenvectors_sub(:, :initial_dimension))

          ! we refresh the projected matrices
          matrix_proj = lapack_matmul('T', 'N', V, lapack_matmul('N', 'N', matrix, V))
          if(gev) then
             second_matrix_proj = lapack_matmul('T', 'N', V, lapack_matmul('N', 'N', second_matrix, V))
          end if

       end if

    end do outer_loop

    !  8. Check convergence
    if (i > max_iterations) then
       iters = i
       print *, "Warning: Algorithm did not converge!!"
    end if
        
    ! Free memory
    call check_deallocate_matrix(correction)
    deallocate(eigenvalues_sub, eigenvectors_sub, V, matrix_proj)

    ! free optional matrix
    if (gev) then
       call check_deallocate_matrix(second_matrix_proj)
    endif
    
  end subroutine generalized_eigensolver_dense

  subroutine update_projection_dense(A, V, A_proj)
   !> update the projected matrices
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

  subroutine check_deallocate_matrix(matrix)
    !> deallocate a matrix if allocated
    real(dp), dimension(:, :), allocatable, intent(inout) ::  matrix
    
    if (allocated(matrix)) then
       deallocate(matrix)
    end if
    
  end subroutine check_deallocate_matrix  
  
end module davidson_dense


module davidson_free

  use numeric_kinds, only: dp
  use lapack_wrapper, only: lapack_generalized_eigensolver, lapack_matmul, lapack_matrix_vector, &
       lapack_qr, lapack_solver
  use array_utils, only: concatenate, eye, generate_preconditioner, norm
  use davidson_dense, only: generalized_eigensolver_dense
  implicit none

  !> \private
  private
  !> \public
  public :: generalized_eigensolver_free, free_matmul
  
contains

 subroutine generalized_eigensolver_free(fun_matrix_gemv, eigenvalues, ritz_vectors, lowest, method, max_iterations, &
       tolerance, iters, max_dim_sub, fun_second_matrix_gemv)
    !> \brief use a pair of functions fun_matrix and fun_second_matrix to compute on the fly the matrices to solve
    !>  the general eigenvalue problem
    !> The current implementation uses a general  davidson algorithm, meaning
    !> that it compute all the eigenvalues simultaneusly using a block approach.
    !> The family of Davidson algorithm only differ in the way that the correction
    !> vector is computed.
    
    !> \param[in] fun_matrix_gemv: Function to apply the matrix to a buncof vectors
    !> \param[in, opt] fun_second_matrix_gemv: (optional) function to apply the pencil to a bunch of vectors.
    !> \param[out] eigenvalues Computed eigenvalues
    !> \param[out] ritz_vectors approximation to the eigenvectors
    !> \param[in] lowest Number of lowest eigenvalues/ritz_vectors to compute
    !> \param[in] method Method to compute the correction vector. Available
    !> methods are,
    !>    DPR: Diagonal-Preconditioned-Residue
    !>    GJD: Generalized Jacobi Davidson
    !> \param[in] max_iterations: Maximum number of iterations
    !> \param[in] tolerance norm-2 error of the eigenvalues
    !> \param[in] method: Method to compute the correction vectors
    !> \param[in, opt] max_dim_sub: maximum dimension of the subspace search   
    !> \param[out] iters: Number of iterations until convergence
    !> \return eigenvalues and ritz_vectors of the matrix `matrix`

    implicit none

    ! input/output variable
    integer, intent(in) :: lowest
    real(dp), dimension(lowest), intent(out) :: eigenvalues
    real(dp), dimension(:, :), intent(out) :: ritz_vectors
    integer, intent(in) :: max_iterations
    integer, intent(in), optional :: max_dim_sub
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: method
    integer, intent(out) :: iters
    
    ! Function to compute the target matrix on the fly
    interface

       function fun_matrix_gemv(input_vect) result(output_vect)
         !> \brief Function to compute the optional matrix on the fly
         !> \param[in] i column/row to compute from matrix
         !> \param vec column/row from matrix
         use numeric_kinds, only: dp
         real (dp), dimension(:,:), intent(in) :: input_vect
         real (dp), dimension(size(input_vect,1),size(input_vect,2)) :: output_vect

       end function fun_matrix_gemv
       
       function fun_second_matrix_gemv(input_vect) result(output_vect)
         !> \brief Fucntion to compute the optional second_matrix matrix on the fly
         !> \param[in] i column/row to compute from second_matrix
         !> \param vec column/row from second_matrix
         use numeric_kinds, only: dp
         real(dp), dimension(:,:), intent(in) :: input_vect
         real (dp), dimension(size(input_vect,1),size(input_vect,2)) :: output_vect
         
       end function fun_second_matrix_gemv

    end interface
    
    !local variables
    integer :: dim_matrix, initial_dimension, max_dim, i, j
    
    ! ! Basis of subspace of approximants
    real(dp), dimension(size(ritz_vectors, 1)  ) :: diag_matrix, diag_second_matrix, copy_d
    real(dp), dimension(lowest):: errors
    
    ! ! Working arrays
    real(dp), dimension(:), allocatable :: eigenvalues_sub
    real(dp), dimension(:, :), allocatable :: correction, eigenvectors_sub, guess, lambda, matrix_proj
    real(dp), dimension(:, :), allocatable :: V, second_matrix_proj, matrixV, second_matrixV, residues
    
    ! Iteration subpsace dimension
    initial_dimension = lowest * 2
    
    ! maximum dimension of the basis for the subspace
    if (present(max_dim_sub)) then
       max_dim  = max_dim_sub
    else
       max_dim = lowest * 10
    endif
    
    ! dimension of the matrix
    dim_matrix = size(ritz_vectors, 1)

    ! extract the diagonals of the matrices
    diag_matrix = extract_diagonal_free(fun_matrix_gemv,dim_matrix)
    diag_second_matrix = extract_diagonal_free(fun_second_matrix_gemv,dim_matrix)
    
    ! 1. Variables initialization
    ! Select the initial ortogonal subspace based on lowest elements
    ! of the diagonal of the matrix
    copy_d = diag_matrix
    V = generate_preconditioner(copy_d, initial_dimension) ! Initial orthonormal basis
    
    ! Outer loop block Davidson schema
    outer_loop: do i=1, max_iterations

       ! 2. Generate subspace matrix problem by projecting into V
       matrixV = fun_matrix_gemv(V)
       second_matrixV = fun_second_matrix_gemv(V)
       matrix_proj = lapack_matmul('T', 'N', V, matrixV)
       second_matrix_proj = lapack_matmul('T', 'N', V, second_matrixV)

       ! 3. compute the eigenvalues and their corresponding ritz_vectors
       ! for the projected matrix using lapack
       call check_deallocate_matrix(eigenvectors_sub)
       
       if (allocated(eigenvalues_sub)) then
          deallocate(eigenvalues_sub)
       end if
       
       allocate(eigenvalues_sub(size(matrix_proj, 1)))
       allocate(eigenvectors_sub(size(matrix_proj, 1), size(matrix_proj, 2)))
       
       call lapack_generalized_eigensolver(matrix_proj, eigenvalues_sub, eigenvectors_sub, second_matrix_proj)
       
       ! 4. Check for convergence
       ritz_vectors = lapack_matmul('N', 'N', V, eigenvectors_sub(:, :lowest))

       ! 4.1 Residue calculation
       !! Matrix with the the eigenvalues in the diagonal
       lambda = eye( size(V, 2), size(V, 2))
       do j= 1, size(V, 2)
          lambda(j, j)= eigenvalues_sub(j)
       enddo

       ! residues
       residues = lapack_matmul('N', 'N', second_matrixV, eigenvectors_sub)
       guess = lapack_matmul('N', 'N', residues, lambda)
       deallocate(residues)
       residues =  lapack_matmul('N', 'N', matrixV, eigenvectors_sub) - guess
       ! errors
       do j=1,lowest
          errors(j) = norm(residues(:, j))
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

          correction = compute_DPR_free(matrixV, eigenvalues_sub, residues, diag_matrix, diag_second_matrix)
          
          ! 6. Increase Basis size
          call concatenate(V, correction)
          
          ! 7. Orthogonalize basis
          call lapack_qr(V)
          
       else
          ! 6. Otherwise reduce the basis of the subspace to the current correction
          V = lapack_matmul('N', 'N', V, eigenvectors_sub(:, :initial_dimension))
       end if
       
    end do outer_loop
    
    !  8. Check convergence
    if (i > max_iterations / initial_dimension) then
       print *, "Warning: Algorithm did not converge!!"
    end if
    
    
    ! Select the lowest eigenvalues and their corresponding ritz_vectors
    ! They are sort in increasing order
    eigenvalues = eigenvalues_sub(:lowest)
    
    ! Free memory
    call check_deallocate_matrix(correction)
    deallocate(eigenvalues_sub, eigenvectors_sub, V, matrix_proj, matrixV, second_matrixV, guess, residues)
    
    ! free optional matrix
    call check_deallocate_matrix(second_matrix_proj)
    
  end subroutine generalized_eigensolver_free
  

function compute_DPR_free(matrixV, eigenvalues, residues, diag_matrix, diag_second_matrix) result(correction)

      !> compute the correction vector using the DPR method for a matrix free diagonalization
      !> See correction_methods submodule for the implementations
      !> \param[in] matrixV: projection matrix * V
      !> \param[in] V: Basis of the iteration subspace
      !> \param[in] eigenvalues: of the reduce problem
      !> \param[in] residues: error for each eigenvalue/eigenvector pair
      !> \return correction matrix
      
      real(dp), dimension(:), intent(in) :: eigenvalues
      real(dp), dimension(:, :), intent(in) :: matrixV, residues
      real(dp), dimension(:), intent(in) :: diag_matrix, diag_second_matrix

      ! local variables
      !real(dp), dimension(size(V, 1),1) :: vector
      real(dp), dimension(size(matrixV, 1), size(matrixV, 2)) :: correction
      integer :: ii, j

      do j=1, size(matrixV, 2)
         do ii=1,size(matrixV, 1)
            correction(ii, j) = residues(ii, j) / (eigenvalues(j) * diag_second_matrix(ii)  - diag_matrix(ii))
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
           !> \brief Function to compute the optional matrix on the fly
           !> \param[in] i column/row to compute from matrix
           !> \param vec column/row from matrix
           use numeric_kinds, only: dp
           real (dp), dimension(:,:), intent(in) :: input_vect
           real (dp), dimension(size(input_vect,1),size(input_vect,2)) :: output_vect

         end function fun_A_gemv
    end interface

    ! local variable
    integer :: ii
    real(dp), dimension(dim,1) :: tmp_array
    
    do ii = 1,dim
      tmp_array = 0.0_dp
      tmp_array(ii,1) = 1.0_dp
      tmp_array = fun_A_gemv(tmp_array)
      out(ii) = tmp_array(ii,1)
    end do

  end function extract_diagonal_free


  function free_matmul(fun, array) result (matrix)
    !> \brief perform a matrix-matrix multiplication by generating a matrix on the fly using `fun`
    !> \param[in] fun function to compute a matrix on the fly
    !> \param[in] array matrix to multiply with fun
    !> \return resulting matrix

    ! input/output
    implicit none
    real(dp), dimension(:, :), intent(in) :: array
    real(dp), dimension(size(array, 1), size(array, 2)) :: matrix

    interface
       function fun(i, dim) result(vec)
         !> \brief Fucntion to compute the matrix `matrix` on the fly
         !> \param[in] i column/row to compute from `matrix`
         !> \param vec column/row from matrix
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
          matrix(i, j) = dot_product(vec, array(:, j))
       end do
    end do
    !$OMP END PARALLEL DO
    
  end function free_matmul

    
  subroutine check_deallocate_matrix(matrix)
    !> deallocate a matrix if allocated
    real(dp), dimension(:, :), allocatable, intent(inout) ::  matrix
    
    if (allocated(matrix)) then
       deallocate(matrix)
    end if
    
  end subroutine check_deallocate_matrix  


end module davidson_free


module davidson

  use numeric_kinds, only: dp
  use lapack_wrapper, only: lapack_generalized_eigensolver, lapack_matmul, lapack_matrix_vector, &
       lapack_qr, lapack_solver
  use array_utils, only: concatenate, eye, norm
  use davidson_dense, only: generalized_eigensolver_dense
  use davidson_free, only: generalized_eigensolver_free
  implicit none

  !> \private
  private
  !> \public
  public :: generalized_eigensolver
  
  interface generalized_eigensolver
  !> \brief Solve a (general) eigenvalue problem using different types of Davidson algorithms.

  !> \param[in] matrix: Matrix to diagonalize
  !> \param[in, opt] second_matrix: Optional matrix for the general eigenvalue problem:
  !> \f$ matrix \lambda = V second_matrix \lambda \f$
     
  !> \param[out] eigenvalues Computed eigenvalues
  !> \param[out] ritz_vectors approximation to the eigenvectors
  !> \param[in] lowest Number of lowest eigenvalues/ritz_vectors to compute
  !> \param[in] method Method to compute the correction vector. Available
  !> methods are,
  !>    DPR: Diagonal-Preconditioned-Residue
  !>    GJD: Generalized Jacobi Davidson
  !> \param[in] max_iterations: Maximum number of iterations
  !> \param[in] tolerance norm-2 error of the eigenvalues
  !> \param[in] method: Method to compute the correction vectors
  !> \param[in, opt] max_dim_sub: maximum dimension of the subspace search   
  !> \param[out] iters: Number of iterations until convergence
  !> \return eigenvalues and ritz_vectors of the matrix `matrix`

     procedure generalized_eigensolver_dense
     procedure generalized_eigensolver_free
     
  end interface generalized_eigensolver

end module davidson  


submodule (davidson_dense) correction_methods_generalized_dense
  !> submodule containing the implementations of different kind
  !> algorithms to compute the correction vectors for the Davidson's diagonalization

  implicit none
  
contains

  module function compute_correction_generalized_dense(matrix, eigenvalues, ritz_vectors, residues, method, second_matrix) &
       result(correction)
    !> see interface in davidson module
    real(dp), dimension(:), intent(in) :: eigenvalues
    real(dp), dimension(:, :), intent(in) :: matrix, ritz_vectors, residues
    real(dp), dimension(:, :), intent(in), optional :: second_matrix
    character(len=*), optional,intent(in) :: method
    logical :: gev 

    ! local variables
    character(len=10) :: opt 
    real(dp), dimension(size(matrix, 1), size(residues, 2)) :: correction

    !check optional arguments
    gev = present(second_matrix)
    opt="DPR"
    if (present(method)) opt=trim(method)
    
    select case (method)
    case ("DPR")
      if(gev) then
       correction = compute_DPR_generalized_dense(matrix, eigenvalues, residues, second_matrix)
      else
        correction = compute_DPR_generalized_dense(matrix, eigenvalues, residues)
      end if
    case ("GJD")
      if(gev) then
       correction = compute_GJD_generalized_dense(matrix, eigenvalues, ritz_vectors, residues, second_matrix)
      else
        correction = compute_GJD_generalized_dense(matrix, eigenvalues, ritz_vectors, residues)
      end if
    end select
    
  end function compute_correction_generalized_dense

  function compute_DPR_generalized_dense(matrix, eigenvalues, residues, second_matrix) result(correction)
    !> compute Diagonal-Preconditioned-Residue (DPR) correction
    real(dp), dimension(:), intent(in) :: eigenvalues
    real(dp), dimension(:, :), intent(in) :: matrix, residues
    real(dp), dimension(:, :), intent(in), optional ::  second_matrix 
    real(dp), dimension(size(matrix, 1), size(residues, 2)) :: correction
    
    ! local variables
    integer :: ii,j, m
    logical :: gev

    ! shape of matrix
    m = size(matrix, 1)
    gev = (present(second_matrix))

    do j=1, size(residues, 2)
       do ii=1,size(correction,1)
          if (gev) then
             correction(ii, j) = residues(ii, j) / (eigenvalues(j) * second_matrix(ii,ii) - matrix(ii, ii))
           else
              correction(ii, j) = residues(ii, j) / (eigenvalues(j)  - matrix(ii, ii))
           endif
        end do
    end do

  end function compute_DPR_generalized_dense

  function compute_GJD_generalized_dense(matrix, eigenvalues, ritz_vectors, residues, second_matrix) result(correction)
    !> Compute the Generalized Jacobi Davidson (GJD) correction
    
    real(dp), dimension(:), intent(in) :: eigenvalues
    real(dp), dimension(:, :), intent(in) :: matrix, ritz_vectors, residues
    real(dp), dimension(:, :), intent(in), optional :: second_matrix
    real(dp), dimension(size(matrix, 1), size(ritz_vectors, 2)) :: correction

    ! local variables
    integer :: k, m
    logical :: gev
    real(dp), dimension(size(matrix, 1), 1) :: rs
    real(dp), dimension(size(matrix, 1), size(matrix, 2)) :: arr, xs, ys
    real(dp), dimension(size(matrix, 1), 1) :: brr

    ! Diagonal matrix
    m = size(matrix, 1)
    gev = present(second_matrix)

    do k=1, size(ritz_vectors, 2)
       rs(:, 1) = ritz_vectors(:, k)
       xs = eye(m, m) - lapack_matmul('N', 'T', rs, rs)
       if(gev) then
        ys = matrix - eigenvalues(k) * second_matrix
       else
         ys = substract_from_diagonal(matrix, eigenvalues(k))
       end if
       arr = lapack_matmul('N', 'N', xs, lapack_matmul('N', 'N', ys, xs))
       brr = - reshape(residues(:, k), (/m, 1/))
       
       call lapack_solver(arr, brr)
       correction(:, k) = brr(:, 1)
    end do
    
  end function compute_GJD_generalized_dense

  function substract_from_diagonal(matrix, alpha) result(arr)
    !> susbstract an scalar from the diagonal of a matrix
    !> \param matrix: square matrix
    !> \param alpha: scalar to substract
    real(dp), dimension(:, :), intent(in) :: matrix
    real(dp), dimension(size(matrix, 1), size(matrix, 2)) :: arr
    real(dp), intent(in) :: alpha
    integer :: i

    arr = matrix
    do i=1,size(matrix, 1)
       arr(i, i) = arr(i, i) - alpha
    end do
    
  end function substract_from_diagonal
  
end submodule correction_methods_generalized_dense
