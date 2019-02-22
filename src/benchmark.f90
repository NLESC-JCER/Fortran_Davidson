module benchmark

  use numeric_kinds, only: dp
  use davidson, only: generalized_eigensolver
  use array_utils, only: generate_diagonal_dominant

  implicit none

contains
  
  subroutine compute_benchmark(dims, lowest, sparsity, times, iters)
    !> Benchmark the Davidson algorithm
    !> \param dims: Vector of integer with the dimension of the matrices to test
    !> \param lowest: number of eigenvalues/eigenvectors to compute
    !> \param sparsity: magnitud of the off-diagonal terms
    !> \return times: resulting times
    !> \param number of iterations
    
    integer, dimension(:), intent(in) :: dims
    integer, intent(in) :: lowest
    integer, dimension(size(dims), 2), intent(out) :: iters
    real(dp), dimension(size(dims), 2), intent(out) :: times
    real(dp) :: sparsity
    
    ! local variable
    integer :: i, iter_i
    real(dp), dimension(:), allocatable :: eigenvalues
    real(dp), dimension(:, :), allocatable :: mtx, eigenvectors
    real(dp) ::  dt


    do i=1, size(dims)       
       call benchmark_method(mtx, eigenvalues, eigenvectors, "DPR", dims(i), lowest, sparsity, dt, iter_i)
       print *, "cycles: ", iter_i
       iters(i, 1) = iter_i
       times(i, 1) = dt
       
       call benchmark_method(mtx, eigenvalues, eigenvectors, "GJD", dims(i), lowest, sparsity, dt, iter_i)
       print *, "cycles: ", iter_i
       iters(i, 2) = iter_i
       times(i, 2) = dt
       
    end do

    call free_space(mtx, eigenvalues, eigenvectors)
    
  end subroutine compute_benchmark

  subroutine benchmark_method(mtx, eigenvalues, eigenvectors, method, dim, lowest, sparsity, dt, iter_i)
    !> Benchmark different correction methods
    real(dp), dimension(:), allocatable :: eigenvalues
    real(dp), dimension(:, :), allocatable :: mtx, eigenvectors
    real(dp), intent(out) :: dt
    integer, intent(in) :: dim, lowest
    integer, intent(out) :: iter_i
    character(len=*), intent(in) :: method
    real(dp) :: t1, t2, sparsity

    call free_space(mtx, eigenvalues, eigenvectors)
    ! generate test matrix
    mtx = generate_diagonal_dominant(dim, sparsity)
    ! diagonalize
    allocate(eigenvalues(lowest))
    allocate(eigenvectors(dim, lowest))
    print *, "Davidson method: ", method, " dimension: ", dim
    call cpu_time(t1)
    call generalized_eigensolver(mtx, eigenvalues, eigenvectors, 3, method, 1000, 1d-8, iter_i)
    call cpu_time(t2)
    dt = t2 - t1
    print *, "time: ", dt

    
  end subroutine benchmark_method
  
  subroutine free_space(mtx, es, vs)
    !> Deallocate variables
    real(dp), dimension(:, :), allocatable, intent(inout) :: mtx, vs
    real(dp), dimension(:), allocatable, intent(inout) :: es

    if (allocated(mtx)) then
       deallocate(mtx)
    end if

    if (allocated(vs)) then
       deallocate(vs)
    end if

    if (allocated(es)) then
       deallocate(es)
    end if

  end subroutine free_space
    
    
end module benchmark
