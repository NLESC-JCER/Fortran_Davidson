module benchmark

  use numeric_kinds, only: dp
  use Davidson, only: eigensolver, generate_diagonal_dominant

  implicit none

contains
  
  subroutine compute_benchmark(dims, lowest, times)
    !> Benchmark the Davidson algorithm
    !> \param dims: Vector of integer with the dimension of the matrices to test
    !> \param lowest: number of eigenvalues/eigenvectors to compute
    !> \return times: resulting times
    
    integer, dimension(:), intent(in) :: dims
    integer, intent(in) :: lowest
    real(dp), dimension(size(dims), 2), intent(out) :: times

    ! local variable
    integer :: i, m
    real(dp), dimension(:), allocatable :: eigenvalues
    real(dp), dimension(:, :), allocatable :: mtx, eigenvectors
    real(dp) ::  dt


    do i=1, size(dims)

       call benchmark_method(mtx, eigenvalues, eigenvectors, "DPR", dims(i), lowest, dt)
       times(i, 1) = dt
       ! call benchmark_method(mtx, eigenvalues, eigenvectors, "GJD", dims(i), lowest, dt)
       ! times(i, 2) = dt
    end do
    
    deallocate(mtx)
    
  end subroutine compute_benchmark

  subroutine benchmark_method(mtx, eigenvalues, eigenvectors, method, dim, lowest, dt)
    !> Benchmark different correction methods
    real(dp), dimension(:), allocatable :: eigenvalues
    real(dp), dimension(:, :), allocatable :: mtx, eigenvectors
    real(dp), intent(out) :: dt
    integer, intent(in) :: dim, lowest
    character(len=*), intent(in) :: method
    real(dp) :: t1, t2

    call free_space(mtx, eigenvalues, eigenvectors)
    ! generate test matrix
    mtx = generate_diagonal_dominant(dim,1d-4)
    ! diagonalize
    allocate(eigenvalues(lowest))
    allocate(eigenvectors(dim, lowest))
    print *, "Davidson method: ", method, " dimension: ", dim
    call cpu_time(t1)
    call eigensolver(mtx, eigenvalues, eigenvectors, 3, method, 1000, 1d-8)
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
