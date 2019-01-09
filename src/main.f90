program main
  use numeric_kinds, only: sp, dp
  
  implicit none

  ! real(dp) :: matrix(4, 4)=0_dp
  ! real(dp) :: diag(4)
  ! integer :: i
  ! integer :: error !error flag
  ! character(len=*), parameter :: path_file = "tests/test_files/hdf5_test.h5"
  ! character(len=*), parameter :: dataset1 = "array"
  ! character(len=*), parameter :: dataset2 = "group0/vector"
  ! real(dp), dimension(4, 4) :: buffer1
  ! real(dp), dimension(10) :: buffer2
  ! real(dp), dimension(10) :: random_values


  ! ! generate a diagonal matrix
  ! forall(i = 1: 4) matrix(i, i) = i ** 2
  ! print *, "diagonal: ", diagonal(matrix)
  ! print *, "trace: ", trace(matrix)
  ! print *, "rank: ", rank(matrix)

  ! ! read data from a hdf5 file
  ! call retrieve_data_in_hdf5(path_file, dataset1, buffer1)
  ! call retrieve_data_in_hdf5(path_file, dataset2, buffer2)
  ! print *, sum(buffer1)
  ! print *, dot_product(buffer2, buffer2)

  
end program main
