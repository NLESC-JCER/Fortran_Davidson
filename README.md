[![Build Status](https://travis-ci.org/NLESC-JCER/Fortran_Davidson.svg?branch=master)](https://travis-ci.org/NLESC-JCER/Fortran_Davidson)

Davidson Eigensolver
===================
This package contains a Modern Fortran implementation of the *Davidson diagonalization algorithms*.
Different schemas are available to compute the correction.

Available correction methods are:
 * **DPR**: Diagonal-Preconditioned-Residue
 * **GJD**: Generalized Jacobi Davidson


### Note:
The Davidson method is suitable for **diagonal-dominant symmetric matrices**, that are quite common
in certain scientific problems like [electronic structure](https://en.wikipedia.org/wiki/Electronic_structure). The Davidson method could be not practical
for other kind of symmetric matrice.

Usage
-----
The following program call the `eigensolver` *subroutine* from the `davidson` module and computes
the lowest 3 eigenvalues and corresponding eigenvectors, using the *GJD* method with a tolerance
of `1e-8` and `100` maximum iteration.
```fortran
program main
  use numeric_kinds, only: dp
  use davidson, only: eigensolver, generate_diagonal_dominant
 
  implicit none

  real(dp), dimension(20, 20) :: mtx
  real(dp), dimension(3) :: eigenvalues
  real(dp), dimension(20, 3) :: eigenvectors

  mtx = generate_diagonal_dominant(20, 1d-4)
  call eigensolver(mtx, eigenvalues, eigenvectors, 3, "GJD", 100, 1d-8)
  print *, eigenvalues
  print *, eigenvectors

end program main
```
The helper  `generate_diagonal_dominant` function generates a diagonal dominant
matrix with entries to the diagonal close to row number `(i=1, number_of_rows)`
and random number of the order `1e-4` on the off-diagonal entries.


### References:
 * [Davidson diagonalization method and its applications to electronic structure calculations](https://pdfs.semanticscholar.org/57811/eaf768d1a006f505dfe24f329874a679ba59.pdf?_ga=2.219777566.664950272.1547548596-1327556406.1547548596)
 * [Numerical Methods for Large Eigenvalue Problem](https://doi.org/10.1137/1.9781611970739)

Installation
------------

To compile execute:
```
cmake -H. -Bbuild
cmake --build build
```

Dependencies
------------
This packages assumes that you have installed the following packages:
 * A Fortran compiler >=  version 6.0 
 * [CMake](https://cmake.org/)
 * [Lapack](http://www.netlib.org/lapack/)
	
Optionally, If an [MKL](https://software.intel.com/en-us/mkl) library is available the package will try to find it.
