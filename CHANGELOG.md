
# Version 0.0.2

## Added
 * Algorithm for the [Generalized eigenvalue problem](https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix#Generalized_eigenvalue_problem) using the Davidson method.
 * Matrix free implementation of the Davidson method using the **DPR** correction.
 * Tests for the lapack wrappers.
 
## Changed
 * Append the orthonormalized correction vector to the previous subspace by computing only the
 matrix elements involved with the new correction vector.
 
 * Used as maximum dimension of the subspace 10 times the number of requested eigenvalues.

# Version 0.0.1

## Dependencies
 
 * Fortran compiler >= 6.0 (supporting `submodules` and `move_alloc`)
 * [lapack](http://www.netlib.org/lapack/) or [MKL](https://software.intel.com/en-us/mkl).
 * [Cmake](https://cmake.org/) >= 3.0
