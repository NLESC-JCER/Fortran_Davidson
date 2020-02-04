# Version 0.1.1 [04/02/2020]

### Changed

* Improved variables names.
* Removed unnecessary ritz_vectors calculations.
* Optimized the DPR for the free version.
* Removed recomputation of the residues.

# Version 0.1.0

### Changed

* Update the whole projection matrix after adding some correction vector. Replace the block update schema.

### Fixed

* Fixed [several bugs](https://github.com/NLESC-JCER/Fortran_Davidson/issues/29) in the matrix free implementation.

# Version 0.0.5

### New

 * Select the initial orthonormal basis set based on the lowest diagonal elements of the matrix

# Version 0.0.4

## changed
 *  For the Matrix-free implementation, the function that computes the target matrix on the fly,
 was replaced by another function representing the action of a matrix-operator over a vector

# Version 0.0.3

## Changed
 * split the `dense` and `matrix` free into modules
 * Moved the `lapack` calls to their own module
 * Moved the array utilities to their own module
 
## Deleted
 * Removed the `DPR` vs `GJD` benchmark

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
