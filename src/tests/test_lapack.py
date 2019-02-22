from subprocess import (PIPE, Popen)
import argparse
import fnmatch
import numpy as np
from scipy import linalg
import os

msg = "test_lapack.py -i executable"

parser = argparse.ArgumentParser(description=msg)
parser.add_argument('-i', required=True, help="program to run")


def check_with_numpy(files: list) -> None:
    """
    Check that the calls to Lapack are the same that the ones get from Numpy.
    """
    names = ["test_lapack_{}.txt".format(x) for x in ("eigenvalues", "eigenvectors", "matrix")]
    names_generalized = ["test_lapack_eigenvalues_gen.txt", "test_lapack_eigenvectors_gen.txt",
                         "test_lapack_matrix.txt", "test_lapack_stx.txt"]
    print("STANDARD EIGENVALUE PROBLEM")
    check_lapack_eigensolver(names)
    print("GENERALIZED EIGENVALUE PROBLEM")
    check_lapack_eigensolver(names_generalized, True)

    # test QR
    check_lapack_qr("test_lapack_qr.txt", "test_lapack_matrix.txt")


def check_lapack_eigensolver(names: list, generalized: bool = False) -> None:
    """
    Compare the library calls of lapack to numpy (They are both wrappers over lapack)
    """
    if generalized:
        es_wrapper, vs_wrapper, mtx, stx = [np.loadtxt(f) for f in names]
    else:
        es_wrapper, vs_wrapper, mtx = [np.loadtxt(f) for f in names]
    dim = int(np.sqrt(mtx.size))
    mtx = mtx.reshape(dim, dim)
    vs_wrapper = vs_wrapper.reshape(dim, dim)
    if generalized:
        stx = stx.reshape(dim, dim)
    else:
        stx = np.eye(dim)

    # compute eigenvalue/eigenvectors with numpy
    es_numpy, vs_numpy = linalg.eigh(mtx, b=stx)

    # compare the results
    assert np.allclose(es_wrapper, es_numpy)
    assert np.allclose(np.abs(vs_wrapper), np.abs(vs_numpy))
    print("Check that the call to Lapack produce the same value that Numpy: True")


def check_lapack_qr(file_qr: str, file_matrix: str) -> None:
    """
    Compare the QR subroutine to the one in Numpy
    """
    qr, mtx = [np.loadtxt(x) for x in [file_qr, file_matrix]]
    dim = int(np.sqrt(mtx.size))
    mtx = mtx.reshape(dim, dim)
    qr = qr.reshape(dim, dim)

    q_numpy, _r = np.linalg.qr(mtx)

    np.allclose(qr, q_numpy)
    print("COMPARE THE QR SUBROUTINE TO THE ONE IN NUMPY")


def main():
    executable = read_cmd_line()
    print("path to executable: ", executable)

    p = Popen(executable, stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=True)
    rs = p.communicate()
    err = rs[1]
    if err:
        raise RuntimeError("Submission Errors: {}".format(err))
    else:
        files = fnmatch.filter(os.listdir('.'), "test_lapack_*.txt")
        check_with_numpy(files)


def read_cmd_line():
    """
    Read the file to run
    """
    args = parser.parse_args()
    return args.i


if __name__ == "__main__":
    main()
