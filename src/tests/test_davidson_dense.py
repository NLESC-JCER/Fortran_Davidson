from subprocess import (PIPE, Popen)
import argparse
import fnmatch
import numpy as np
from scipy import linalg
import os

msg = "test_davidson_dense.py -i executable"

parser = argparse.ArgumentParser(description=msg)
parser.add_argument('-i', required=True, help="program to run")


def check_eigenvalues(files, generalized: bool = False):
    """
    Check that the eigenvalues/eigenvectors are the same that the ones
    computed with numpy
    """
    files.sort()
    if generalized:
        es_DPR, es_GJD, vs_DPR, vs_GJD, mtx, stx = [np.loadtxt(x) for x in files]
    else:
        es_DPR, es_GJD, vs_DPR, vs_GJD, mtx = [np.loadtxt(x) for x in files]
    dim = int(np.sqrt(mtx.size))
    mtx = mtx.reshape(dim, dim)
    ncols = es_DPR.size
    vs_DPR = vs_DPR.reshape(dim, ncols)
    vs_GJD = vs_GJD.reshape(dim, ncols)

    if generalized:
        stx = stx.reshape(dim, dim)
    else:
        stx = np.eye(dim)
    # compute the eigenvalues/eigenvectors of the test matrix
    es, vs = linalg.eigh(mtx, b=stx)

    # eigenvalues are the same
    assert np.allclose(es_DPR, es_GJD)
    assert np.allclose(es[:ncols], es_DPR)

    # Precision Eigenvectos numpy
    for i in range(ncols):
        print("precision eigenvalue ", i, ":")
        x = np.linalg.norm(np.dot(mtx, vs[:, i]) - (es[i] * np.dot(stx, vs[:, i])))
        y = np.linalg.norm(np.dot(mtx, vs_DPR[:, i]) - (es_DPR[i] * np.dot(stx, vs_DPR[:, i])))
        z = np.linalg.norm(np.dot(mtx, vs_GJD[:, i]) - (es_GJD[i] * np.dot(stx, vs_GJD[:, i])))

        msg = "\tnumpy: {:5.2e} DPR: {:5.2e} GJD: {:5.2e}".format(x, y, z)
        print(msg)


def main():
    executable = read_cmd_line()
    print("path to executable: ", executable)

    p = Popen(executable, stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=True)
    rs = p.communicate()
    err = rs[1]
    if err:
        raise RuntimeError("Submission Errors: {}".format(err))
    else:
        files = fnmatch.filter(os.listdir('.'), "test_dense_spec_*.txt")
        print("testing normal eigenvalue solver")
        check_eigenvalues(files)

        # generalized case
        files_generalized = fnmatch.filter(os.listdir('.'), "test_dense_gen_*.txt")
        print("testing generalized eigenvalue solver")
        check_eigenvalues(files_generalized, True)


def read_cmd_line():
    """
    Read the file to run
    """
    args = parser.parse_args()
    return args.i


if __name__ == "__main__":
    main()
