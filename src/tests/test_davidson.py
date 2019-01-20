from subprocess import (PIPE, Popen)
import argparse
import fnmatch
import numpy as np
import os

msg = "namd.py -i input"

parser = argparse.ArgumentParser(description=msg)
parser.add_argument('-i', required=True, help="program to run")


def check_eigenvalues(files):
    """
    Check that the eigenvalues/eigenvectors are the same that the ones
    computed with numpy
    """
    files.sort()
    es_DPR, es_GJD, vs_DPR, vs_GJD, mtx = [np.loadtxt(x) for x in files]
    dim = int(np.sqrt(mtx.size))
    mtx = mtx.reshape(dim, dim)
    ncols = es_DPR.size
    vs_DPR = vs_DPR.reshape(dim, ncols)

    # compute the eigenvalues/eigenvectors of the test matrix
    es, vs = np.linalg.eigh(mtx)

    # eigenvalues are the same
    assert np.allclose(es_DPR, es_GJD)
    assert np.allclose(es[:ncols], es_DPR)


def main():
    executable = read_cmd_line()
    print("path to executable: ", executable)

    p = Popen(executable, stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=True)
    rs = p.communicate()
    err = rs[1]
    if err:
        raise RuntimeError("Submission Errors: {}".format(err))
    else:
        files = fnmatch.filter(os.listdir('.'), "*.txt")
        check_eigenvalues(files)


def read_cmd_line():
    """
    Read the file to run
    """
    args = parser.parse_args()
    return args.i


if __name__ == "__main__":
    main()
