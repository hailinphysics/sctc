import numpy as np
import scipy.stats as ss
from scipy.sparse.linalg import eigs


def complexity_index(mcg):
    """
    Calculates the complexity index based on input gene expression matrix `mcg`.

    Args:
        mcg (numpy.ndarray): The input 2D gene expression matrix.

    Returns:
        tuple: A tuple containing two arrays:
            - cci (numpy.ndarray): Cell Complexity Index (CCI).
            - gci (numpy.ndarray): Gene Complexity Index (CCI).
    """
    mgc = np.transpose(mcg)

    kc0 = np.sum(mcg, axis=1)
    kg0 = np.sum(mcg, axis=0)

    # Compute MCC matrix
    mcc = np.dot(mcg / np.transpose([kc0]), mgc / np.transpose([kg0]))

    # Compute eigenvalues and eigenvectors of MCC
    e_val, e_vec = eigs(mcc)

    # Extract the second eigenvector (real part) as CCI
    cci = np.real(e_vec[:, 1])

    # Determine the direction of CCI
    scc = ss.spearmanr(kc0, cci)[0]
    if scc < 0:
        cci = -cci

    # Compute the GCI based on the CCI
    gci = np.dot(mgc, np.transpose([cci])).flatten() / kg0

    # Normalize CCI and GCI to [0, 1]
    cci = (cci - np.min(cci)) / (np.max(cci) - np.min(cci))
    gci = (gci - np.min(gci)) / (np.max(gci) - np.min(gci))

    return cci, gci
