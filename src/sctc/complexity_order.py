import numpy as np


def complexity_order(mcg, nmax=30):
    """
    Calculates the N-order complexity of a network.

    Args:
        mcg (numpy.ndarray): The input 2D gene expression matrix.
        nmax (int, optional): Maximum number of iterations. Defaults to 30.

    Returns:
        tuple: Nth-order cell complexity (kc) and Nth-order gene complexity (kg).
               kc (kg) is 2D array where rows represent cells (genes) and columns correspond to complexity order N.
    """
    # Convert input matrix to float64
    mcg = mcg.astype(np.float64)
    mgc = np.transpose(mcg)

    kc = []
    kg = []

    # Calculate 0th-order complexity
    kc.append(np.sum(mcg, axis=1))
    kg.append(np.sum(mcg, axis=0))

    for n in range(1, nmax):
        # Calculate Nth-order complexity iteratively
        kc.append(np.dot(mcg, np.transpose([kg[n - 1]])).flatten() / kc[0])
        kg.append(np.dot(mgc, np.transpose([kc[n - 1]])).flatten() / kg[0])

    # Normalize kc and kg
    kc = np.array(kc)
    kg = np.array(kg)
    max_kc = np.max(kc, axis=1).reshape(kc.shape[0], 1)
    min_kc = np.min(kc, axis=1).reshape(kc.shape[0], 1)
    kc = (kc - min_kc) / (max_kc - min_kc)
    max_kg = np.max(kg, axis=1).reshape(kg.shape[0], 1)
    min_kg = np.min(kg, axis=1).reshape(kg.shape[0], 1)
    kg = (kg - min_kg) / (max_kg - min_kg)

    return kc, kg
