import numpy as np
import igraph as ig
import scipy.sparse as sp

def gene_proximity(mcg, rca_th=1):
    """
    Calculates gene proximity based on the gene expression matrix.

    Args:
        mcg (np.ndarray): The input 2D gene expression matrix.
        rca_th (float, optional): Threshold for Revealed Comparative Advantage (RCA). Defaults to 1.

    Returns:
        np.ndarray: Gene proximity matrix.
    """
    # Calculate RCA matrix
    mcg_c = np.sum(mcg, axis=1)
    mcg_c = mcg_c.reshape(-1, 1)
    mcg_g = np.sum(mcg, axis=0)
    mcg_s = np.sum(mcg)

    rca = mcg / mcg_c / mcg_g * mcg_s

    # Threshold RCA values
    rca[rca < rca_th] = 0
    rca[rca >= rca_th] = 1

    # Calculate gene proximity matrix base on RCA matrix
    n_cell, n_gene = mcg.shape
    gene_proxim = np.dot(rca.transpose(), rca)

    for i in range(n_gene - 1):
        for j in range(i + 1, n_gene):
            gene_proxim[i, j] /= max(gene_proxim[i, i], gene_proxim[j, j])
            gene_proxim[j, i] = gene_proxim[i, j]
    np.fill_diagonal(gene_proxim, 0)

    return gene_proxim


def gene_space(gene_proximity):
    """
    Constructs a gene space network using a spanning tree approach.

    Args:
        gene_proximity (np.ndarray): The input 2D Gene proximity matrix.

    Returns:
        ig.Graph: Gene space network.
    """
    # Calculate gene proximity matrix which server as the adjacency matrix of the gene space network.
    gene_adjmat = sp.csc_matrix(gene_proximity)

    sources, targets = gene_adjmat.nonzero()
    edgelist = list(zip(sources.tolist(), targets.tolist()))
    g = ig.Graph(edgelist, edge_attrs={'weight': gene_adjmat.data.tolist()})
    g.vs['id'] = g.vs.indices
    g.to_undirected()

    # Calculate the maximum spanning tree as the gene space network
    inv_weight = [1. / w for w in g.es['weight']]
    gene_space = g.spanning_tree(weights=inv_weight)

    return gene_space

