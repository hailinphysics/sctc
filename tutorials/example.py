import scanpy as sc
import numpy as np
import igraph as ig
import sctc
import matplotlib.pyplot as plt

# Load scRNA-seq data and ensure the count matrix is a numpy array
adata = sc.read_h5ad('./data/hnd.h5ad')
if not isinstance(adata.X, np.ndarray):
    adata.X = adata.X.toarray()

# Preprocessing
sc.pp.filter_cells(adata, min_genes=1)
sc.pp.filter_genes(adata, min_cells=1)

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# Calculate Cell Complexity Index (CCI) and Gene Complexity Index (GCI)
cci, gci = sctc.complexity_index(adata.X)

# Calculate Nth-order complexity of cells (kcn) and genes (kgn)
kcn, gcn = sctc.complexity_order(adata.X)

# Choosing a subset of cells for visualization of cell ranking
n_choice = 50
choice = np.random.choice(kcn.shape[1], n_choice, replace=False)
choice_kc = kcn[:, choice]
choice_cci = cci[choice]

# Converting complexity scores to rankings
n_list = list(range(0, 16, 2))
complexity_array = choice_kc[n_list]
complexity_array = np.vstack([complexity_array, choice_cci])
complexity_array = np.transpose(complexity_array)
ranking_lists = sctc.convert_to_ranking(complexity_array)

# Creating a color map
cmap = plt.cm.get_cmap('Spectral', len(set(adata.obs['Day'])))
cmap = [cmap(i) for i in range(len(set(adata.obs['Day'])))]
cmap = dict(zip(adata.obs['Day'].unique(), cmap))
colors = [cmap[day] for day in adata.obs['Day'][choice]]

# Plotting the ranking
fig, ax = sctc.ranking_plot(ranking_lists, colors, marker_size=20, line_width=2)
fig.set_size_inches(7, 19)

xticks = []
xticks.extend([str(i) for i in n_list])
xticks.append('CCI')
ax.set_xticks(range(len(xticks)))
ax.set_xticklabels(xticks)
ax.set_xlabel('N')
ax.set_ylabel('Index')
plt.savefig('./results/ranking.png')

# Selecting a subset of genes for visualization of gene space
adata = sc.read_h5ad('./data/hnd.h5ad')
if not isinstance(adata.X, np.ndarray):
    adata.X = adata.X.toarray()

sc.pp.filter_cells(adata, min_genes=1)
sc.pp.filter_genes(adata, min_cells=1)

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

selected_indices = np.random.choice(adata.n_vars, 1500, replace=False)
adata = adata[:, selected_indices].copy()

# Creating a gene space network
gene_proxim = sctc.gene_proximity(adata.X)
gene_space = sctc.gene_space(gene_proxim)
layout = gene_space.layout_fruchterman_reingold()
gene_space.vs['size'] = 8
gene_space.vs['color'] = 'green'

# Plotting the gene space network
plot = ig.plot(gene_space, layout=layout)
plot.save('./results/gene_space.png')
