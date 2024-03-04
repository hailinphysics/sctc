# SCTC: Single-Cell Transcriptome Complexity

The **SCTC** Python package provides tools for calculating Single-Cell Transcriptome Complexity (SCTC) from scRNA-seq data. This metric helps to characterize cell developmental potential and infer single-cell pseudotime.

## Installation

You can install the SCTC package using `pip`:

```bash
pip install sctc
```

## Functions

### 1. Complexity Index

The `complexity_index` function computes the complexity index (CCI) and gene contribution index (GCI) for a given single-cell transcriptome dataset.

#### Usage

```python
import sctc

# Load your scRNA-seq gene expression matrix as a numpy array 'mcg'
# mcg should have cells as rows and genes as columns

cci, gci = sctc.complexity_index(mcg)
```

#### Output

- `cci`: CCI for each cell (normalized between 0 and 1).
- `gci`: GCI for each gene (normalized between 0 and 1).

### 2. Nth-complexity

The `complexity_order` function calculates the Nth-complexity of cells and genes.

#### Usage

```python
import sctc

# Load your scRNA-seq gene expression matrix as a numpy array 'mcg'
# mcg should have cells as rows and genes as columns

kc, kg = sctc.complexity_order(mcg, nmax=30)
```

#### Output

- `kc`: CCI for each cell (normalized between 0 and 1).
- `kg`: GCI for each gene (normalized between 0 and 1).

### 3. Complexity Ranking

The `convert_to_ranking` function convert a 2D complexity array into a 2D ranking array.\
The `ranking_plot` function visualize the rankings.

#### Usage

```python
import sctc

# Load a 2D array of complexities.

ranking_lists = convert_to_ranking(complexity_array)
fig, ax = ranking_plot(ranking_lists)
```

#### Output

- `fig`: matplotlib.figure.Figure.
- `ax`: matplotlib.axes._subplots.AxesSubplot.

### 4. Gene Space

The `gene_proximity` function calculates gene proximity based on the gene expression matrix.\
The `gene_space` function constructs a gene space network using a spanning tree approach.

#### Usage

```python
import sctc

# Load your scRNA-seq gene expression matrix as a numpy array 'mcg'
# mcg should have cells as rows and genes as columns

gene_proxim = sctc.gene_proximity(mcg)
gene_space = sctc.gene_space(gene_proxim)
```

#### Output

- `gene_space`: An igraph object representing the maximum spanning tree.

## Example

```python
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
N_list = list(range(0, 16, 2))
complexity_array = choice_kc[N_list]
complexity_array = np.vstack([complexity_array, choice_cci])
complexity_array = np.transpose(complexity_array)
ranking = sctc.convert_to_ranking(complexity_array)

# Creating a color map
cmap = plt.cm.get_cmap('Spectral', len(set(adata.obs['Day'])))
cmap = [cmap(i) for i in range(len(set(adata.obs['Day'])))]
cmap = dict(zip(adata.obs['Day'].unique(), cmap))
colors = [cmap[day] for day in adata.obs['Day'][choice]]

# Plotting the ranking
fig, ax = sctc.ranking_plot(ranking, colors, marker_size=20, line_width=2)
fig.set_size_inches(7, 19)

xticks = []
xticks.extend([str(i) for i in N_list])
xticks.append('CCI')
ax.set_xticks(range(len(xticks)))
ax.set_xticklabels(xticks)
ax.set_xlabel('N')
ax.set_ylabel('Index')
plt.savefig('./results/ranking.png')

# Selecting a subset of genes for visualization
selected_indices = np.random.choice(adata.n_vars, 1500, replace=False)
adata = adata[:, selected_indices]

# Creating a gene space network
gene_proxim = sctc.gene_proximity(adata.X)
gene_space = sctc.gene_space(gene_proxim)
layout = gene_space.layout_fruchterman_reingold()
gene_space.vs['size'] = 8
gene_space.vs['color'] = 'green'

# Plotting the gene space network
plot = ig.plot(gene_space, layout=layout)
plot.save('./results/gene_space.png')
```

## Tutorials

For more examples, please refer to the [tutorials](https://github.com/hailinphysics/sctc/tree/main/tutorials)
