import anndata as ad
import pytest
import scanpy as sc
from scipy import sparse
import pyucell


def test_ranks_from_anndata():
    adata = sc.datasets.pbmc3k()
    ranks = pyucell.get_rankings(adata, max_rank=500)
    return ranks


def test_ranks_from_matrix():        
    X = sparse.random(1000, 20000, density=0.1, format='csr')
    ranks = pyucell.get_rankings(X, max_rank=200)
    return ranks

def test_compute_ucell():        
    adata = sc.datasets.pbmc3k()
    signatures = {
        'T_cell': ['CD3D', 'CD3E', 'CD2'],
        'B_cell': ['MS4A1', 'CD79A', 'CD79B']
    }
    compute_ucell_scores(adata, signatures=signatures, chunk_size=500)
    return adata   