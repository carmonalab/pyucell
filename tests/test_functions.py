import anndata as ad
import numpy as np
import pytest
from scipy import sparse

@pytest.fixture
def ranks_from_anndata():
    adata = sc.datasets.pbmc3k()
    ranks2 = pyucell.tl.get_rankings(adata, max_rank=500)

@pytest.fixture
def ranks_from_matrix():        
    X = sparse.random(1000, 20000, density=0.1, format='csr')
    ranks = pyucell.tl.get_rankings(X, max_rank=200)