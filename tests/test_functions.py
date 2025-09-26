import anndata as ad
import pytest
import scanpy as sc
from scipy import sparse
import pyucell

adata = sc.datasets.pbmc3k() 
signatures = {
    'T_cell': ['CD3D', 'CD3E', 'CD2'],
    'B_cell': ['MS4A1', 'CD79A', 'CD79B']
}

def signature_columns_exist(adata, signatures, suffix="_UCell"):
    missing_cols = [f"{sig}{suffix}" for sig in signatures if f"{sig}{suffix}" not in adata.obs.columns]
    assert not missing_cols, f"Missing columns in adata.obs: {missing_cols}"


def test_ranks_from_anndata():   
    ranks = pyucell.get_rankings(adata)
    assert isinstance(ranks, sparse.spmatrix)

def test_ranks_from_matrix():        
    X = sparse.random(1000, 20000, density=0.1, format='csr')
    ranks = pyucell.get_rankings(X, max_rank=500)
    assert isinstance(ranks, sparse.spmatrix)

def test_compute_ucell():          
    pyucell.compute_ucell_scores(adata, signatures=signatures)
    signature_columns_exist(adata, signatures)

def test_chunk():          
    pyucell.compute_ucell_scores(adata, signatures=signatures, chunk_size=100)
    signature_columns_exist(adata, signatures) 

def test_wneg():          
    pyucell.compute_ucell_scores(adata, signatures=signatures, w_neg=0.5)
    signature_columns_exist(adata, signatures) 

def skip_missing():          
    pyucell.compute_ucell_scores(adata, signatures=signatures, missing_genes="skip")
    signature_columns_exist(adata, signatures)

def test_serial():          
    pyucell.compute_ucell_scores(adata, signatures=signatures, n_jobs=1)
    signature_columns_exist(adata, signatures)


def test_neg_signatures():        
    signatures_neg = {
        'T_cell': ['CD3D+', 'CD3E+', 'CD2+', 'LYZ-'],
        'B_cell': ['MS4A1+', 'CD79A+', 'CD2-']
    }
    pyucell.compute_ucell_scores(adata, signatures=signatures_neg)
    signature_columns_exist(adata, signatures)  

def test_missing_genes():        
    signatures_miss = {
        'T_cell': ['CD3D', 'CD3E', 'CD2'],
        'B_cell': ['MS4A1', 'CD79A', 'notagene']
    }
    pyucell.compute_ucell_scores(adata, signatures=signatures_miss)
    signature_columns_exist(adata, signatures)  

def all_missing():        
    signatures_miss = {
        'T_cell': ['CD3D', 'CD3E', 'CD2'],
        'B_cell': ['notagene1','notagene2']
    }
    pyucell.compute_ucell_scores(adata, signatures=signatures_miss)
    signature_columns_exist(adata, signatures)  

def layers():        
    adata.layers["newlayer"] = adata.X.copy()
    pyucell.compute_ucell_scores(adata, signatures=signatures, layer="newlayer")
    signature_columns_exist(adata, signatures)



