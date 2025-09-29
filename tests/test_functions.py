import scanpy as sc
from scipy import sparse

import pyucell

adata = sc.datasets.pbmc3k()
signatures = {"Tcell": ["CD3D", "CD3E", "CD2"], "Bcell": ["MS4A1", "CD79A", "CD79B"]}


def signature_columns_exist(adata, signatures, suffix="_UCell"):
    missing_cols = [f"{sig}{suffix}" for sig in signatures if f"{sig}{suffix}" not in adata.obs.columns]
    assert not missing_cols, f"Missing columns in adata.obs: {missing_cols}"


def test_ranks_from_anndata():
    ranks = pyucell.get_rankings(adata)
    assert isinstance(ranks, sparse.spmatrix)


def test_ranks_from_matrix():
    X = sparse.random(1000, 20000, density=0.1, format="csr")
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
    signatures_neg = {"Tcell": ["CD3D+", "CD3E+", "CD2+", "LYZ-"], "Bcell": ["MS4A1+", "CD79A+", "CD2-"]}
    pyucell.compute_ucell_scores(adata, signatures=signatures_neg)
    signature_columns_exist(adata, signatures_neg)


def test_missing_genes():
    signatures_miss = {"Tcell": ["CD3D", "CD3E", "CD2"], "Bcell": ["MS4A1", "CD79A", "notagene"]}
    pyucell.compute_ucell_scores(adata, signatures=signatures_miss)
    signature_columns_exist(adata, signatures_miss)


def all_missing():
    signatures_miss = {"Tcell": ["CD3D", "CD3E", "CD2"], "Bcell": ["notagene1", "notagene2"]}
    pyucell.compute_ucell_scores(adata, signatures=signatures_miss)
    signature_columns_exist(adata, signatures_miss)


def layers():
    adata.layers["newlayer"] = adata.X.copy()
    pyucell.compute_ucell_scores(adata, signatures=signatures, layer="newlayer")
    signature_columns_exist(adata, signatures)


def knn_basic():
    sc.pp.log1p(adata)
    sc.tl.pca(adata, svd_solver="arpack", n_comps=10)

    suffix1 = "_UCell"
    obs_cols = [s + suffix1 for s in signatures.keys()]
    pyucell.compute_ucell_scores(adata, signatures=signatures, suffix=suffix1)

    suffix2 = "_kNN"
    knn_cols = [s + suffix2 for s in obs_cols]
    pyucell.smooth_knn_scores(adata, obs_columns=obs_cols, suffix=suffix2)

    signature_columns_exist(adata, knn_cols)


def knn_from_graph():
    sc.pp.log1p(adata)
    sc.tl.pca(adata, svd_solver="arpack", n_comps=10)
    sc.pp.neighbors(adata, n_neighbors=10, use_rep="X_pca", key_added="customgraph")

    suffix1 = "_UCell"
    obs_cols = [s + suffix1 for s in signatures.keys()]
    pyucell.compute_ucell_scores(adata, signatures=signatures, suffix=suffix1)

    suffix2 = "_kNN"
    knn_cols = [s + suffix2 for s in obs_cols]
    pyucell.smooth_knn_scores(adata, obs_columns=obs_cols, graph_key="customgraph_connectivities")

    signature_columns_exist(adata, knn_cols)
