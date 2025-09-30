import scanpy as sc
from scipy import sparse
import pytest
import pyucell

@pytest.fixture(scope="session")
def base_adata():
    adata = sc.datasets.pbmc3k()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    return adata

@pytest.fixture
def adata(base_adata):
    # each test gets a fresh copy
    return base_adata.copy()

@pytest.fixture(scope="session")
def signatures():
    signatures = {"Tcell": ["CD3D", "CD3E", "CD2"], "Bcell": ["MS4A1", "CD79A", "CD79B"]}
    return signatures

@pytest.fixture
def adata_with_scores(adata, signatures):    
    sc.tl.pca(adata, svd_solver="arpack", n_comps=10)
    suffix1 = "_UCell"
    suffix2 = "_kNN"
    obs_cols = [s + suffix1 for s in signatures.keys()]
    pyucell.compute_ucell_scores(adata, signatures=signatures, suffix=suffix1)
    return adata.copy()


def signature_columns_exist(adata, signatures, suffix="_UCell"):
    missing_cols = [f"{sig}{suffix}" for sig in signatures if f"{sig}{suffix}" not in adata.obs.columns]
    assert not missing_cols, f"Missing columns in adata.obs: {missing_cols}"


def test_ranks_from_anndata(adata):
    ranks = pyucell.get_rankings(adata)
    assert isinstance(ranks, sparse.spmatrix)


def test_ranks_from_matrix():
    X = sparse.random(1000, 20000, density=0.1, format="csr")
    ranks = pyucell.get_rankings(X, max_rank=500)
    assert isinstance(ranks, sparse.spmatrix)


def test_compute_ucell(adata, signatures):
    pyucell.compute_ucell_scores(adata, signatures=signatures)
    signature_columns_exist(adata, signatures)


def test_chunk(adata, signatures):
    pyucell.compute_ucell_scores(adata, signatures=signatures, chunk_size=100)
    signature_columns_exist(adata, signatures)


def test_wneg(adata, signatures):
    pyucell.compute_ucell_scores(adata, signatures=signatures, w_neg=0.5)
    signature_columns_exist(adata, signatures)


def test_skip_missing(adata, signatures):
    pyucell.compute_ucell_scores(adata, signatures=signatures, missing_genes="skip")
    signature_columns_exist(adata, signatures)


def test_serial(adata, signatures):
    pyucell.compute_ucell_scores(adata, signatures=signatures, n_jobs=1)
    signature_columns_exist(adata, signatures)


def test_neg_signatures(adata):
    signatures_neg = {"Tcell": ["CD3D+", "CD3E+", "CD2+", "LYZ-"], "Bcell": ["MS4A1+", "CD79A+", "CD2-"]}
    pyucell.compute_ucell_scores(adata, signatures=signatures_neg)
    signature_columns_exist(adata, signatures_neg)


def test_missing_genes(adata):
    signatures_miss = {"Tcell": ["CD3D", "CD3E", "CD2"], "Bcell": ["MS4A1", "CD79A", "notagene"]}
    pyucell.compute_ucell_scores(adata, signatures=signatures_miss)
    signature_columns_exist(adata, signatures_miss)


def test_all_missing(adata):
    signatures_miss = {"Tcell": ["CD3D", "CD3E", "CD2"], "Bcell": ["notagene1", "notagene2"]}
    pyucell.compute_ucell_scores(adata, signatures=signatures_miss)
    signature_columns_exist(adata, signatures_miss)


def test_layers(adata, signatures):
    adata.layers["newlayer"] = adata.X.copy()
    pyucell.compute_ucell_scores(adata, signatures=signatures, layer="newlayer")
    signature_columns_exist(adata, signatures)


def test_knn_basic(adata_with_scores, signatures):

    suffix1 = "_UCell"
    suffix2 = "_kNN"
    obs_cols = [s + suffix1 for s in signatures.keys()]
    pyucell.smooth_knn_scores(adata_with_scores, obs_columns=obs_cols, suffix=suffix2)
    signature_columns_exist(adata_with_scores, obs_cols, suffix=suffix2)

def test_knn_from_graph(adata_with_scores, signatures):

    suffix1 = "_UCell"
    suffix2 = "_kNN"
    obs_cols = [s + suffix1 for s in signatures.keys()]
    sc.pp.neighbors(adata_with_scores, n_neighbors=10, use_rep="X_pca", key_added="customgraph")
    pyucell.smooth_knn_scores(adata_with_scores, obs_columns=obs_cols, graph_key="customgraph_connectivities")
    signature_columns_exist(adata_with_scores, obs_cols, suffix=suffix2)

def test_knn_uponly(adata_with_scores, signatures):

    suffix1 = "_UCell"
    suffix2 = "_kNN"
    obs_cols = [s + suffix1 for s in signatures.keys()]
    pyucell.smooth_knn_scores(adata_with_scores, obs_columns=obs_cols, up_only=True)
    signature_columns_exist(adata_with_scores, obs_cols, suffix=suffix2)

def test_knn_invalud(adata_with_scores, signatures):

    suffix1 = "_UCell"
    suffix2 = "_kNN"
    obs_cols = [s + suffix1 for s in signatures.keys()]

    with pytest.raises(ValueError, match="decay must be between 0 and 1"):
        pyucell.smooth_knn_scores(adata_with_scores, obs_columns=obs_cols, decay=-1.0)

    with pytest.raises(ValueError, match="Graph 'not_a_graph' not found in adata.obsp"):
        pyucell.smooth_knn_scores(adata_with_scores, obs_columns=obs_cols, graph_key="not_a_graph")



