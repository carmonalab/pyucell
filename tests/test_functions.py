import numpy as np
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



# ---------------------------------------------------------------------------
# compute_scores_from_ranks
# ---------------------------------------------------------------------------

def test_scores_from_ranks_matches_compute_ucell(adata, signatures):
    ranks = pyucell.get_rankings(adata)
    pyucell.compute_scores_from_ranks(adata, ranks, signatures=signatures)
    ref = adata.copy()
    pyucell.compute_ucell_scores(ref, signatures=signatures)
    for sig in signatures:
        col = f"{sig}_UCell"
        np.testing.assert_allclose(
            adata.obs[col].to_numpy(),
            ref.obs[col].to_numpy(),
            atol=1e-5,
        )


def test_scores_from_ranks_columns_exist(adata, signatures):
    ranks = pyucell.get_rankings(adata)
    pyucell.compute_scores_from_ranks(adata, ranks, signatures=signatures)
    signature_columns_exist(adata, signatures)


def test_scores_from_ranks_custom_suffix(adata, signatures):
    ranks = pyucell.get_rankings(adata)
    pyucell.compute_scores_from_ranks(adata, ranks, signatures=signatures, suffix="_score")
    signature_columns_exist(adata, signatures, suffix="_score")


def test_scores_from_ranks_reuse(adata):
    sigs_a = {"Tcell": ["CD3D", "CD3E", "CD2"]}
    sigs_b = {"Bcell": ["MS4A1", "CD79A", "CD79B"]}
    ranks = pyucell.get_rankings(adata)
    pyucell.compute_scores_from_ranks(adata, ranks, signatures=sigs_a)
    pyucell.compute_scores_from_ranks(adata, ranks, signatures=sigs_b)
    ref = adata.copy()
    pyucell.compute_ucell_scores(ref, signatures={**sigs_a, **sigs_b})
    for sig in {**sigs_a, **sigs_b}:
        col = f"{sig}_UCell"
        np.testing.assert_allclose(
            adata.obs[col].to_numpy(),
            ref.obs[col].to_numpy(),
            atol=1e-5,
        )


def test_scores_from_ranks_neg_signatures(adata):
    sigs = {"Tcell": ["CD3D+", "CD3E+", "LYZ-"]}
    ranks = pyucell.get_rankings(adata)
    pyucell.compute_scores_from_ranks(adata, ranks, signatures=sigs, w_neg=0.5)
    ref = adata.copy()
    pyucell.compute_ucell_scores(ref, signatures=sigs, w_neg=0.5)
    np.testing.assert_allclose(
        adata.obs["Tcell_UCell"].to_numpy(),
        ref.obs["Tcell_UCell"].to_numpy(),
        atol=1e-5,
    )


def test_scores_from_ranks_missing_genes_skip(adata):
    sigs = {"Tcell": ["CD3D", "CD3E", "CD2"], "Bcell": ["MS4A1", "CD79A", "notagene"]}
    ranks = pyucell.get_rankings(adata)
    pyucell.compute_scores_from_ranks(adata, ranks, signatures=sigs, missing_genes="skip")
    signature_columns_exist(adata, sigs)


def test_scores_from_ranks_all_missing_skip(adata):
    # Regression test: empty index lists with missing_genes="skip" previously
    # raised TypeError when trying to subscript a Python float.
    sigs = {"Tcell": ["CD3D", "CD3E"], "Ghost": ["notagene1", "notagene2"]}
    ranks = pyucell.get_rankings(adata)
    pyucell.compute_scores_from_ranks(adata, ranks, signatures=sigs, missing_genes="skip")
    signature_columns_exist(adata, sigs)
    assert (adata.obs["Ghost_UCell"] == 0.0).all()


def test_scores_from_ranks_shape_mismatch(adata, signatures):
    wrong_ranks = sparse.csr_matrix((adata.n_vars + 1, adata.n_obs))
    with pytest.raises(ValueError, match="ranks shape"):
        pyucell.compute_scores_from_ranks(adata, wrong_ranks, signatures=signatures)


# ---------------------------------------------------------------------------
# Optional torch backend
# ---------------------------------------------------------------------------

try:
    import torch
    HAS_TORCH = True
except ImportError:
    HAS_TORCH = False

# Create a reusable decorator for torch-dependent tests
requires_torch = pytest.mark.skipif(not HAS_TORCH, reason="torch not installed; skipping GPU backend tests")


def _torch_devices():
    if not HAS_TORCH:
        return []
    devs = ["cpu"]
    if torch.cuda.is_available():
        devs.append("cuda")
    mps = getattr(torch.backends, "mps", None)
    if mps is not None and mps.is_available():
        devs.append("mps")
    return devs


# Apply the decorator to your torch test cases
@requires_torch
@pytest.mark.parametrize("device", _torch_devices())
def test_compute_ucell_torch_matches_cpu(adata, signatures, device):
    # CPU reference uses ties='min' so it matches the torch backend's ordinal/min semantics.
    cpu_ad = adata.copy()
    pyucell.compute_ucell_scores(cpu_ad, signatures=signatures, ties_method="min")
    gpu_ad = adata.copy()
    pyucell.compute_ucell_scores(gpu_ad, signatures=signatures, ties_method="min", device=device)
    for sig in signatures:
        col = f"{sig}_UCell"
        np.testing.assert_allclose(
            gpu_ad.obs[col].to_numpy(),
            cpu_ad.obs[col].to_numpy(),
            atol=1e-5,
        )


@requires_torch
def test_compute_ucell_torch_rejects_average_ties(adata, signatures):
    with pytest.raises(ValueError, match="ties_method"):
        pyucell.compute_ucell_scores(adata, signatures=signatures, ties_method="average", device="cpu")


@requires_torch
@pytest.mark.parametrize("device", _torch_devices())
def test_smooth_knn_torch_matches_cpu(adata_with_scores, signatures, device):
    cols = [f"{s}_UCell" for s in signatures]
    cpu_ad = adata_with_scores.copy()
    pyucell.smooth_knn_scores(cpu_ad, obs_columns=cols, suffix="_kNN_cpu")
    gpu_ad = adata_with_scores.copy()
    pyucell.smooth_knn_scores(gpu_ad, obs_columns=cols, suffix="_kNN_gpu", device=device)
    for col in cols:
        np.testing.assert_allclose(
            gpu_ad.obs[f"{col}_kNN_gpu"].to_numpy(),
            cpu_ad.obs[f"{col}_kNN_cpu"].to_numpy(),
            atol=1e-5,
        )