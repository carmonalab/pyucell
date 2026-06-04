# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog][],
and this project adheres to [Semantic Versioning][].

[keep a changelog]: https://keepachangelog.com/en/1.0.0/
[semantic versioning]: https://semver.org/spec/v2.0.0.html

## Version 0.3.0

### Added

	- First stable implementation of the UCell algorithm
	- Implements gene ranking and calculation of signature scores
	- Compared to the R version, we also include two different ways of
	  handling missing genes ("impute" or "skip", see the missing_genes parameter)

## Version 0.4.0

### Added

	- Smoothing of UCell scores by k-neareast neighbors. Implemented
	  in new function `smooth_knn_scores()`

## Version 0.5.0

### Added

	- Fixed a bug in `get_rankings()` where ties spanning max_rank could cause broadcasting errors.


## Version 0.7.0

### Added

	- Optional PyTorch GPU backend for `get_rankings`, `compute_ucell_scores`, and
	  `smooth_knn_scores`. Enable via a new `device` keyword argument
	  (`"auto" | "cuda" | "mps" | "cpu"`). Default behavior is unchanged.
	- New optional install extra: `pip install pyucell[gpu]` (adds `torch>=2.1`).
	- Vectorised neighbour-weight construction in `smooth_knn_scores`: the per-cell
	  Python loop is replaced by a single sparse matmul, faster even on CPU.

### Notes

	- Credit to Erick Armingol (https://orcid.org/0000-0002-1546-9165) for GPU
	  support.
	- The torch backend supports `ties_method="min"` and `"ordinal"` (for CUDA), only
	  `"ordinal"` for MPS. Passing
	  `"average"` together with `device=` raises a clear `ValueError`.
	- When `device` is set, chunks run serially (no joblib subprocesses) and
	  the default `chunk_size` is 5000 to better saturate the GPU.

## Version 0.7.3

### Added

	- Implement scoring from pre-computed rank matrix. This can be useful to test
	  new signature without recalculating ranks, but can demand more memory as it
	  requires storing the full rank matrix (unlike the regular pipeline, which
	  utilizes chunks of cells).
	- Fixed some edge cases `s_max == s_min`
	  (e.g. a 1-gene signature with `max_rank=1`) and `missing_genes="skip"` and all
	  signature genes are absent from the dataset (both pos and neg index lists empty).
	- Fixed `smooth_knn_scores` silently overwriting `adata.obsp["connectivities"]`
	  when called without `graph_key`. The internally computed graph is now stored
	  under the dedicated key `"smooth_knn_connectivities"` to avoid clobbering
	  pre-existing neighbor graphs used for UMAP/clustering.
	- Fixed self-loop double-counting in `_build_weight_matrix`: diagonal entries
	  are now stripped from the input graph before computing neighbor weights.
