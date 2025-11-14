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
