# pyUCell

[![Tests][badge-tests]][tests]
[![Documentation][badge-docs]][documentation]

[badge-tests]: https://img.shields.io/github/actions/workflow/status/carmonalab/pyucell/test.yaml?branch=master
[badge-docs]: https://img.shields.io/readthedocs/pyucell

⚠️ **Under Development** ⚠️

Python implementation for the UCell algorithm

## Getting started

Please refer to the [documentation][],
in particular, the [API documentation][].

## Installation

You need to have Python 3.10 or newer installed on your system.
If you don't have Python installed, we recommend installing [uv][].


<!--
1) Install the latest release of `pyUCell` from [PyPI][]:

```bash
pip install pyUCell
```
-->

1. Install the latest development version:

```bash
pip install git+ssh://git@github.com/carmonalab/pyucell.git@master
```


#### Test the installation
```python
import pyucell
import scanpy as sc

adata = sc.datasets.pbmc3k()

signatures = {
    'T_cell': ['CD3D', 'CD3E', 'CD2'],
    'B_cell': ['MS4A1', 'CD79A', 'CD79B']
}

pyucell.compute_ucell_scores(adata, signatures=signatures, chunk_size=500)
```

3. Visualize results e.g. on UMAP

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack', n_comps=50)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

sc.pl.umap(adata, color='T_cell', cmap='viridis', ax=axes[0], size=20, show=False)
sc.pl.umap(adata, color='B_cell', cmap='viridis', ax=axes[1], size=20, show=False)
```

## Developer guide for scverse tools

https://github.com/scverse/cookiecutter-scverse?tab=readme-ov-file



[uv]: https://github.com/astral-sh/uv
[scverse discourse]: https://discourse.scverse.org/
[issue tracker]: https://github.com/mass-a/pyUCell/issues
[tests]: https://github.com/mass-a/pyUCell/actions/workflows/test.yaml
[documentation]: https://pyUCell.readthedocs.io
[changelog]: https://pyUCell.readthedocs.io/en/latest/changelog.html
[api documentation]: https://pyUCell.readthedocs.io/en/latest/api.html
[pypi]: https://pypi.org/project/pyUCell
