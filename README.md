# pyUCell: Robust and scalable single-cell signature scoring 

[![PyPI](https://img.shields.io/pypi/v/pyucell.svg)](https://pypi.org/project/pyucell/)
[![Tests][badge-tests]][tests]
[![Documentation][badge-docs]][documentation]
[![codecov](https://codecov.io/gh/carmonalab/pyucell/graph/badge.svg?token=J4DES60KDX)](https://codecov.io/gh/carmonalab/pyucell)
[![Downloads](https://static.pepy.tech/badge/pyucell)](https://pepy.tech/project/pyucell)
[![Downloads](https://static.pepy.tech/badge/pyucell/month)](https://pepy.tech/project/pyucell)


[badge-tests]: https://img.shields.io/github/actions/workflow/status/carmonalab/pyucell/test.yaml?branch=master
[badge-docs]: https://img.shields.io/readthedocs/pyucell


In single-cell RNA-seq analysis, gene signature (or “module”) scoring constitutes a simple yet powerful approach to evaluate the strength of biological signals – typically associated to a specific cell type or biological process – in a transcriptome.

UCell is a computational method for evaluating gene signatures in single-cell datasets. UCell signature scores, based on the Mann-Whitney U statistic, are robust to dataset size and heterogeneity, and their calculation demands less computing time and memory than other available methods, enabling the processing of large datasets in a few minutes even on machines with limited computing power.

[pyUCell](https://github.com/carmonalab/pyucell) is a python implementation for the UCell algorithm, also available for the R programming language ([Bioconductor](https://bioconductor.org/packages/UCell/) and [GitHub](https://github.com/carmonalab/UCell))

### Getting started

Please see installation instructions below, and refer to the [documentation][].

### Installation

Install the latest release of `pyUCell` from [PyPI][]:

```bash
pip install pyucell
```


or, for the latest development version:

```bash
pip install git+ssh://git@github.com/carmonalab/pyucell.git@master
```


### Test the installation
```python
import pyucell as uc
import scanpy as sc

adata = sc.datasets.pbmc3k()

signatures = {
    'T_cell': ['CD3D', 'CD3E', 'CD2'],
    'B_cell': ['MS4A1', 'CD79A', 'CD79B']
}

uc.compute_ucell_scores(adata, signatures=signatures)
```

### Tutorials and how-to

Have a look at the [documentation][] section; you may start from a [basic tutorial][] or explore [some important pyUCell parameters][]

### Get help

Please address your questions and bug reports at: [UCell issues](https://github.com/carmonalab/pyucell/issues).

### Citation

UCell: robust and scalable single-cell gene signature scoring. Massimo Andreatta & Santiago J Carmona **(2021)** *CSBJ* https://doi.org/10.1016/j.csbj.2021.06.043


### Developer guide for scverse tools

https://github.com/scverse/cookiecutter-scverse?tab=readme-ov-file


[scverse discourse]: https://discourse.scverse.org/
[issue tracker]: https://github.com/carmonalab/pyucell/issues
[tests]: https://github.com/carmonalab/pyucell/actions/workflows/test.yaml
[documentation]: https://pyucell.readthedocs.io
[changelog]: https://pyucell.readthedocs.io/en/latest/changelog.html
[api documentation]: https://pyucell.readthedocs.io/en/latest/api.html
[pypi]: https://pypi.org/project/pyucell
[basic tutorial]: https://pyucell.readthedocs.io/en/latest/notebooks/basic.html
[some important pyUCell parameters]: https://pyucell.readthedocs.io/en/latest/notebooks/parameters.html
