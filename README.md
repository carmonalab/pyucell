# pyUCell: gene signature scoring for scRNA-seq 

[![Tests][badge-tests]][tests]
[![Documentation][badge-docs]][documentation]

[badge-tests]: https://img.shields.io/github/actions/workflow/status/carmonalab/pyucell/test.yaml?branch=master
[badge-docs]: https://img.shields.io/readthedocs/pyucell

⚠️ **Under Development** ⚠️

In single-cell RNA-seq analysis, gene signature (or “module”) scoring constitutes a simple yet powerful approach to evaluate the strength of biological signals – typically associated to a specific cell type or biological process – in a transcriptome.

UCell is a computational method for evaluating gene signatures in single-cell datasets. UCell signature scores, based on the Mann-Whitney U statistic, are robust to dataset size and heterogeneity, and their calculation demands less computing time and memory than other available methods, enabling the processing of large datasets in a few minutes even on machines with limited computing power.

pyUCell is the python implementation for the UCell algorithm, also available for R - [Bioconductor](https://bioconductor.org/packages/UCell/) and [GitHub](https://github.com/carmonalab/UCell)

## Getting started

Please see installation instructions below, and refer to the [documentation][].

## Installation


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

## Tutorials and how-to

Have a look at the [documentation][] section; you may start from a [basic tutorial][]

### Get help

Please address your questions and bug reports at: [UCell issues](https://github.com/carmonalab/pyucell/issues).

### Citation

UCell: robust and scalable single-cell gene signature scoring. Massimo Andreatta & Santiago J Carmona **(2021)** *CSBJ* https://doi.org/10.1016/j.csbj.2021.06.043


## Developer guide for scverse tools

https://github.com/scverse/cookiecutter-scverse?tab=readme-ov-file


[scverse discourse]: https://discourse.scverse.org/
[issue tracker]: https://github.com/carmonalab/pyucell/issues
[tests]: https://github.com/carmonalab/pyucell/actions/workflows/test.yaml
[documentation]: https://pyucell.readthedocs.io
[changelog]: https://pyucell.readthedocs.io/en/latest/changelog.html
[api documentation]: https://pyucell.readthedocs.io/en/latest/api.html
[pypi]: https://pypi.org/project/pyucell
[basic tutorial]: https://pyucell.readthedocs.io/en/latest/notebooks/basic.html
