"""Lazy torch import and device helpers for the optional GPU backend.

Importing this module never imports torch. torch is only resolved when one of
the helpers is actually called, so users without the ``gpu`` extra installed
incur no overhead.

Credit to Erick Armingol (https://orcid.org/0000-0002-1546-9165) for GPU support.
"""

from __future__ import annotations

import numpy as np
from scipy import sparse

_TORCH = None


def import_torch():
    """Import torch lazily, raising a clear error if it isn't installed."""
    global _TORCH
    if _TORCH is not None:
        return _TORCH
    try:
        import torch
    except ImportError as e:
        raise ImportError(
            "pyucell's GPU backend requires PyTorch. "
            "Install it with `pip install pyucell[gpu]` or `pip install torch`."
        ) from e
    _TORCH = torch
    return torch


def resolve_device(device):
    """Map a device string ('auto', 'cuda', 'mps', 'cpu', ...) to a torch.device."""
    torch = import_torch()
    if device is None:
        return None
    if isinstance(device, torch.device):
        return device
    if device == "auto":
        if torch.cuda.is_available():
            return torch.device("cuda")
        mps = getattr(torch.backends, "mps", None)
        if mps is not None and mps.is_available():
            return torch.device("mps")
        return torch.device("cpu")
    dev = torch.device(device)
    if dev.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("device='cuda' requested but torch.cuda.is_available() is False.")
    if dev.type == "mps":
        mps = getattr(torch.backends, "mps", None)
        if mps is None or not mps.is_available():
            raise RuntimeError("device='mps' requested but torch.backends.mps.is_available() is False.")
    return dev


def to_torch_dense(X, device, dtype=None):
    """Convert a 2D numpy/scipy.sparse matrix into a dense torch tensor on ``device``."""
    torch = import_torch()
    if dtype is None:
        dtype = torch.float32
    if sparse.issparse(X):
        X = X.toarray()
    elif hasattr(X, "toarray"):
        X = X.toarray()
    else:
        X = np.asarray(X)
    return torch.as_tensor(np.ascontiguousarray(X), dtype=dtype, device=device)
