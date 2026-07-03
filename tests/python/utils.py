"""Reusable helpers and fixtures for Python tests."""

from pathlib import Path
from types import SimpleNamespace
import importlib.util
import os
import sys

import numpy as np
import pytest


# torch + sklearn/scipy can double-load OpenMP on macOS; tolerate it for this test process.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

KERNEL_PATH = Path(__file__).resolve().parents[2] / "bin" / "nmf_gpu.py"


def load_kernel_module(module_name="nmf_gpu", kernel_path=KERNEL_PATH):
    """Load the standalone kernel script as an importable module for tests."""
    kernel_path = Path(kernel_path)
    if not kernel_path.exists():
        pytest.fail(f"Required NMF GPU kernel file is missing: {kernel_path}", pytrace=False)

    if module_name in sys.modules:
        return sys.modules[module_name]

    spec = importlib.util.spec_from_file_location(module_name, kernel_path)
    if spec is None or spec.loader is None:
        pytest.fail(f"Unable to load {module_name!r} from {kernel_path}", pytrace=False)

    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="session")
def kernel():
    """Loaded `bin/nmf_gpu.py` module under test."""
    return load_kernel_module()


def require_nmf_runtime():
    """Skip tests that need the lazy torch + sklearn runtime dependencies."""
    torch = pytest.importorskip("torch")
    pytest.importorskip("sklearn")
    return torch


def low_rank_matrix(seed=0, cells=40, genes=20, rank=3):
    """Build an exact non-negative rank-k matrix for reconstruction tests."""
    rng = np.random.default_rng(seed)
    W = rng.random((cells, rank)) + 0.2
    H = rng.random((rank, genes)) + 0.2
    return W @ H


def small_nonnegative_matrix(seed=1, cells=10, genes=8):
    """Build a small dense non-negative matrix for quick smoke tests."""
    rng = np.random.default_rng(seed)
    return rng.random((cells, genes)) + 0.1


def assert_valid_nmf_output(X, H, W, k):
    """Assert the shared cNMF-compatible output contract for NMF factors."""
    assert H.shape == (k, X.shape[1])
    assert W.shape == (X.shape[0], k)
    assert H.dtype == np.float64
    assert W.dtype == np.float64
    assert np.isfinite(H).all()
    assert np.isfinite(W).all()
    assert (H >= 0).all()
    assert (W >= 0).all()


def fake_torch_backend(cuda_available=False, mps_available=False, bf16_supported=True):
    """Create a tiny fake torch surface for backend policy tests."""
    return SimpleNamespace(
        float32="float32",
        float64="float64",
        bfloat16="bfloat16",
        cuda=SimpleNamespace(
            is_available=lambda: cuda_available,
            is_bf16_supported=lambda: bf16_supported,
        ),
        backends=SimpleNamespace(
            mps=SimpleNamespace(is_available=lambda: mps_available),
        ),
    )
