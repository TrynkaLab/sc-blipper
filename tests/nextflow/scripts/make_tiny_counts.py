#!/usr/bin/env python
"""Generate a tiny synthetic raw-count .h5ad for the cnmf_smoke nf-test."""

import sys

import anndata as ad
import numpy as np
import scipy.sparse as sp


out = sys.argv[1] if len(sys.argv) > 1 else "counts.h5ad"
rng = np.random.default_rng(42)
X = rng.poisson(
    rng.gamma(1.0, size=(200, 3)) @ rng.gamma(1.0, size=(3, 60))
).astype("float32")

a = ad.AnnData(sp.csr_matrix(X))
a.obs_names = [f"c{i}" for i in range(200)]
a.var_names = [f"g{j}" for j in range(60)]
a.write(out)
