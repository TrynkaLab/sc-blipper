#!/usr/bin/env python

import argparse
import h5py
import anndata as ad
import pandas as pd
import numpy as np

def main(var_file, usage_file, annot_file, output):
    # Read the obs of the input h5ad file
    with h5py.File(annot_file, "r") as f:
        obs = ad.experimental.read_elem(f["obs"])

    # Read the usages
    X = pd.read_csv(usage_file, sep="\t", compression="infer", index_col=0)

    # Align the index of usages with obs
    X.index = obs.index

    # Read the var file (spectra score)
    var = pd.read_csv(var_file, sep="\t", compression="infer", index_col=0)
    var.index = X.columns
    
    
    # Create AnnData object
    varm_values={}
    varm_values['spectra'] = var.to_numpy()
    
    #gene_names = var.columns.to_frame() 
    uns = {}
    uns['spectra_names'] = var.columns.to_frame()
    uns['spectra_names'].columns = ["gene"]
    adata = ad.AnnData(X=X, obs=obs, varm=varm_values, uns=uns)

    # Write output to h5ad
    adata.write_h5ad(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create AnnData object from usage, spectra score, and original h5ad files.")
    parser.add_argument("--spectra", required=True, help="Path to spectra score file (TSV, possibly compressed)")
    parser.add_argument("--usage", required=True, help="Path to usage file (TSV, possibly compressed)")
    parser.add_argument("--obs", required=True, help="Path to input annotation h5ad file")
    parser.add_argument("--output", required=True, help="Output path for resulting h5ad file")

    args = parser.parse_args()

    main(args.spectra, args.usage, args.obs, args.output)
