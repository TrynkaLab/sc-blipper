#!/usr/bin/env python

import argparse
import sys
import anndata as ad
from cnmf import Preprocess


if __name__ == "__main__":
    

    parser = argparse.ArgumentParser(description="Process and integrate h5ad file")

    parser.add_argument('-i', '--input', type=str, required=True, help='Input single h5ad file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output prefix')
    parser.add_argument('--harmony_vars', type=str, nargs='+', default=None, help='List of metadata variable names to integrate on (space-separated)')
    parser.add_argument('--n_variable', type=int, default=2000, help='Number of variable genes to select (default: 2000)')
    parser.add_argument('--feature_type_col', type=str, default=None, help='Column name in .var specifying if a feature is GEX or ADT')
    parser.add_argument('--seed', type=int, default=42, help='Seed for cnmf preprocess')

    args = parser.parse_args()
    
    print("Input file:", args.input)
    print("Output prefix:", args.output)
    print("Harmony variables:", args.harmony_vars)
    print("Number of variable genes:", args.n_variable)
    print("Feature type column:", args.feature_type_col)

    # Read the anndata
    adata = ad.read_h5ad(args.input)
    
    p = Preprocess(random_seed=args.seed)
    
    # Batch correct the data and save the corrected high-variance gene data to adata_c, and the TPM normalized data to adata_tpm 
    (adata_c, adata_tpm, hvgs) = p.preprocess_for_cnmf(adata,
                                                       harmony_vars=args.harmony_vars,
                                                       n_top_rna_genes=args.n_variable,
                                                       librarysize_targetsum=1e6,
                                                       feature_type_col=args.feature_type_col,
                                                        save_output_base=args.output)
