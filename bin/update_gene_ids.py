#!/usr/bin/env python3
"""
Update gene names in a .h5ad file (optionally using a mapping TSV),
optionally filter genes, and save the result.

Usage:
    python update_gene_ids.py \
        --input input.h5ad \
        --output output.h5ad \
        [--mapping gene_mapping.tsv] \
        [--unmapped unmapped.txt] \
        [--filter-genes keep_genes.txt] \
        [--filter-mode include|exclude]

The mapping TSV should have no header and contain:
    old_id    new_id
"""

import argparse
import pandas as pd
import anndata as ad


def update_gene_ids(adata, mapping_file, unmapped_file=None):
    """Update gene IDs (var_names) in AnnData using a mapping file."""
    print(f"Reading mapping from {mapping_file} ...")
    mapping_df = pd.read_csv(mapping_file, sep="\t", header=None, names=["old", "new"])
    mapping_dict = pd.Series(mapping_df["new"].values,
                             index=mapping_df["old"]).to_dict()

    updated_var_names = []
    unmapped_genes = []

    for gene in adata.var_names:
        if gene in mapping_dict:
            updated_var_names.append(mapping_dict[gene])
        else:
            updated_var_names.append(gene)
            unmapped_genes.append(gene)

    adata.var_names = updated_var_names

    if unmapped_genes:
        print(f"Warning: {len(unmapped_genes)} gene IDs could not be mapped and were kept as original.")
        if unmapped_file:
            pd.Series(unmapped_genes).to_csv(unmapped_file, index=False, header=False)
            print(f"Unmapped gene IDs written to: {unmapped_file}")

    return adata, mapping_dict


def filter_genes(adata, filter_file, mapping_dict=None, mode="include"):
    """Filter AnnData object by genes; mode = 'include' or 'exclude'."""
    print(f"Filtering genes using {filter_file} (mode={mode}) ...")
    filter_genes = pd.read_csv(filter_file, header=None, squeeze=True).astype(str).tolist()

    # If mapping provided, update gene IDs in filter list too
    if mapping_dict is not None:
        filter_genes = [mapping_dict.get(g, g) for g in filter_genes]

    if mode == "include":
        keep_mask = adata.var_names.isin(filter_genes)
    elif mode == "exclude":
        keep_mask = ~adata.var_names.isin(filter_genes)
    else:
        raise ValueError("filter-mode must be 'include' or 'exclude'")

    print(f"Keeping {keep_mask.sum()} genes out of {adata.n_vars}")
    adata = adata[:, keep_mask]
    return adata


def main():
    parser = argparse.ArgumentParser(
        description="Update and/or filter gene IDs in an AnnData (.h5ad) file."
    )
    parser.add_argument("-i", "--input", required=True, help="Input AnnData file (.h5ad)")
    parser.add_argument("-o", "--output", required=True, help="Output AnnData file (.h5ad)")
    parser.add_argument("-m", "--mapping", default=None, help="Optional TSV with old_id and new_id (no header)")
    parser.add_argument("-u", "--unmapped", default=None, help="Optional file to write unmapped gene IDs (only if mapping used)")
    parser.add_argument("-f", "--filter-genes", default=None, help="Optional file with list of gene IDs to filter")
    parser.add_argument("--filter-mode", choices=["include", "exclude"], default="include",
                        help="Filter mode: include only listed genes, or exclude them [default: include]")

    args = parser.parse_args()

    # Load AnnData
    print(f"Loading {args.input} ...")
    adata = ad.read_h5ad(args.input)

    mapping_dict = None
    if args.mapping:
        adata, mapping_dict = update_gene_ids(adata, args.mapping, args.unmapped)
    else:
        print("No mapping file provided. Gene IDs will be kept unchanged.")

    # Optional filtering
    if args.filter_genes:
        adata = filter_genes(adata, args.filter_genes, mapping_dict, mode=args.filter_mode)

    # Save output
    print(f"Saving AnnData to {args.output} ...")
    adata.write_h5ad(args.output)
    print("Done!")


if __name__ == "__main__":
    main()
