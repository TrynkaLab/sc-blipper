#!/usr/bin/env python3
"""
Update gene names in a .h5ad file based on a mapping TSV.

Usage:
    python update_gene_ids.py input.h5ad gene_mapping.tsv output.h5ad

The TSV should have no header, and contain:
    old_id    new_id

- If an old ID is in the mapping, it gets replaced with the new ID.
- If an old ID is not mapped, the original ID is kept and a warning is shown.
- Optionally, unmapped gene IDs are written to a file (if filename provided).
"""

import sys
import pandas as pd
import anndata as ad

def update_gene_ids(h5ad_file, mapping_file, unmapped_file, output_file):
    # Load AnnData
    print(f"Loading {h5ad_file} ...")
    adata = ad.read_h5ad(h5ad_file)

    # Read mapping TSV
    print(f"Reading mapping from {mapping_file} ...")
    mapping_df = pd.read_csv(mapping_file, sep="\t", header=None, names=["old", "new"])
    mapping_dict = pd.Series(mapping_df["new"].values, index=mapping_df["old"]).to_dict()

    # Apply mapping
    updated_var_names = []
    unmapped_genes = []

    for gene in adata.var_names:
        if gene in mapping_dict:
            updated_var_names.append(mapping_dict[gene])
        else:
            updated_var_names.append(gene)  # keep original
            unmapped_genes.append(gene)

    adata.var_names = updated_var_names

    # Output warning + file
    if unmapped_genes:
        print(f"Warning: {len(unmapped_genes)} gene IDs could not be mapped and were kept as original.")
        if unmapped_file:
            pd.Series(unmapped_genes).to_csv(unmapped_file, index=False, header=False)
            print(f"Unmapped gene IDs written to: {unmapped_file}")

    # Save new file
    print(f"Saving updated AnnData to {output_file} ...")
    adata.write_h5ad(output_file)
    print("Done!")

if __name__ == "__main__":
    if len(sys.argv) < 4 or len(sys.argv) > 5:
        print("Usage: python update_gene_ids.py <input.h5ad> <mapping.tsv> <unmapped.txt> <output.h5ad>")
        sys.exit(1)

    #unmapped_file = sys.argv[4] if len(sys.argv) == 5 else None
    update_gene_ids(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
