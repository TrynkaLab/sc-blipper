#!/usr/bin/env python

import anndata as ad
import sys
import os

def basename_no_ext(filepath):
    base = os.path.basename(filepath)  # Get the file name with extension(s)
    while True:
        base, ext = os.path.splitext(base)
        if not ext:
            break
    return base

def merge_h5ad_files(output_file, input_files):

    # Read all the input h5ad files into AnnData objects
    adatas = [ad.read_h5ad(f) for f in input_files]

    # Append batch names as prefix to obs_ids per batch before concatenation
    batch_names=[]
    for batch_name, adata in zip(input_files, adatas):
        batch_name = basename_no_ext(batch_name)
        adata.obs_names = [f"{batch_name}_{obs_id}" for obs_id in adata.obs_names]
        batch_names.append(batch_name)

    # Concatenate the AnnData objects along the observation axis (cells)
    merged = ad.concat(adatas, join='outer', label='orig_h5ad', keys=batch_names)

    # Write the merged AnnData object to the output h5ad file
    merged.write_h5ad(output_file)
    print(f"Merged {len(input_files)} files into {output_file}")


if __name__ == "__main__":
    # Example usage: python merge_h5ad.py merged.h5ad file1.h5ad file2.h5ad file3.h5ad
    if len(sys.argv) < 3:
        print("Usage: python merge_h5ad.py output.h5ad input1.h5ad input2.h5ad ...")
        sys.exit(1)
    output = sys.argv[1]
    inputs = sys.argv[2:]
    merge_h5ad_files(output, inputs)
