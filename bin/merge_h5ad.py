#!/usr/bin/env python3
import anndata as ad
import sys
import os
import argparse

def basename_no_ext(filepath):
    base = os.path.basename(filepath)
    while True:
        base, ext = os.path.splitext(base)
        if not ext:
            break
    return base

def merge_h5ad_files(output_file, input_files, overlap_only=False, nonoverlap_file=None, check=True):
    # Read all the input h5ad files
    adatas = [ad.read_h5ad(f) for f in input_files]

    common_genes = None
    non_overlapping_records = []

    if overlap_only:
        print("Finding overlapping genes across all datasets...")
        # Get the intersection
        common_genes = set(adatas[0].var_names)
        all_genes = set(adatas[0].var_names)
        for a in adatas[1:]:
            common_genes &= set(a.var_names)
            all_genes.update(a.var_names)

        if not common_genes:
            print("No overlapping genes found! Exiting.")
            sys.exit(1)
        
        
        perc_ol = (len(common_genes) / len(all_genes))*100
        
        if check:
            if perc_ol < 25:
                print(f"Fewer then 25% of input genes overlapped {len(common_genes)}/{len(all_genes)} ({perc_ol:.2f}%), exiting. Add --no_check to ignore this")
                sys.exit(1)
        
        common_genes = sorted(list(common_genes))
        print(f"Found {len(common_genes)} ({perc_ol:.2f}%) overlapping genes.")
        

        # Capture non-overlapping genes by dataset
        for fname, a in zip(input_files, adatas):
            unique_genes = sorted(set(a.var_names) - set(common_genes))
            if unique_genes:
                for gene in unique_genes:
                    non_overlapping_records.append((basename_no_ext(fname), gene))
                print(f"{len(unique_genes)} non-overlapping genes in {fname}")

        # Save to file if requested
        if non_overlapping_records and nonoverlap_file:
            with open(nonoverlap_file, "w") as fh:
                fh.write("file\tgene\n")
                for src, gene in non_overlapping_records:
                    fh.write(f"{src}\t{gene}\n")
            print(f"Non-overlapping genes written to {nonoverlap_file}")

        # Subset to only overlapping genes
        adatas = [a[:, common_genes] for a in adatas]

    # Append batch names as prefix to obs_names
    batch_names = []
    for batch_name, adata in zip(input_files, adatas):
        batch_name = basename_no_ext(batch_name)
        adata.obs_names = [f"{batch_name}_{obs_id}" for obs_id in adata.obs_names]
        batch_names.append(batch_name)

    # Merge AnnData
    join_type = 'outer' if not overlap_only else 'inner'
    merged = ad.concat(adatas, join=join_type, label='orig_h5ad', keys=batch_names)

    merged.write_h5ad(output_file)
    print(f"Merged {len(input_files)} files into {output_file}")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge multiple .h5ad files.")
    parser.add_argument("output", help="Output merged h5ad file")
    parser.add_argument("inputs", nargs="+", help="Input .h5ad files")
    parser.add_argument("--overlap", action="store_true",
                        help="Only keep genes present in ALL input files before merging. Otherwise they are put as NA in the count matrix")
    parser.add_argument("--nonoverlap_file", default="non_overlapping_genes.tsv",
                        help="Path to save non-overlapping genes table (TSV)")
    parser.add_argument("--no_check", default=False,
                        help="Skip the overlapping gene check")
    args = parser.parse_args()

    merge_h5ad_files(args.output,
                     args.inputs,
                     overlap_only=args.overlap,
                     nonoverlap_file=args.nonoverlap_file,
                     check=not args.no_check)
