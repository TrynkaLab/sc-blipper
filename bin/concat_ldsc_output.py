#!/usr/bin/env python3
import sys
import os
import re
import pandas as pd

# Filename format: {pheno}__{condition}__ws{window_size}__top{binarize_top|NA}.results
PATTERN = re.compile(r'^(.+?)__(.+?)__ws(\d+)__top(\w+)\.results$')


def parse_filename(base):
    m = PATTERN.match(base)
    if m:
        pheno_id     = m.group(1)
        condition    = m.group(2)
        window_size  = int(m.group(3))
        bt_raw       = m.group(4)
        binarize_top = bt_raw if bt_raw == 'NA' else int(bt_raw)
    else:
        parts = base.replace(".results", "").split("__", 1)
        pheno_id     = parts[0]
        condition    = parts[1] if len(parts) > 1 else "unknown"
        window_size  = None
        binarize_top = None
    return pheno_id, condition, window_size, binarize_top


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <file1.results> [file2.results ...]")
        sys.exit(1)

    input_files = sys.argv[1:]

    records = []
    for fpath in sorted(input_files):
        base = os.path.basename(fpath)
        pheno_id, condition, window_size, binarize_top = parse_filename(base)
        df = pd.read_csv(fpath, sep="\t")
        df.insert(0, "binarize_top", binarize_top)
        df.insert(0, "window_size", window_size)
        df.insert(0, "condition", condition)
        df.insert(0, "phenotype", pheno_id)
        records.append(df)

    if not records:
        raise RuntimeError("No .results files provided")

    combined = pd.concat(records, ignore_index=True)
    combined.to_csv("ldsc_results_aggregated.tsv", sep="\t", index=False)
    print(f"Aggregated {len(records)} result files into ldsc_results_aggregated.tsv")

    annot_indices = combined["Category"].str.extract(r'_(\d+)$', expand=False).dropna().unique()
    for idx in sorted(annot_indices, key=int):
        subset = combined[combined["Category"].str.endswith(f"_{idx}")]
        out = f"ldsc_results_aggregated_annot{idx}.tsv"
        subset.to_csv(out, sep="\t", index=False)
        print(f"Wrote {len(subset)} rows to {out}")


if __name__ == "__main__":
    main()
