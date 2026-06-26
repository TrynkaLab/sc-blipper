#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Permute each column of a matrix independently N times, outputting a single combined matrix. "
                    "Output columns are named {col}_perm_{i} for condition col and permutation index i."
    )
    parser.add_argument("--input", required=True, help="Input matrix (TSV/CSV, rows=genes, cols=conditions)")
    parser.add_argument("--output", required=True, help="Output permuted matrix (TSV)")
    parser.add_argument("--n-permutations", type=int, required=True, help="Total number of permutations per column")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    parser.add_argument("--start-perm", type=int, default=1, help="First permutation index to generate (1-based, inclusive)")
    parser.add_argument("--end-perm", type=int, default=None, help="Last permutation index to generate (1-based, inclusive; defaults to --n-permutations)")
    args = parser.parse_args()

    start = args.start_perm
    end = args.end_perm if args.end_perm is not None else args.n_permutations

    # Offset seed by block start so each block is independently reproducible
    rng = np.random.default_rng(args.seed + (start - 1))

    fname = args.input.replace(".gz", "")
    sep = "," if fname.endswith(".csv") else "\t"
    df = pd.read_csv(args.input, sep=sep, index_col=0)
    print(f"Loaded {df.shape[0]} x {df.shape[1]} matrix ({df.shape[1]} conditions)")

    result = {}
    for col in df.columns:
        col_safe = col.replace(" ", "_")
        vals = df[col].values
        for i in range(start, end + 1):
            result[f"{col_safe}_perm_{i}"] = rng.permutation(vals)

    out = pd.DataFrame(result, index=df.index)
    out.to_csv(args.output, sep="\t", index=True, index_label="rowid")
    print(f"Written {out.shape[1]} permuted columns ({df.shape[1]} conditions x {end - start + 1} permutations, indices {start}-{end})")


if __name__ == "__main__":
    main()
