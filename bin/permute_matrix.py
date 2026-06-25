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
    parser.add_argument("--n-permutations", type=int, required=True, help="Number of permutations per column")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)

    fname = args.input.replace(".gz", "")
    sep = "," if fname.endswith(".csv") else "\t"
    df = pd.read_csv(args.input, sep=sep, index_col=0)
    print(f"Loaded {df.shape[0]} x {df.shape[1]} matrix ({df.shape[1]} conditions)")

    result = {}
    for col in df.columns:
        col_safe = col.replace(" ", "_")
        vals = df[col].values
        for i in range(1, args.n_permutations + 1):
            result[f"{col_safe}_perm_{i}"] = rng.permutation(vals)

    out = pd.DataFrame(result, index=df.index)
    out.to_csv(args.output, sep="\t", index=True, index_label="rowid")
    print(f"Written {out.shape[1]} permuted columns ({df.shape[1]} conditions x {args.n_permutations} permutations)")


if __name__ == "__main__":
    main()
