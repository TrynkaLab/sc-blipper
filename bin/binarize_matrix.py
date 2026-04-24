#!/usr/bin/env python3
"""
Binarize a numeric genes x conditions matrix based on top-N rank per column.

For each column (condition), genes are ranked and the top N are set to 1;
all others are set to 0.  Ranking is done independently per column.

Options
-------
--absolute   Take absolute values before ranking, so the top N genes by
             magnitude are selected regardless of sign.
--ascending  Select genes with the *lowest* scores (rank 1 = smallest).
             Default (flag absent) selects genes with the *highest* scores
             (rank 1 = largest).
"""
import argparse
import sys
import pandas as pd


def read_table(filename):
    compression = "gzip" if filename.endswith(".gz") else None
    base = filename[:-3] if filename.endswith(".gz") else filename
    sep = "," if base.endswith(".csv") else "\t"
    return pd.read_csv(filename, sep=sep, compression=compression, index_col=0).astype(float)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input",     required=True, help="Input matrix (TSV/CSV, optionally .gz, genes as rows)")
    parser.add_argument("--output",    required=True, help="Output binarized matrix (TSV)")
    parser.add_argument("--top",       required=True, type=int, help="Number of top genes to select per condition")
    parser.add_argument("--absolute",  action="store_true", help="Take absolute values before ranking")
    parser.add_argument("--ascending", action="store_true", help="Select lowest-scoring genes (default: highest)")
    args = parser.parse_args()

    log_path = args.output.rsplit(".", 1)[0] + ".log"
    lf = open(log_path, "w")

    def log(msg):
        print(msg, file=lf, flush=True)
        print(msg, flush=True)

    df = read_table(args.input)
    log(f"Loaded matrix: {df.shape[0]} genes x {df.shape[1]} conditions")

    if args.absolute:
        df = df.abs()
        log("Applied absolute values before ranking")

    direction = "ascending (lowest scores)" if args.ascending else "descending (highest scores)"
    log(f"Selecting top {args.top} genes per condition, ranking {direction}")

    result = pd.DataFrame(0, index=df.index, columns=df.columns, dtype=int)
    for col in df.columns:
        ranked = df[col].rank(method="first", ascending=args.ascending, na_option="bottom")
        result.loc[ranked <= args.top, col] = 1
        log(f"  {col}: {int(result[col].sum())} genes selected")

    result.to_csv(args.output, sep="\t", index=True, index_label="rowid")
    log(f"Written binarized matrix ({result.shape[0]} genes x {result.shape[1]} conditions) to {args.output}")
    lf.close()


if __name__ == "__main__":
    main()
