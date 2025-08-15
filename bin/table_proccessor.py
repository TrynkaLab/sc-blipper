#!/usr/bin/env python3
import argparse
import gzip
import os
import pandas as pd

def read_table(filename):
    """Read TSV or CSV, supports .gz files."""
    compression = 'gzip' if filename.endswith('.gz') else None
    base_name = filename[:-3] if filename.endswith('.gz') else filename
    if base_name.endswith('.tsv'):
        sep = '\t'
    elif base_name.endswith('.csv'):
        sep = ','
    else:
        sep = '\t'  # default
    df = pd.read_csv(filename, sep=sep, compression=compression, dtype=str, index_col=0)
    return df, sep

def read_list(file_path):
    """Read a list of items from a file (one per line)."""
    with open(file_path) as f:
        return [line.strip() for line in f if line.strip()]

def update_index_from_mapping(df, mapping_file):
    """Update df index based on old->new ID mapping in mapping_file."""
    mapping_df = pd.read_csv(mapping_file, sep="\t", dtype=str, index_col=None)
    mapping_df.columns = ['old_id', 'new_id']
    id_map = dict(zip(mapping_df['old_id'], mapping_df['new_id']))
    df.index = df.index.to_series().replace(id_map)
    return df

def update_columns_from_mapping(df, mapping_file):
    """Update column names based on old->new ID mapping in mapping_file."""
    mapping_df = pd.read_csv(mapping_file, sep="\t", dtype=str, index_col=None)
    mapping_df.columns = ['old_id', 'new_id']
    id_map = dict(zip(mapping_df['old_id'], mapping_df['new_id']))
    df.rename(columns=id_map, inplace=True)
    return df

def main():
    parser = argparse.ArgumentParser(description="Process TSV/CSV files with gzip support, transposition, subsetting, updating row/column IDs, and format conversion.")
    parser.add_argument("--input", help="Input file (TSV, CSV, optionally .gz)")
    parser.add_argument("--output", help="Output file (TSV or CSV)")
    parser.add_argument("--output.log", help="Log file for command printouts", default=None)
    parser.add_argument("--transpose", "-t", action="store_true", help="Transpose the table")
    parser.add_argument("--out-sep", choices=["tsv", "csv"], default="tsv", help="Output format: tsv (default) or csv")
    parser.add_argument("--row-file", help="File containing row IDs to keep (matches first column). Runs after --update-rows")
    parser.add_argument("--col-file", help="File containing column names to keep (matches header). Runs after --update-cols")
    parser.add_argument("--update-rows", help="Two-column TSV/CSV file (old_id,new_id) to update DataFrame row index. Runs before transpose")
    parser.add_argument("--update-cols", help="Two-column TSV/CSV file (old_id,new_id) to update DataFrame column names. Runs before transpose")

    args = parser.parse_args()

    # Function to print to log
    def log_print(*msg):
        if args.output:
            
            logfile=args.output.replace(".tsv", "").replace(".csv", "")
            with open(f"{logfile}.log", "a") as log_fh:
                print(*msg, file=log_fh)
        else:
            print(*msg)

    # Read table
    df, in_sep = read_table(args.input)
    log_print(f"Loaded matrix of shape {df.shape}")

    # Update row index if requested
    if args.update_rows:
        df = update_index_from_mapping(df, args.update_rows)
    
    # Update column names if requested
    if args.update_cols:
        df = update_columns_from_mapping(df, args.update_cols)

    # Subset rows
    if args.row_file:
        row_ids = set(read_list(args.row_file))
        log_print(f"Loaded {len(row_ids)} row_ids in set.")
        df = df.loc[df.index.isin(row_ids)]
    
    # Subset columns
    if args.col_file:
        col_names = set(read_list(args.col_file))
        log_print(f"Loaded {len(col_names)} col_names in set.")
        keep_cols = [df.columns[0]] + [c for c in df.columns if c in col_names]
        log_print(f"{len(keep_cols)} overlap")
        df = df.loc[:, keep_cols]

    log_print(f"After subsetting, retaining matrix of shape {df.shape}")
    
    # Optionally transpose
    if args.transpose:
        df = df.T
        log_print(f"Transposed to {df.shape}")

    # Output separator
    out_sep = '\t' if args.out_sep == "tsv" else ','

    # Write file
    df.to_csv(args.output, sep=out_sep, index=True, index_label="rowid")
    log_print(f"Written output to {args.output} ({args.out_sep.upper()})")

if __name__ == "__main__":
    main()
