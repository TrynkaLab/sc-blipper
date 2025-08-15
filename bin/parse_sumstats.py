#!/usr/bin/env python3
import pandas as pd
import argparse
import sys
import os

def extract_columns_streaming(input_file, col1, col2, output_file, chunksize=100000):
    compression = 'infer'  # guess from extension

    first_write = True  # whether to write header
    total_rows_written = 0

    print("🔄 Streaming reading of sumstats...")

    try:
        for chunk in pd.read_csv(input_file, 
                                 sep='\t', 
                                 compression=compression, 
                                 chunksize=chunksize,
                                 dtype=str):  # keep as string to handle "NA"
            # Check if columns exist in this chunk (only check once, early)
            if first_write:
                missing = [c for c in (col1, col2) if c not in chunk.columns]
                if missing:
                    print(f"❌ Missing column(s): {', '.join(missing)}")
                    print(f"Available columns: {', '.join(chunk.columns)}")
                    sys.exit(1)

            # Filter
            filtered = chunk[[col1, col2]].dropna(subset=[col1, col2])
            for col in [col1, col2]:
                filtered = filtered[filtered[col] != 'NA']

            # Write out (append mode after first chunk)
            filtered.to_csv(
                output_file,
                sep='\t',
                index=False,
                header=False,
                mode='w' if first_write else 'a'
            )

            total_rows_written += len(filtered)
            first_write = False

    except FileNotFoundError:
        print(f"❌ Error: File '{input_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error reading file: {e}")
        sys.exit(1)

    print(f"✅ Finished. Wrote {total_rows_written} rows to '{output_file}'.")

def main():
    parser = argparse.ArgumentParser(description="Stream-extract two columns from a TSV/TSV.GZ file.")
    parser.add_argument("--input", required=True, help="Path to the input TSV/TSV.GZ file")
    parser.add_argument("--snp-col", required=True, help="First column name to extract")
    parser.add_argument("--pval-col", required=True, help="Second column name to extract")
    parser.add_argument("--output", required=True, help="Path to the output TSV file")
    parser.add_argument("--chunksize", type=int, default=100000, help="Number of rows to process per chunk")
    args = parser.parse_args()

    extract_columns_streaming(args.input, args.snp_col, args.pval_col, args.output, args.chunksize)

if __name__ == "__main__":
    main()
