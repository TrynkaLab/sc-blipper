#!/usr/bin/env python3
import argparse


def main():
    parser = argparse.ArgumentParser(
        description="Merge per-block MAGMA .gsa.out files into a single combined file. "
                    "Comment lines (#) and the column header are taken from the first file only; "
                    "data rows are appended from all files."
    )
    parser.add_argument("--inputs", nargs="+", required=True, help="Block .gsa.out files to merge")
    parser.add_argument("--output", required=True, help="Output merged .gsa.out file")
    args = parser.parse_args()

    # Sort by filename so block order is deterministic regardless of channel ordering
    block_files = sorted(args.inputs)

    with open(args.output, "w") as out:
        header_written = False
        total_rows = 0
        for path in block_files:
            with open(path) as f:
                found_col_header = False
                for line in f:
                    stripped = line.strip()
                    if not stripped:
                        continue
                    if stripped.startswith("#"):
                        # Write comment/metadata lines only from the first block
                        if not header_written:
                            out.write(line)
                        continue
                    if not found_col_header:
                        # First non-comment line is the column header
                        found_col_header = True
                        if not header_written:
                            out.write(line)
                            header_written = True
                        continue
                    # Data row
                    out.write(line)
                    total_rows += 1

    print(f"Merged {len(block_files)} block file(s) into {args.output} ({total_rows} data rows)")


if __name__ == "__main__":
    main()
