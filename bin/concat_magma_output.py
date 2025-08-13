#!/usr/bin/env python3
import sys
import os
import re

def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} output.tsv input1.out input2.out ...")
        sys.exit(1)

    output_file = sys.argv[1]
    input_files = sys.argv[2:]

    header = None
    merged_rows = []

    for file_path in input_files:
        
        read_header=False
        with open(file_path, "r") as f:
            for line in f:
                if line.strip().startswith("#") or not line.strip():
                    continue
                parts = line.strip().split()
                
                # First non-comment line in the first file becomes header
                if not read_header:
                    read_header = True
                    if header is None:
                        header = parts + ["FILENAME", "TRAIT", "DATABASE"]
                    continue
                
                filename = os.path.basename(file_path)
                parts.append(filename)
                
                pattern = re.compile(r'^(.+?)__(.+?)\.gsa\.out$')
                match = pattern.match(filename)
                if match:
                    trait, database = match.groups()
                    parts.append(trait)
                    parts.append(database)
                else:
                    parts.append("NA")
                    parts.append("NA")
               
                merged_rows.append(parts)

    # Write to TSV
    with open(output_file, "w") as out:
        out.write("\t".join(header) + "\n")
        for row in merged_rows:
            out.write("\t".join(row) + "\n")

    print(f"Merged {len(input_files)} files into {output_file}")

if __name__ == "__main__":
    main()
