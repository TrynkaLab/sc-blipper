#!/usr/bin/env bash
# Usage: ./txt_to_gmt.sh input.txt column_number output.gmt [gene_set_name]
# Example: ./txt_to_gmt.sh genes.txt 1 myset.gmt MyGeneSet

infile="$1"
colnum="$2"
outfile="$3"
setname="${4:-GeneSet}"   # default set name if not provided

# ------------------------------------------------------------
# Sanity checks
if [[ -z "$infile" || -z "$colnum" || -z "$outfile" ]]; then
    echo "Usage: $0 input.txt column_number output.gmt [gene_set_name]"
    exit 1
fi

if [[ ! -f "$infile" ]]; then
    echo "Error: Input file '$infile' not found."
    exit 1
fi

if ! [[ "$colnum" =~ ^[0-9]+$ ]]; then
    echo "Error: column_number must be a positive integer."
    exit 1
fi

# ------------------------------------------------------------
# Extract specified column and convert to GMT format
{
    echo -n -e "${setname}\tNA\t"
    awk -v col="$colnum" -F'\t' '{print $col}' "$infile" | paste -sd '\t' -
    echo
} > "$outfile"

echo "Saved GMT file to $outfile (column $colnum from $infile, set name '$setname')"
