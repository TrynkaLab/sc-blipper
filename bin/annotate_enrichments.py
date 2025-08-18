#!/usr/bin/env python 
import pandas as pd
import numpy as np
import argparse
from statsmodels.stats.multitest import multipletests

def merge_with_conditions(main_tsv, output_tsv, mapping_tsv=None):
    """
    Reads the main TSV file, performs Benjamini-Hochberg correction on the pvalue column
    to add 'padj_global', then optionally merges with mapping TSV on 'condition',
    and saves the output.
    """
    # Load main TSV
    if main_tsv.endswith(".gz"):
        df_main = pd.read_csv(main_tsv, sep="\t", compression='gzip')
    else:
        df_main = pd.read_csv(main_tsv, sep="\t")

    # Check pvalue column exists
    if "pvalue" not in df_main.columns:
        raise ValueError("Input TSV must contain a 'pvalue' column for correction")
    
    if "test" not in df_main.columns:
        raise ValueError("Input TSV must contain a 'test' column for correction")
    
    # ----------------------------------------------------------------------
    # Perform Benjamini-Hochberg correction for False Discovery Rate
    # multipletests returns four things; we only need the corrected p-values here
    # Separate rows with valid and missing pvalues
    df_valid = df_main.dropna(subset=["pvalue"]).copy()
    df_na = df_main[df_main["pvalue"].isna()].copy()

    # BH correction on valid pvalues only
    _, padj, _, _ = multipletests(df_valid["pvalue"], method="fdr_bh")
    df_valid["padj_global"] = padj

    # Assign NaN padj_global to rows with missing pvalue
    df_na["padj_global"] = pd.NA
    df_na["padj_group"] = pd.NA
    df_na["padj_test"] = pd.NA

    # ----------------------------------------------------------------------
    # Perform per test FDR correction
    # Prepare a column for adjusted p-values

    if "filename" in df_main.columns:
        df_valid["padj_group"] =  df_valid["test"].astype(str) + ";" + df_valid["filename"].astype(str)
    else:
        df_valid["padj_group"] =  df_valid["test"]

    df_valid["padj_test"] = pd.NA

    # Perform BH correction within each 'test' group separately
    for test_val in df_valid['padj_group'].unique():
        mask = df_valid['padj_group'] == test_val
        pvals = df_valid.loc[mask, "pvalue"]
        _, padj, _, _ = multipletests(pvals, method="fdr_bh")
        df_valid.loc[mask, "padj_test"] = padj
    
    # ----------------------------------------------------------------------
    # Recombine and restore original order by index
    df_main = pd.concat([df_valid, df_na]).sort_index()

    # Optional merge with mapping file on 'condition'
    if mapping_tsv:
        
        df_map = pd.read_csv(mapping_tsv, sep="\t")

        df_map.rename(columns=lambda x: "annot_" + x if x not in ["condition"] else x, inplace=True)

        if "condition" not in df_main.columns:
            raise ValueError("Input TSV must contain a 'condition' column for annotating")

        if "condition" not in df_map.columns:
            raise ValueError("Annotation TSV must contain a 'condition' column for annotating")
    
        df_main['condition'] = df_main['condition'].astype(str)
        df_map['condition'] = df_map['condition'].astype(str)
        df_merged = df_main.merge(df_map, on="condition", how="left")
    else:
        df_merged = df_main

    # Save result
    df_merged.to_csv(output_tsv, sep="\t", index=False, na_rep='NA')

    return df_merged

def main():
    parser = argparse.ArgumentParser(description="Add BH corrected pvalues (padj_global) to TSV and optionally merge extra columns by 'condition'.")
    parser.add_argument("--input", help="Path to the main TSV file", required=True)
    parser.add_argument("--output", help="Path to save the output TSV", required=True)
    parser.add_argument("--annot", help="Optional mapping TSV file with extra columns to merge on 'condition'", default=None)

    args = parser.parse_args()

    merged_df = merge_with_conditions(
        main_tsv=args.input,
        output_tsv=args.output,
        mapping_tsv=args.annot
    )

    print("✅ Processing complete. Preview of output:")
    print(merged_df.head())


if __name__ == "__main__":
    main()
