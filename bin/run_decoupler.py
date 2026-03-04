#!/usr/bin/env python3
"""Translate of run_decoupler.r into Python using decoupler.

Inputs/outputs mirror the R version:
- <output_prefix>_progeny_activities.tsv
- <output_prefix>_progeny_activity_matrix.tsv
- <output_prefix>_pathways.pdf
- <output_prefix>_pathways_scaled.pdf
- <output_prefix>_collectri_activities.tsv
- <output_prefix>_collectri_activity_matrix.tsv
- <output_prefix>_signif_tfs.pdf
- <output_prefix>_signif_tfs_scaled.pdf

This script expects a CSV with rownames in the first column (genes in rows).
"""
import argparse
import os
import sys
import math
import pandas as pd
import numpy as np
import decoupler as dc
import matplotlib.pyplot as plt
import seaborn as sns


def read_matrix(path):
    # Try reading with flexible separator detection first, fallbacks to common formats
    df = pd.DataFrame()
    readers = [
        {"sep": None, "engine": "python"},
        {"sep": "\t"},
        {"delim_whitespace": True},
    ]
    for kw in readers:
        try:
            df = pd.read_csv(path, header=0, index_col=0, compression="infer", **kw)
            if df.shape[0] > 0 and df.shape[1] > 0:
                break
        except Exception:
            df = pd.DataFrame()

    if df.empty or df.shape[1] == 0:
        # Try reading without header (some files may lack header row)
        try:
            df2 = pd.read_csv(path, header=None, compression="infer", sep=None, engine="python")
            if df2.shape[1] >= 2:
                df2.columns = range(df2.shape[1])
                df = df2.set_index(0)
            else:
                raise ValueError("Input matrix appears to have a single column or is malformed")
        except Exception as e:
            raise ValueError(f"Could not read matrix from {path}: {e}")

    # Coerce all data columns to numeric (like R fread + as.matrix)
    df = df.apply(pd.to_numeric, errors="coerce")
    # Drop rows with any NA to emulate R's na.omit(mat)
    df = df.dropna()
    if df.empty:
        raise ValueError(f"Matrix from {path} contains no numeric data after coercion/dropna.")
    return df

def maybe_map_targets(net_df, mapping):
    if mapping is None:
        return net_df
    if "target" in net_df.columns:
        net_df["target"] = net_df["target"].map(lambda x: mapping.get(x, x))
    return net_df


def read_network_tsv(path):
    return pd.read_csv(path, sep="\t")


def melt_acts(acts_df, pvals_df=None):
    # acts_df: index = sources, columns = conditions
    acts_long = acts_df.reset_index().melt(id_vars=acts_df.index.name or "source", var_name="condition", value_name="score")
    # ensure consistent column name for source
    acts_long = acts_long.rename(columns={acts_df.index.name or "source": "source"})
    if pvals_df is not None:
        p_long = pvals_df.reset_index().melt(id_vars=pvals_df.index.name or "source", var_name="condition", value_name="p_value")
        p_long = p_long.rename(columns={pvals_df.index.name or "source": "source"})
        merged = pd.merge(acts_long, p_long, on=["source", "condition"], how="left")
        return merged
    return acts_long


def plot_heatmap(df, outfile, scale_rows=False):
    # df: rows=sources, cols=conditions
    if scale_rows:
        df_plot = df.apply(lambda r: (r - r.mean()) / (r.std() if r.std() != 0 else 1), axis=1)
    else:
        df_plot = df.copy()

    nrows, ncols = df_plot.shape
    figsize = (2 + (ncols * 0.5), 1 + (nrows * 0.5))
    plt.figure(figsize=figsize)
    sns.heatmap(df_plot, cmap="RdBu_r", center=0, cbar_kws={"shrink": 0.5})
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()

def main(args):

    #if args.cache_dir:
        # try to set omnipath cache env var used by omnipathpy/decoupler
    #    os.environ["OMNIPATH_CACHE_DIR"] = args.cache_dir

    gene_replacement = None
    gene_replacement_rev = None
    if args.id_linker:
        map_df = pd.read_csv(args.id_linker, header=None, sep="\t", names=["old", "new"], dtype=str)
        gene_replacement = dict(zip(map_df["old"], map_df["new"]))
        gene_replacement_rev = dict(zip(map_df["new"], map_df["old"]))

    mat = read_matrix(args.matrix)

    if not args.transpose:
        mat = mat.T

    print("Matrix loaded. Now running PROGENy...")
    # Load PROGENy
    if args.progeny_network_file:
        progeny = read_network_tsv(args.progeny_network_file)
    else:
        try:
            progeny = dc.op.progeny(organism="human", top=500)
        except TypeError:
            # some versions may not support top
            progeny = dc.op.progeny(organism="human")

    progeny = maybe_map_targets(progeny, gene_replacement)

    print("PROGENy network loaded. Now running pathway scoring...")
    # run pathway scoring; in R they used run_mlm for progeny; prefer mlm when available
    acts_pw, p_pw = dc.mt.mlm(data=mat, net=progeny)

    acts_pw.index.name = "condition"
    p_pw.index.name = "condition"
    acts_pw.columns.name = "source"
    p_pw.columns.name = "source"

    # Write outputs
    # activity matrix as wide (rows = source)
    pw_matrix = acts_pw.copy()
    pw_matrix.to_csv(f"{args.output_prefix}_progeny_activity_matrix.tsv", sep="\t")

    # melt long activities and include p-values
    acts_pw_long = melt_acts(acts_pw.T, p_pw.T)
    acts_pw_long['statistic'] = 'mlm'
    acts_pw_long = acts_pw_long[['statistic','source', 'condition', 'score', 'p_value']]
    acts_pw_long.to_csv(f"{args.output_prefix}_progeny_activities.tsv", sep="\t", index=False)

    print("PROGENy analysis completed. Now plotting...")
    # Visualization
    df_plot = pw_matrix.copy()
    if gene_replacement_rev:
        df_plot.index = [gene_replacement_rev.get(x, x) for x in df_plot.index]

    plot_heatmap(df_plot, f"{args.output_prefix}_pathways.pdf", scale_rows=False)
    plot_heatmap(df_plot, f"{args.output_prefix}_pathways_scaled.pdf", scale_rows=True)
  
    print("PROGENy analysis completed. Now running DoRothEA/CollectTRI...")
    # CollectTRI network
    if args.collectri_network_file:
        collectri = read_network_tsv(args.collectri_network_file)
    else:
        try:
            collectri = dc.op.collectri(organism="human")
        except Exception:
            # fallback via op.resource
            collectri = dc.op.collectri(organism="human")

    collectri = maybe_map_targets(collectri, gene_replacement)

    print("CollectTRI network loaded. Now running TF activity scoring...")
    acts_tf, p_tf = dc.mt.ulm(data=mat, net=collectri)

    acts_tf.index.name = "condition"
    p_tf.index.name = "condition"
    acts_tf.columns.name = "source"
    p_tf.columns.name = "source"

    # Write outputs
    # activity matrix as wide (rows = source)
    pw_matrix = acts_tf.copy()
    pw_matrix.to_csv(f"{args.output_prefix}_collectri_activity_matrix.tsv", sep="\t")

    # melt long activities and include p-values
    acts_tf_long = melt_acts(acts_tf.T, p_tf.T)
    acts_tf_long['statistic'] = 'ulm'
    acts_tf_long = acts_tf_long[['statistic','source', 'condition', 'score', 'p_value']]
    acts_tf_long.to_csv(f"{args.output_prefix}_collectri_activities.tsv", sep="\t", index=False)

    tf_matrix = acts_tf.copy()
    tf_matrix.index.name = "source"
    tf_matrix.to_csv(f"{args.output_prefix}_collectri_activity_matrix.tsv", sep="\t")

    # Select significant TFs similar to R code
    # compute activity matrix p (p-values) and bonferroni threshold
    try:
        activity_pvals = p_tf.copy()
    except Exception:
        activity_pvals = None

    if activity_pvals is not None:
        # Bonferroni threshold using number of rows in long activities (sources*conditions)
        n_acts = acts_tf_long.shape[0]
        thresh = 0.05 / max(1, n_acts)
        sig_mask_row = (activity_pvals < thresh).any(axis=1)
        sig_mask_col = (activity_pvals < thresh).any(axis=0)        
        df_plot2 = tf_matrix.copy()
        df_plot2 = df_plot2.loc[sig_mask_row, sig_mask_col]

    else:
        df_plot2 = tf_matrix

    if gene_replacement_rev is not None:
        df_plot2.index = [gene_replacement_rev.get(x, x) for x in df_plot2.index]

    plot_heatmap(df_plot2, f"{args.output_prefix}_signif_tfs.pdf", scale_rows=False)
    plot_heatmap(df_plot2, f"{args.output_prefix}_signif_tfs_scaled.pdf", scale_rows=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--matrix", required=True, help="Path to numeric matrix CSV file (genes in rows by default).")
    parser.add_argument("-o", "--output_prefix", required=True, help="Prefix for output files.")
    parser.add_argument("-t", "--transpose", action="store_true", help="Whether to transpose the matrix before processing. Default: FALSE.")
    #parser.add_argument("--cache_dir", default=None, help="Cache dir for OmniPath (best effort).")
    parser.add_argument("--id_linker", default=None, help="TSV file with gene name to ensembl (old, new) to convert omnipath with Default: NULL")
    parser.add_argument("--progeny_network_file", default=None, help="Optional TSV for PROGENy/decoupler network.")
    parser.add_argument("--collectri_network_file", default=None, help="Optional TSV for CollectTRI network.")

    args = parser.parse_args()
    
    main(args)